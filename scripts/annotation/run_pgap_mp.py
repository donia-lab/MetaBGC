#!/usr/bin/env python3
"""
Multiprocessing-based PGAP annotation pipeline.

Reads a multi-sequence FASTA file and runs PGAP on each sequence in parallel
using a pool of worker processes.

Usage:
    python run_pgap_mp.py --fasta <fasta> --outdir <outdir> \
        --submol <submol.yaml> --pgap-exe <pgap.py> --container <pgap.sif> \
        --workers 4

Requirements:
    biopython
"""

import argparse
import io
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run PGAP on a multi-sequence FASTA file using multiprocessing."
    )
    parser.add_argument("--fasta", required=True,
                        help="Input multi-sequence FASTA file.")
    parser.add_argument("--outdir", required=True,
                        help="Base output directory for PGAP results.")
    parser.add_argument("--submol", required=True,
                        help="Path to submol.yaml template.")
    parser.add_argument("--pgap-exe", required=True, dest="pgap_exe",
                        help="Path to pgap.py executable.")
    parser.add_argument("--container", required=True,
                        help="Path to PGAP Singularity .sif container.")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel worker processes. Default: 4")
    parser.add_argument("--cpus-per-worker", type=int, default=None,
                        dest="cpus_per_worker",
                        help="CPUs passed to PGAP (--cpu) per worker. "
                             "Default: total CPU count / --workers, rounded down (min 1).")
    parser.add_argument("--run-dir-base", dest="run_dir_base", default=None,
                        help="Base directory for temporary per-sample run folders. "
                             "Default: <outdir>/run_folders")
    parser.add_argument("--keep-run-dirs", dest="keep_run_dirs",
                        action="store_true", default=False,
                        help="Keep temporary run directories after completion.")
    parser.add_argument("--pgap-extra", dest="pgap_extra",
                        default="-n --no-internet --ignore-all-errors --no-self-update",
                        help="Extra flags passed verbatim to pgap.py.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def get_logger(name="pgap_mp"):
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                              datefmt="%Y-%m-%d %H:%M:%S")
        )
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def make_sample_name(seq_id):
    return re.sub(r"[^\w.\-]", "_", seq_id)


def is_already_done(outdir, sample_name):
    result_dir = os.path.join(outdir, f"{sample_name}_results")
    annot_faa = os.path.join(result_dir, "annot.faa")
    if os.path.isdir(result_dir):
        if os.path.isfile(annot_faa) and os.path.getsize(annot_faa) > 0:
            return True
        shutil.rmtree(result_dir)
    return False


def write_input_yaml(run_dir, fasta_filename):
    input_yaml_path = os.path.join(run_dir, "input.yaml")
    content = (
        "fasta:\n"
        "  class: File\n"
        f"  location: {fasta_filename}\n"
        "submol:\n"
        "  class: File\n"
        "  location: submol.yaml\n"
    )
    with open(input_yaml_path, "w") as fh:
        fh.write(content)
    return input_yaml_path


# ---------------------------------------------------------------------------
# Worker function (called in subprocess)
# ---------------------------------------------------------------------------

def process_sequence(task):
    """
    Process a single sequence: create run dir, run PGAP, return result dict.
    task is a dict produced by build_tasks().
    """
    seq_id = task["seq_id"]
    fasta_str = task["fasta_str"]
    sample_name = task["sample_name"]
    run_dir = task["run_dir"]
    output_dir = task["output_dir"]
    log_file = task["log_file"]
    submol_path = task["submol"]
    pgap_exe = task["pgap_exe"]
    container = task["container"]
    pgap_extra = task["pgap_extra"]
    keep_run_dirs = task["keep_run_dirs"]

    pid = os.getpid()
    logger = get_logger(f"pgap_mp.{pid}")

    logger.info("Processing %s (pid=%d)", seq_id, pid)

    # Create run directory
    os.makedirs(run_dir, exist_ok=True)
    fasta_filename = f"{sample_name}.fasta"
    with open(os.path.join(run_dir, fasta_filename), "w") as fh:
        fh.write(fasta_str)
    shutil.copy2(submol_path, os.path.join(run_dir, "submol.yaml"))
    input_yaml_path = write_input_yaml(run_dir, fasta_filename)

    # Run PGAP
    cmd = [pgap_exe] + pgap_extra.split() + [
        "--cpu", str(task["cpus_per_worker"]),
        "--container-path", container,
        "-o", output_dir,
        input_yaml_path,
    ]
    env = os.environ.copy()
    #env["PGAP_INPUT_DIR"] = run_dir

    with open(log_file, "w") as log_fh:
        result = subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT, env=env)

    returncode = result.returncode
    success = (returncode == 0)

    if success:
        logger.info("Finished %s (returncode=0).", seq_id)
        if not keep_run_dirs:
            shutil.rmtree(run_dir, ignore_errors=True)
    else:
        logger.warning("PGAP failed for %s (returncode=%d). Run dir: %s",
                       seq_id, returncode, run_dir)

    return {"seq_id": seq_id, "returncode": returncode, "success": success}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    logger = get_logger()

    # Validate inputs
    for label, path in [("--fasta", args.fasta), ("--submol", args.submol),
                        ("--pgap-exe", args.pgap_exe), ("--container", args.container)]:
        if not os.path.exists(path):
            logger.error("Path for %s does not exist: %s", label, path)
            sys.exit(1)

    run_dir_base = args.run_dir_base or os.path.join(args.outdir, "run_folders")
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "logs"), exist_ok=True)
    os.makedirs(run_dir_base, exist_ok=True)

    logger.info("Loading sequences from %s", args.fasta)
    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        logger.warning("Input FASTA is empty. Nothing to do.")
        return

    logger.info("Found %d sequences.", len(records))

    # Resolve parallelism settings before building task list
    workers = min(args.workers, len(records))
    total_cpus = os.cpu_count() or 1
    if args.cpus_per_worker is not None:
        cpus_per_worker = max(1, args.cpus_per_worker)
    else:
        cpus_per_worker = max(1, total_cpus // workers)
    logger.info("Workers: %d | CPUs per PGAP call: %d | Total CPUs on node: %d",
                workers, cpus_per_worker, total_cpus)

    # Build task list, skipping already-done sequences
    tasks = []
    skipped = 0
    for record in records:
        sample_name = make_sample_name(record.id)
        if is_already_done(args.outdir, sample_name):
            logger.info("Skipping %s (already complete).", sample_name)
            skipped += 1
            continue
        buf = io.StringIO()
        SeqIO.write(record, buf, "fasta")
        tasks.append({
            "seq_id": record.id,
            "fasta_str": buf.getvalue(),
            "sample_name": sample_name,
            "run_dir": os.path.join(run_dir_base, sample_name),
            "output_dir": os.path.join(args.outdir, f"{sample_name}_results"),
            "log_file": os.path.join(args.outdir, "logs", f"{sample_name}.log"),
            "submol": args.submol,
            "pgap_exe": args.pgap_exe,
            "container": args.container,
            "pgap_extra": args.pgap_extra,
            "keep_run_dirs": args.keep_run_dirs,
            "cpus_per_worker": cpus_per_worker,
        })

    total = len(tasks)
    logger.info("%d sequences to process, %d skipped.", total, skipped)

    if not tasks:
        return

    workers = min(workers, total)

    with multiprocessing.Pool(processes=workers) as pool:
        results = pool.map(process_sequence, tasks)

    completed = sum(1 for r in results if r["success"])
    failed = sum(1 for r in results if not r["success"])
    logger.info("Done. %d succeeded, %d failed, %d skipped.", completed, failed, skipped)

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()

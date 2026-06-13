#!/usr/bin/env python3
"""
MPI-based PGAP annotation pipeline.

Rank 0 acts as a controller: reads a multi-sequence FASTA file and dispatches
one sequence at a time to worker ranks. Workers create per-sample run directories,
run PGAP via subprocess, and report back when done. Work continues until all
sequences are processed.

Usage:
    mpirun -n <N> python run_pgap_mpi.py --fasta <fasta> --outdir <outdir> \
        --submol <submol.yaml> --pgap-exe <pgap.py> --container <pgap.sif>

Requirements:
    mpi4py, biopython
"""

import argparse
import io
import logging
import os
import re
import shutil
import subprocess
import sys

from mpi4py import MPI
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# MPI message tags
# ---------------------------------------------------------------------------
TAG_READY = 1  # Worker -> Controller: ready for work
TAG_WORK = 2   # Controller -> Worker: work item or STOP sentinel
TAG_DONE = 3   # Worker -> Controller: completion report

STOP_SENTINEL = None  # Sent by controller to signal worker shutdown


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run PGAP on a multi-sequence FASTA file using MPI."
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
    parser.add_argument("--run-dir-base", dest="run_dir_base", default=None,
                        help="Base directory for temporary per-sample run folders. "
                             "Default: <outdir>/run_folders")
    parser.add_argument("--keep-run-dirs", dest="keep_run_dirs",
                        action="store_true", default=False,
                        help="Keep temporary run directories after completion. "
                             "Default: delete on success.")
    parser.add_argument("--pgap-extra", dest="pgap_extra",
                        default="-n --no-internet --ignore-all-errors --no-self-update",
                        help="Extra flags passed verbatim to pgap.py.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def get_logger(rank):
    logger = logging.getLogger(f"pgap_mpi.rank{rank}")
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter(f"%(asctime)s [rank {rank}] %(levelname)s %(message)s",
                              datefmt="%Y-%m-%d %H:%M:%S")
        )
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def make_sample_name(seq_id):
    """Derive a filesystem-safe directory name from a sequence ID."""
    name = re.sub(r"[^\w.\-]", "_", seq_id)
    return name


def seq_record_to_fasta_str(record):
    """Serialize a SeqRecord to a FASTA-formatted string."""
    buf = io.StringIO()
    SeqIO.write(record, buf, "fasta")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Controller-side helpers (rank 0)
# ---------------------------------------------------------------------------

def load_sequences(fasta_path):
    """Parse multi-FASTA and return a list of SeqRecord objects."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    return records


def setup_output_dirs(outdir, run_dir_base):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    os.makedirs(run_dir_base, exist_ok=True)


def is_already_done(outdir, sample_name):
    """
    Mirror the completion check from create_pgap_batch.py.
    Returns True if annot.faa exists and is non-empty (skip).
    Removes the result directory if it exists but annot.faa is missing/empty.
    """
    result_dir = os.path.join(outdir, f"{sample_name}_results")
    annot_faa = os.path.join(result_dir, "annot.faa")
    if os.path.isdir(result_dir):
        if os.path.isfile(annot_faa) and os.path.getsize(annot_faa) > 0:
            return True
        # Incomplete — remove so PGAP can re-run cleanly
        shutil.rmtree(result_dir)
    return False


# ---------------------------------------------------------------------------
# Worker-side helpers (ranks 1+)
# ---------------------------------------------------------------------------

def create_run_directory(run_dir, fasta_str, submol_path):
    """
    Create a per-sample run directory containing:
      - single-sequence FASTA file
      - submol.yaml (copied from template)
      - input.yaml referencing those files
    Returns the path to input.yaml.
    """
    os.makedirs(run_dir, exist_ok=True)

    # Parse sequence ID from FASTA string to derive filename
    record = SeqIO.read(io.StringIO(fasta_str), "fasta")
    fasta_filename = f"{make_sample_name(record.id)}.fasta"
    fasta_dest = os.path.join(run_dir, fasta_filename)
    with open(fasta_dest, "w") as fh:
        fh.write(fasta_str)

    shutil.copy2(submol_path, os.path.join(run_dir, "submol.yaml"))

    input_yaml_path = write_input_yaml(run_dir, fasta_filename)
    return input_yaml_path


def write_input_yaml(run_dir, fasta_filename):
    """Write the two-field input.yaml expected by PGAP."""
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


def run_pgap(pgap_exe, container_path, input_yaml_path, output_dir, log_file, pgap_extra):
    """
    Run PGAP via subprocess, capturing all output to log_file.
    Returns the subprocess returncode.
    """
    cmd = [pgap_exe] + pgap_extra.split() + [
        "--container-path", container_path,
        "-o", output_dir,
        input_yaml_path,
    ]
    env = os.environ.copy()
    # PGAP may need PGAP_INPUT_DIR pointing to the run directory
    env["PGAP_INPUT_DIR"] = os.path.dirname(input_yaml_path)

    with open(log_file, "w") as log_fh:
        result = subprocess.run(
            cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            env=env,
        )
    return result.returncode


# ---------------------------------------------------------------------------
# Controller (rank 0)
# ---------------------------------------------------------------------------

def run_controller(comm, args):
    logger = get_logger(0)
    size = comm.Get_size()
    num_workers = size - 1

    if num_workers < 1:
        logger.error("Need at least 2 MPI ranks (1 controller + 1 worker). Got %d.", size)
        comm.Abort(1)

    run_dir_base = args.run_dir_base or os.path.join(args.outdir, "run_folders")

    # Validate inputs
    for label, path in [("--fasta", args.fasta), ("--submol", args.submol),
                        ("--pgap-exe", args.pgap_exe), ("--container", args.container)]:
        if not os.path.exists(path):
            logger.error("Path for %s does not exist: %s", label, path)
            comm.Abort(1)

    setup_output_dirs(args.outdir, run_dir_base)

    logger.info("Loading sequences from %s", args.fasta)
    records = load_sequences(args.fasta)
    if not records:
        logger.warning("Input FASTA is empty. Nothing to do.")
        # Send STOP to all workers immediately
        for _ in range(num_workers):
            comm.recv(source=MPI.ANY_SOURCE, tag=TAG_READY)
        for rank in range(1, size):
            comm.send(STOP_SENTINEL, dest=rank, tag=TAG_WORK)
        return

    logger.info("Found %d sequences.", len(records))

    # Build work queue, skipping already-completed sequences
    work_queue = []
    skipped = 0
    for idx, record in enumerate(records):
        sample_name = make_sample_name(record.id)
        if is_already_done(args.outdir, sample_name):
            logger.info("Skipping %s (already complete).", sample_name)
            skipped += 1
        else:
            fasta_str = seq_record_to_fasta_str(record)
            work_queue.append({"seq_id": record.id, "fasta_str": fasta_str,
                                "seq_index": idx})

    total = len(work_queue)
    logger.info("%d sequences to process, %d skipped.", total, skipped)

    completed = 0
    failed = 0
    queue_idx = 0
    workers_stopped = 0

    while workers_stopped < num_workers:
        # Wait for any worker to become ready
        status = MPI.Status()
        msg = comm.recv(source=MPI.ANY_SOURCE, tag=TAG_READY, status=status)
        worker_rank = status.Get_source()

        if queue_idx < len(work_queue):
            work_item = work_queue[queue_idx]
            queue_idx += 1
            comm.send(work_item, dest=worker_rank, tag=TAG_WORK)
            logger.info("[%d/%d] Dispatched %s to rank %d.",
                        queue_idx, total, work_item["seq_id"], worker_rank)
        else:
            # No more work — send STOP
            comm.send(STOP_SENTINEL, dest=worker_rank, tag=TAG_WORK)
            workers_stopped += 1

    # Collect DONE messages for any in-flight work items before workers exit
    # (Workers send DONE before sending the final READY that triggers STOP)
    # Summary already updated inside the DONE receive below.

    logger.info("All workers stopped. Summary: %d completed, %d failed, %d skipped.",
                completed, failed, skipped)


# ---------------------------------------------------------------------------
# Worker (ranks 1+)
# ---------------------------------------------------------------------------

def run_worker(comm, args):
    rank = comm.Get_rank()
    logger = get_logger(rank)
    run_dir_base = args.run_dir_base or os.path.join(args.outdir, "run_folders")

    while True:
        try:
            # Signal readiness
            comm.send({"rank": rank}, dest=0, tag=TAG_READY)
            work = comm.recv(source=0, tag=TAG_WORK)

            if work is STOP_SENTINEL:
                logger.info("Received STOP. Exiting.")
                break

            seq_id = work["seq_id"]
            fasta_str = work["fasta_str"]
            sample_name = make_sample_name(seq_id)

            run_dir = os.path.join(run_dir_base, sample_name)
            output_dir = os.path.join(args.outdir, f"{sample_name}_results")
            log_file = os.path.join(args.outdir, "logs", f"{sample_name}.log")

            logger.info("Processing %s", seq_id)

            input_yaml_path = create_run_directory(run_dir, fasta_str, args.submol)

            returncode = run_pgap(
                pgap_exe=args.pgap_exe,
                container_path=args.container,
                input_yaml_path=input_yaml_path,
                output_dir=output_dir,
                log_file=log_file,
                pgap_extra=args.pgap_extra,
            )

            success = (returncode == 0)
            if success:
                logger.info("Finished %s (returncode=0).", seq_id)
                if not args.keep_run_dirs:
                    shutil.rmtree(run_dir, ignore_errors=True)
            else:
                logger.warning("PGAP failed for %s (returncode=%d). "
                               "Run dir preserved: %s", seq_id, returncode, run_dir)

            comm.send(
                {"seq_id": seq_id, "returncode": returncode,
                 "success": success, "skipped": False},
                dest=0, tag=TAG_DONE,
            )

        except Exception as exc:
            logger.error("Unhandled exception for rank %d: %s", rank, exc, exc_info=True)
            # Best-effort DONE so controller is not left waiting
            try:
                comm.send(
                    {"seq_id": "unknown", "returncode": -1,
                     "success": False, "skipped": False},
                    dest=0, tag=TAG_DONE,
                )
            except Exception:
                pass
            # Re-signal ready so controller can send STOP cleanly
            try:
                comm.send({"rank": rank}, dest=0, tag=TAG_READY)
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    args = parse_args()

    if rank == 0:
        run_controller(comm, args)
    else:
        run_worker(comm, args)


if __name__ == "__main__":
    main()

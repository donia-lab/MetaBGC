#!/usr/bin/env python3
"""
For each <id>_seg.fasta in blastx_result/, write a subset fasta containing
only records whose IDs are NOT present in the corresponding <id>_seg.txt
(first column, tab-delimited).  Output files go to blastx_result/unmatched/.
"""

import os
import glob

INPUT_DIR = "blastx_result"
OUTPUT_DIR = os.path.join(INPUT_DIR, "unmatched")
COMMANDS_FILE = os.path.join(OUTPUT_DIR, "blastx_commands.sh")
os.makedirs(OUTPUT_DIR, exist_ok=True)

BLAST_CMD = (
    'blastx -num_threads 4 -query {query} -db nr -max_target_seqs 1 '
    '-outfmt "6 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle" '
    '-out {out}'
)

fasta_files = glob.glob(os.path.join(INPUT_DIR, "*_seg.fasta"))

for fasta_path in sorted(fasta_files):
    base = os.path.basename(fasta_path).replace(".fasta", "")  # e.g. 0_seg
    txt_path = os.path.join(INPUT_DIR, base + ".txt")

    if not os.path.exists(txt_path):
        print(f"WARNING: no txt file for {fasta_path}, skipping")
        continue

    # collect IDs present in the txt file (first tab-delimited column)
    matched_ids = set()
    with open(txt_path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                matched_ids.add(line.split("\t")[0])

    # stream through fasta, writing records whose ID is not in matched_ids
    out_path = os.path.join(OUTPUT_DIR, base + "_unmatched.fasta")
    written = 0
    current_header = None
    current_seq_lines = []
    in_keep = False

    def flush(out_fh, header, seq_lines):
        out_fh.write(header + "\n")
        out_fh.writelines(seq_lines)

    with open(fasta_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            if line.startswith(">"):
                # flush previous record if kept
                if in_keep and current_header is not None:
                    flush(out_fh, current_header, current_seq_lines)
                    written += 1
                record_id = line[1:].strip().split()[0]
                in_keep = record_id not in matched_ids
                current_header = line.rstrip()
                current_seq_lines = []
            else:
                if in_keep:
                    current_seq_lines.append(line)
        # flush last record
        if in_keep and current_header is not None:
            flush(out_fh, current_header, current_seq_lines)
            written += 1

    print(f"{base}: {written} unmatched records -> {out_path}")

commands = []
for fasta_path in sorted(glob.glob(os.path.join(OUTPUT_DIR, "*_seg_unmatched.fasta"))):
    fname = os.path.basename(fasta_path)
    out_name = fname.replace(".fasta", ".txt")
    commands.append(BLAST_CMD.format(
        query=os.path.join(OUTPUT_DIR, fname),
        out=os.path.join(OUTPUT_DIR, out_name),
    ))

with open(COMMANDS_FILE, "w") as fh:
    fh.write("\n".join(commands) + "\n")

print(f"\n{len(commands)} blastx commands written to {COMMANDS_FILE}")

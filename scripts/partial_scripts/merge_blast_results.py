#!/usr/bin/env python3
"""
For each <id>_seg.txt in blastx_result/, concatenate it with the
corresponding <id>_seg_unmatched.txt in blastx_result/unmatched/
(the tblastn results for reads that had no original blastx hit).
Merged output goes to blastx_result/merged/<id>_seg.txt.
"""

import os
import glob

INPUT_DIR = "blastx_result"
UNMATCHED_DIR = os.path.join(INPUT_DIR, "unmatched")
OUTPUT_DIR = os.path.join(INPUT_DIR, "merged")
os.makedirs(OUTPUT_DIR, exist_ok=True)

txt_files = glob.glob(os.path.join(INPUT_DIR, "*_seg.txt"))

for txt_path in sorted(txt_files):
    base = os.path.basename(txt_path).replace(".txt", "")  # e.g. 0_seg
    unmatched_path = os.path.join(UNMATCHED_DIR, base + "_unmatched.txt")
    out_path = os.path.join(OUTPUT_DIR, base + ".txt")

    n_orig = 0
    n_unmatched = 0

    with open(out_path, "w") as out_fh:
        with open(txt_path) as fh:
            for line in fh:
                out_fh.write(line)
                n_orig += 1

        if os.path.exists(unmatched_path):
            with open(unmatched_path) as fh:
                for line in fh:
                    out_fh.write(line)
                    n_unmatched += 1
        else:
            print(f"WARNING: no unmatched txt for {base}")

    print(f"{base}: {n_orig} original + {n_unmatched} unmatched -> {out_path}")

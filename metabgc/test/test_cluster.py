"""End-to-end unit test for the clustering step.

Builds a tiny abundance matrix in which two reads have perfectly
correlated abundance profiles (and a third is uncorrelated), then checks
that DBSCAN groups the correlated pair together and that the expected
output files are produced.
"""
import json

import pandas as pd

from metabgc.src.metabgccluster import mbgccluster


def _write_inputs(tmp_path):
    # Read-level abundance matrix: first column is the read id, the rest are
    # per-sample abundances. r2 == 2 * r1 (correlation 1); r3 is uncorrelated.
    matrix = tmp_path / "abundance.txt"
    matrix.write_text(
        "qseqid\tS1\tS2\tS3\tS4\n"
        "r1\t10\t0\t5\t1\n"
        "r2\t20\t0\t10\t2\n"
        "r3\t0\t10\t0\t5\n"
    )

    # Long/wide abundance table: one row per (read, sample) with a count.
    wide = tmp_path / "abundance_wide.txt"
    wide.write_text(
        "qseqid\tSample\tcount\n"
        "r1\tS1\t10\n"
        "r1\tS3\t5\n"
        "r1\tS4\t1\n"
        "r2\tS1\t20\n"
        "r2\tS3\t10\n"
        "r2\tS4\t2\n"
        "r3\tS2\t10\n"
        "r3\tS4\t5\n"
    )

    fasta = tmp_path / "identified-biosynthetic-reads.fasta"
    fasta.write_text(">r1\nACGTACGT\n>r2\nTTTTGGGG\n>r3\nCCCCAAAA\n")
    return matrix, wide, fasta


def test_mbgccluster_groups_correlated_reads(tmp_path):
    matrix, wide, fasta = _write_inputs(tmp_path)

    summary_file, abund_file = mbgccluster(
        str(matrix),
        str(wide),
        str(fasta),
        max_dist=0.1,
        min_samples=1,
        readThresh=2,
        abundThresh=1,
        cpu=1,
    )

    labels = json.loads((tmp_path / "abundance_DBSCAN.json").read_text())
    # r1 and r2 are perfectly correlated -> same bin; r3 is on its own.
    assert labels["r1"] == labels["r2"]
    assert labels["r3"] != labels["r1"]

    summary_text = (tmp_path / "BinSummary.txt").read_text()
    assert "Number of Bins: 2" in summary_text

    # The correlated pair (>= readThresh reads) is written to gt_thresh.
    gt_thresh = tmp_path / "bin_fasta" / "gt_thresh"
    assert gt_thresh.is_dir()
    assert any(gt_thresh.iterdir())

    assert summary_file and abund_file

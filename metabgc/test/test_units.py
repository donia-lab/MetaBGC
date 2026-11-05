"""Self-contained unit tests for MetaBGC.

These tests build their own tiny fixtures and exercise the pure-Python
data-transformation helpers only. They do NOT require the external
bioinformatics binaries (hmmsearch, blastn, muscle, transeq, cd-hit, art)
or the large integration datasets, so they run anywhere the Python
dependencies are installed.
"""
import os

import pandas as pd
import pytest

from metabgc.src import metabgcbuild as build
from metabgc.src.metabgcidentify import filter_spHMM_data, runidentify
from metabgc.src.metabgcquantify import combine_blast_results, create_clustering_file
from metabgc.src.evaluate_sphmms import getOverlap
from metabgc.src.utils import interleave, createPandaDF
from metabgc.src.hmmrecord import HMMRecord


# ---------------------------------------------------------------------------
# metabgcbuild.ungappedseqsearch
# ---------------------------------------------------------------------------
def test_ungappedseqsearch_finds_match():
    assert build.ungappedseqsearch("--MSE-HDTDV---LVGGSM", "TDV-LV") == [8, 16]
    assert build.ungappedseqsearch(
        "HLAGEPLETCLRYGAIAGAYACTIPATRAGAIDRAALLRPAA--------",
        "ACTIPATRAGAIDRAALLRPAA--------",
    ) == [20, 42]


def test_ungappedseqsearch_no_match():
    assert build.ungappedseqsearch("ACDE", "WWWW") == [-1, -1]


# ---------------------------------------------------------------------------
# utils.interleave
# ---------------------------------------------------------------------------
def test_interleave_alternates_records():
    assert list(interleave([1, 3, 5], [2, 4, 6])) == [1, 2, 3, 4, 5, 6]


def test_interleave_stops_at_shortest():
    assert list(interleave(["f1", "f2", "f3"], ["r1"])) == ["f1", "r1"]


# ---------------------------------------------------------------------------
# utils.createPandaDF
# ---------------------------------------------------------------------------
def test_createPandaDF_sorts_by_score_descending(tmp_path):
    out = tmp_path / "hmm.txt"
    hmm_dict = {
        "readA": HMMRecord("readA", "T", "S1", "P", 50.0, "30", "0_30"),
        "readB": HMMRecord("readB", "T", "S1", "P", 80.0, "30", "0_30"),
    }
    createPandaDF(hmm_dict, str(out))

    df = pd.read_csv(out, sep="\t", header=None)
    assert len(df) == 2
    # highest HMMScore (column index 4) must come first
    assert df.iloc[0, 0] == "readB"
    assert df.iloc[0, 4] == 80.0
    assert df.iloc[1, 0] == "readA"


def test_createPandaDF_empty_dict_writes_empty_file(tmp_path):
    out = tmp_path / "empty.txt"
    createPandaDF({}, str(out))
    assert out.exists()
    assert out.stat().st_size == 0


# ---------------------------------------------------------------------------
# metabgcidentify.filter_spHMM_data
# ---------------------------------------------------------------------------
def _make_hmm_df():
    return pd.DataFrame(
        {
            "readID": ["r1", "r2", "r3"],
            "HMMScore": [100.0, 40.0, 90.0],
            "interval": ["0_30", "0_30", "30_60"],
        }
    )


def test_filter_spHMM_data_applies_per_interval_cutoffs():
    spHMM_df = _make_hmm_df()
    cutoff_df = pd.DataFrame({"interval": ["0_30", "30_60"], "cutoff": [50.0, 50.0]})
    out = filter_spHMM_data(spHMM_df, cutoff_df)
    # r1 (100 >= 50) kept, r2 (40 < 50) dropped, r3 (90 >= 50) kept
    assert set(out["readID"]) == {"r1", "r3"}


def test_filter_spHMM_data_empty_cutoffs_returns_empty_frame():
    spHMM_df = _make_hmm_df()
    out = filter_spHMM_data(spHMM_df, spHMM_df.iloc[0:0][["interval"]].assign(cutoff=[]))
    assert len(out) == 0
    assert list(out.columns) == list(spHMM_df.columns)


# ---------------------------------------------------------------------------
# metabgcidentify.runidentify (end-to-end file IO)
# ---------------------------------------------------------------------------
def test_runidentify_dereplicates_frames_and_applies_cutoff(tmp_path):
    hmm_file = tmp_path / "CombinedHmmSearch.txt"
    # columns: readID sampleType Sample protType HMMScore window interval
    hmm_file.write_text(
        "read1_1\tALL\tS1\tP\t100\t30\t0_30\n"
        "read1_2\tALL\tS1\tP\t120\t30\t0_30\n"
        "read2_1\tALL\tS1\tP\t20\t30\t0_30\n"
    )
    cutoff_file = tmp_path / "cutoff.tsv"
    cutoff_file.write_text("interval\tcutoff\n0_30\t50\n")

    filtered = tmp_path / "spHMM-filtered-results.txt"
    read_ids = tmp_path / "CombinedReadIds.txt"
    runidentify(str(hmm_file), str(cutoff_file), str(filtered), str(read_ids))

    # The two frames of read1 collapse to a single read id; read2 is below cutoff.
    ids = [line.strip() for line in read_ids.read_text().splitlines() if line.strip()]
    assert ids == ["read1"]

    df = pd.read_csv(filtered, sep="\t")
    assert set(df["readID"]) == {"read1"}
    # highest HMM score across frames is retained
    assert df["HMMScore"].max() == 120


# ---------------------------------------------------------------------------
# evaluate_sphmms.getOverlap
# ---------------------------------------------------------------------------
def test_getOverlap_overlapping_interval():
    # overlap of [10,20] with [5,15] is positions 10..15 -> 6 (inclusive count)
    assert getOverlap({"qstart": 10, "qend": 20}, [5, 15]) == 6


# ---------------------------------------------------------------------------
# metabgcquantify.combine_blast_results + create_clustering_file
# ---------------------------------------------------------------------------
def _write_blast_hit(qseqid, sample_line_count):
    # 11 whitespace-separated columns:
    # sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue
    line = "subj 100 1 100 {q} 150 1 100 90 99.0 1e-50\n".format(q=qseqid)
    return line * sample_line_count


def test_combine_blast_results_concatenates_and_tags_samples(tmp_path):
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    (blast_dir / "sampleA.txt").write_text(_write_blast_hit("read1", 2))
    (blast_dir / "sampleB.txt").write_text(_write_blast_hit("read1", 1))

    combined = tmp_path / "combined.txt"
    n_rows = combine_blast_results(str(blast_dir), str(combined), "cohortX")
    assert n_rows == 3

    df = pd.read_csv(combined, sep="\t", header=None)
    assert df.shape[0] == 3
    # last two appended columns are Sample and cohort
    assert set(df.iloc[:, -2]) == {"sampleA", "sampleB"}
    assert set(df.iloc[:, -1]) == {"cohortX"}


def test_create_clustering_file_builds_abundance_matrix(tmp_path):
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    (blast_dir / "sampleA.txt").write_text(_write_blast_hit("read1", 2))
    (blast_dir / "sampleB.txt").write_text(_write_blast_hit("read1", 1))
    combined = tmp_path / "combined.txt"
    combine_blast_results(str(blast_dir), str(combined), "cohortX")

    abund = tmp_path / "abundance.txt"
    abund_wide = tmp_path / "abundance_wide.txt"
    create_clustering_file(str(combined), str(abund), str(abund_wide))

    matrix = pd.read_csv(abund, sep="\t", index_col=0)
    # read1 appears twice in sampleA and once in sampleB
    assert matrix.loc["read1", "sampleA"] == 2
    assert matrix.loc["read1", "sampleB"] == 1

    wide = pd.read_csv(abund_wide, sep="\t")
    assert set(wide.columns) == {"Sample", "qseqid", "count"}
    assert wide["count"].sum() == 3

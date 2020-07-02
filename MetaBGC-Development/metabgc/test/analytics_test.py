import pytest
from metabgc.src.metabgcanalytics import mbgcanalytics


def test1():
    data_dir = "AbcK/output"
    cohort_metadata_file = "../metabgc/metadata/combined_cohort.csv"
    BLASTDB = "nr"
    mbgcanalytics(BLASTDB, data_dir, cohort_metadata_file, data_dir, 1)
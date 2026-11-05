import os as _os
import pytest as _pytest

# These are integration tests that require the large example datasets
# (AbcK / RdhA / OxyN, distributed as zip archives in this folder) to be
# extracted next to this file, plus the external bioinformatics binaries.
# They are skipped automatically when that data is not present.
if not _os.path.isdir(_os.path.join(_os.path.dirname(__file__), "AbcK")):
    _pytest.skip("integration datasets not extracted", allow_module_level=True)

from metabgc.src.metabgccluster import *

def test_cluster():
    cohortStr = "MetaHit"
    table = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table.txt"
    tableAbundance = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table-wide.txt"
    identifiedReadFile = "AbcK/output/"+cohortStr+"/identify/identified-biosynthetic-reads.fasta"
    mbgccluster(table, tableAbundance, identifiedReadFile, 0.1, 1, 10, 10, 4)

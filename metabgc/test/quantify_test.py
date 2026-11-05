import os as _os
import pytest as _pytest

# These are integration tests that require the large example datasets
# (AbcK / RdhA / OxyN, distributed as zip archives in this folder) to be
# extracted next to this file, plus the external bioinformatics binaries.
# They are skipped automatically when that data is not present.
if not _os.path.isdir(_os.path.join(_os.path.dirname(__file__), "AbcK")):
    _pytest.skip("integration datasets not extracted", allow_module_level=True)

from metabgc.src.metabgcquantify import *

def test_combine():
    cohortStr = "MetaHit"
    blastn_search_directory = "AbcK/data/"+cohortStr+"/quantify_blastn_result"
    abundFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table.txt"
    abundWideFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table-wide.txt"
    combinedBLASTPath = "AbcK/output/"+cohortStr+"/CombinedQuantifyBLAST.txt"
    combine_blast_results(blastn_search_directory, combinedBLASTPath, cohortStr)

def test_clustering():
    cohortStr = "MetaHit"
    abundFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table.txt"
    abundWideFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table-wide.txt"
    combinedBLASTPath = "AbcK/output/"+cohortStr+"/CombinedQuantifyBLAST.txt"
    create_clustering_file(combinedBLASTPath,abundFile,abundWideFile)
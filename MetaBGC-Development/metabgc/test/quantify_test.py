import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcquantify import *

def test_combine():
    cohortStr = "MetaHit"
    blastn_search_directory = "AbcK/data/"+cohortStr+"/quantify_blastn_result"
    abundFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table.txt"
    abundWideFile = "AbcK/output/"+cohortStr+"/unique-biosynthetic-reads-abundance-table-wide.txt"
    combinedBLASTPath = "AbcK/output/"+cohortStr+"/CobminedQuantifyBLAST.txt"
    combine_blast_results(blastn_search_directory, combinedBLASTPath, cohortStr)
    create_clustering_file(combinedBLASTPath,abundFile,abundWideFile)
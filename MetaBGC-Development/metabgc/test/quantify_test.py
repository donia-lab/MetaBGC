import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcquantify import *

def test_combine_2():
    sampleStr = "combined_2"
    blastn_search_directory="C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr+"/results"
    quant_op_dir = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr
    combinedBLASTPath = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr + "/CobminedQuantifyBLAST.txt"
    combinedBLASTPath = combine_blast_results(blastn_search_directory, quant_op_dir, sampleStr)
    abund_file, abund_wide_file = create_clustering_file(quant_op_dir, combinedBLASTPath)

def test_combine():
    cohortStr = "MetaHit"
    blastn_search_directory = "AbcK/data/"+cohortStr+"/quantify_blastn_result"
    quant_op_dir = "AbcK/output/"+cohortStr
    combinedBLASTPath = combine_blast_results(blastn_search_directory, quant_op_dir, cohortStr)
    abund_file, abund_wide_file = create_clustering_file(quant_op_dir, combinedBLASTPath)
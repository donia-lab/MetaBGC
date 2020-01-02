import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcquantify import *

def test_combine():
    sampleStr = "MetaHit"
    blastn_search_directory="C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr
    quant_op_dir = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr
    combinedBLASTPath = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr + "/CobminedQuantifyBLAST.txt"
    combinedBLASTPath = combine_blast_results(blastn_search_directory, quant_op_dir, sampleStr)
    abund_file = create_clustering_file(quant_op_dir, combinedBLASTPath)
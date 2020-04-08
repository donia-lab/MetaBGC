import pytest
import os
from metabgc.src.metabgccluster import *
import pandas as pd
import numpy as np
import statistics

#def test_cluster():
#    cohortStr = "MetaHit"
#    table = "AbcK/output/" + cohortStr +"/unique-biosynthetic-reads-abundance-table.txt"
#    tableAbundance = "AbcK/output/" + cohortStr + "/unique-biosynthetic-reads-abundance-table-wide.txt"
#    identifiedReadFile = "AbcK/output/" + cohortStr +"/identify/identified-biosynthetic-reads.fasta"
#    mbgccluster(table, tableAbundance, identifiedReadFile, 0.1, 1, 10, 10, 4)

def test_cluster():
    table = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/combined_2/unique-biosynthetic-reads-abundance-table.txt"
    tableAbundance = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/combined_2/unique-biosynthetic-reads-abundance-table-wide.txt"
    identifiedReadFile = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/combined_2/identified-biosynthetic-reads.fasta"
    mbgccluster(table, tableAbundance, identifiedReadFile, 0.1, 1, 10, 10, 4)

import pytest
import os
from metabgc.src.metabgccluster import *


def test_cluster():
    sampleStr = "HMP-1_2"
    table="C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/"+sampleStr+"/unique-biosynthetic-reads-abundance-table.txt"
    mbgccluster(table, 0.1, 1, 4)
import pytest
import os
from metabgc.src.metabgccluster import *


def test_cluster():
    cohortStr = "MetaHit"
    table = "AbcK/output/" + cohortStr +"/unique-biosynthetic-reads-abundance-table.txt"
    tableAbundance = "AbcK/output/" + cohortStr + "/unique-biosynthetic-reads-abundance-table-wide.txt"
    identifiedReadFile = "AbcK/output/" + cohortStr +"/identify/identified-biosynthetic-reads.fasta"
    mbgccluster(table, tableAbundance, identifiedReadFile, 0.1, 1, 10, 10, 4)
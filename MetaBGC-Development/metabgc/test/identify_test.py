import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcidentify import *
from metabgc.src.extractfastaseq import *

def test_parseHMM():
    cohortStr = "MetaHit"
    hmmPathFile = "AbcK/data/"+cohortStr+"/hmmsearch_result/G30365-merged__160_190.tbl"
    ouputDir = "AbcK/output/"+cohortStr
    parsedOuputDir = "AbcK/output/"+cohortStr+"/parsed_tables"
    os.makedirs(parsedOuputDir, 0o777, True)
    protType="AbcK"
    hmm_rec_dict = parseHMM(hmmPathFile, "hmmer3-tab", "ALL", "G30365-merged", protType, "30", "10")
    assert len(hmm_rec_dict) == 7

def test_identify():
    cohortStr = "MetaHit"
    allHMMResult = "AbcK/data/CombinedHmmSearch.txt"
    cutoff_file = "AbcK/output/build/HiPer_spHMMs/AbcK_F1_Cutoff.tsv"
    filteredTableFile="AbcK/output/"+cohortStr+"/identify/spHMM-filtered-results.txt"
    identifyReadIdFile="AbcK/output/"+cohortStr+"/identify/CombinedReadIds.txt"
    runidentify(allHMMResult, cutoff_file, filteredTableFile,identifyReadIdFile)

def test_extract():
    cohortStr = "MetaHit"
    filteredTableFile = "AbcK/data/"+ cohortStr +"/identify/spHMM-filtered-results.txt"
    readsDir = "AbcK/data/"+cohortStr+"/nucl"
    outputDir = "AbcK/output/"+cohortStr+"/identify/reads"
    identifyOutFile = "AbcK/output/"+cohortStr+"/identify/identified-biosynthetic-reads.fasta"
    os.makedirs(outputDir, 0o777, True)
    RunExtractDirectoryPar(readsDir, filteredTableFile, outputDir, identifyOutFile, "fasta", ncpus=4)


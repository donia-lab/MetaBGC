import os as _os
import pytest as _pytest

# These are integration tests that require the large example datasets
# (AbcK / RdhA / OxyN, distributed as zip archives in this folder) to be
# extracted next to this file, plus the external bioinformatics binaries.
# They are skipped automatically when that data is not present.
if not _os.path.isdir(_os.path.join(_os.path.dirname(__file__), "AbcK")):
    _pytest.skip("integration datasets not extracted", allow_module_level=True)

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
    cutoff_file = "AbcK/output/build/HiPerf/AbcK_F1_Cutoff.tsv"
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
    RunExtractDirectoryPar(readsDir, filteredTableFile, outputDir, identifyOutFile, "fasta", True, True, ncpus=4)


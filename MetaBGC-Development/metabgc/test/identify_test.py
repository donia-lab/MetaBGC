import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcidentify import *
from metabgc.src.extractfastaseq import *


def test_parseHMM():
    cohortStr = "MetaHit"
    input_dir = "AbcK/data/"+cohortStr+"/hmmsearch_result"
    ouputDir = "AbcK/output/"+cohortStr
    parsedOuputDir = "AbcK/output/"+cohortStr+"/parsed_tables"
    os.makedirs(parsedOuputDir, 0o777, True)
    window = "30_10"
    protType="AbcK"
    for subdir, dirs, files in os.walk(input_dir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.tbl$", file) and os.path.getsize(filePath) > 0 and "parsed_tables" not in filePath:
                hmm_tok = os.path.splitext(file)[0].split("_")
                interval = hmm_tok[1] + "_" + hmm_tok[2]
                read_sample = hmm_tok[0]
                result_dict = parseHMM(filePath, 'hmmer3-tab',cohortStr, read_sample, protType, window, interval)
                hmmSearchFileName = read_sample + "__" + interval + ".txt"
                hmmSearchFilePath = os.path.join(parsedOuputDir, hmmSearchFileName)
                createPandaDF(result_dict, hmmSearchFilePath)
    allHMMResult = os.path.join(ouputDir, "CombinedHmmSearch.txt")
    with open(allHMMResult, 'w') as outfile:
        for subdir, dirs, files in os.walk(parsedOuputDir):
            for file in files:
                filePath = os.path.join(subdir, file)
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)

def test_identify():
    cohortStr = "MetaHit"
    allHMMResult = "AbcK/output/"+cohortStr+"/CombinedHmmSearch.txt"
    cutoff_file = "AbcK/output/build/HiPer_spHMMs/AbcK_F1_Cutoff.tsv"
    filteredTableFile="AbcK/output/"+cohortStr+"/identify/spHMM-filtered-results.txt"
    identifyReadIdFile="AbcK/output/"+cohortStr+"/identify/CombinedReadIds.txt"
    runidentify(allHMMResult, cutoff_file, filteredTableFile,identifyReadIdFile)

def test_extract():
    cohortStr = "MetaHit"
    filteredTableFile = "AbcK/output/"+cohortStr+"/identify/spHMM-filtered-results.txt"
    readsDir = "AbcK/data/"+cohortStr+"/nucl"
    outputDir = "AbcK/output/"+cohortStr+"/identify/reads"
    identifyOutFile = "AbcK/output/"+cohortStr+"/identify/identified-biosynthetic-reads.fasta"
    os.makedirs(outputDir, 0o777, True)
    RunExtractDirectoryPar(readsDir, filteredTableFile, outputDir, identifyOutFile, ncpus=4)

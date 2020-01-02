import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcidentify import *
from metabgc.src.extractfastaseq import *


def test_parseHMM():
    sampleStr = "MetaHit"
    input_dir = "AbcK/data/identify/hmmsearch_result"
    ouputDir = "AbcK/output/MetaHit/identify/parsed_tables"
    os.makedirs(ouputDir, 0o777, True)
    window = "30_10"
    protType="AbcK"
    for subdir, dirs, files in os.walk(input_dir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.txt$", file) and os.path.getsize(filePath) > 0 and "parsed_tables" not in filePath:
                hmm_dir = os.path.dirname(filePath)
                interval = hmm_dir.split("__")[2].split('-')[0]
                read_sample = file.split("-against-")[0]
                result_dict = parseHMM(filePath, 'hmmer3-text',sampleStr, read_sample, protType, window, interval)
                hmmSearchFileName = read_sample + "_" + interval + ".txt"
                hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
                createPandaDF(result_dict, hmmSearchFilePath)
    allHMMResult = "AbcK/output/MetaHit/identify/CombinedHmmSearch.txt"
    with open(allHMMResult, 'w') as outfile:
        for subdir, dirs, files in os.walk(ouputDir):
            for file in files:
                filePath = os.path.join(subdir, file)
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)

def test_identify():
    sampleStr = "MetaHit"
    allHMMResult = "AbcK/output/MetaHit/identify/CombinedHmmSearch.txt"
    cutoff_file = "AbcK/output/build/HiPer_spHMMs/AbcK_F1_Cutoff.tsv"
    filteredTableFile="AbcK/output/MetaHit/identify/spHMM-filtered-results.txt"
    identifyReadIdFile="AbcK/output/MetaHit/identify/CombinedReadIds.txt"
    runidentify(allHMMResult, cutoff_file, filteredTableFile,identifyReadIdFile)


def test_extract():
    filteredTableFile = "AbcK/output/MetaHit/identify/spHMM-filtered-results.txt"
    readsDir = "AbcK/data/reads"
    outputDir = "AbcK/output/MetaHit/identify/reads"
    os.makedirs(outputDir, 0o777, True)
    RunExtractDirectoryPar(readsDir, filteredTableFile, outputDir, "identified-biosynthetic-reads.fasta", ncpus=4)

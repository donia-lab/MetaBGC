import pytest
import os
import re
from metabgc.src.utils import *
from metabgc.src.metabgcidentify import *
from metabgc.src.extractfastaseq import *


def test_parseHMM():
    sampleStr = "MetaHit"
    # input_dir = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/identify/"+sampleStr+"/hmmsearch_result"
    # ouputDir = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/identify/"+sampleStr+"/parsed_tables"
    # os.makedirs(ouputDir, 0o777, True)
    # window = "30_10"
    # protType="AbcK"
    # for subdir, dirs, files in os.walk(input_dir):
    #     for file in files:
    #         filePath = os.path.join(subdir, file)
    #         if re.match(r".*.txt$", file) and os.path.getsize(filePath) > 0 and "parsed_tables" not in filePath:
    #             hmm_dir = os.path.dirname(filePath)
    #             interval = hmm_dir.split("__")[2].split('-')[0]
    #             read_sample = file.split("-against-")[0]
    #             result_dict = parseHMM(filePath, 'hmmer3-text',sampleStr, read_sample, protType, window, interval)
    #             hmmSearchFileName = read_sample + "_" + interval + ".txt"
    #             hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
    #             createPandaDF(result_dict, hmmSearchFilePath)
    allHMMResult = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/identify/"+sampleStr+"/CombinedHmmSearch.txt"
    # with open(allHMMResult, 'w') as outfile:
    #     for subdir, dirs, files in os.walk(ouputDir):
    #         for file in files:
    #             filePath = os.path.join(subdir, file)
    #             if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
    #                 with open(filePath) as infile:
    #                     for line in infile:
    #                         outfile.write(line)

    cutoff_file = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/build/HiPer_spHMMs/AbcK_F1_Cutoff.tsv"
    filteredTableFile="C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/identify/"+sampleStr+"/spHMM-filtered-results.txt"
    identifyReadIdFile="C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/identify/"+sampleStr+"/CombinedReadIds.txt"
    runidentify(allHMMResult, cutoff_file, filteredTableFile,identifyReadIdFile)
    RunExtractDirectoryPar("asdasda", filteredTableFile, "xxx", "xxx", ncpus=4)
#    input_dir = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/MetaHIT/"
#    ouputDir = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/MetaHIT/parsed_tables"
#    os.makedirs(ouputDir, 0o777, True)
#    window = "30_10"
#    sampleType="HMP"
#    protType="AbcK"
#    for subdir, dirs, files in os.walk(input_dir):
#        for file in files:
#            filePath = os.path.join(subdir, file)
#            if re.match(r".*.txt$", file) and os.path.getsize(filePath) > 0 and "parsed_tables" not in filePath:
#                hmm_dir = os.path.dirname(filePath)
#                interval = hmm_dir.split("__")[2].split('-')[0]
#                sampleStr = file.split("-against-")[0]
#                result_dict = parseHMM(filePath, 'hmmer3-text',sampleType, sampleStr, protType, window, interval)
#                hmmSearchFileName = sampleStr + "_" + interval + ".txt"
#                hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
#                createPandaDF(result_dict, hmmSearchFilePath)
#
#    allHMMResult = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/MetaHIT/CombinedHmmSearch.txt"
#    with open(allHMMResult, 'w') as outfile:
#        for subdir, dirs, files in os.walk(ouputDir):
#            for file in files:
#                filePath = os.path.join(subdir, file)
#                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
#                    with open(filePath) as infile:
#                        for line in infile:
#                            outfile.write(line)
#
#    sampleStr = "MetaHIT"
#    identify_directory = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/identify"
#    fasta_dir = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/read_files"
#    allHMMResult = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/CombinedHmmSearch.txt"
#    os.makedirs(identify_directory, 0o777, True)
#    os.makedirs(fasta_dir, 0o777, True)
#    cutoff_file = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/build/HiPer_spHMMs/AbcK_F1_Cutoff.tsv"
#    runidentify(allHMMResult, identify_directory, cutoff_file, fasta_dir)

#    identifyReadIds = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/identify/CombinedReadIds.txt"
#    with open(identifyReadIds, 'w') as outfile:
#        for filename in os.listdir(fasta_dir):
#            if filename.endswith(".txt"):
#                filePath = os.path.join(fasta_dir, filename)
#                with open(filePath) as infile:
#                   for line in infile:
#                        outfile.write(line)
#    multiFastaFile = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/identify/identified-biosynthetic-reads.fasta"
#    nucl_seq_directory = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/nucl"
#    extract_read_directory = "/tigress/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/"+sampleStr+"/identify/reads"
#    RunExtractDirectoryPar(nucl_seq_directory, identifyReadIds, extract_read_directory)

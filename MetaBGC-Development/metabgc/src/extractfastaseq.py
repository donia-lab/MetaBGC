#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

from Bio import SeqIO
import os
import re
from multiprocessing import Pool, freeze_support
from itertools import repeat
import pandas as pd
import logging
"""
Function to parse fasta file based on text file with fasta header ids
"""
def ExtractFASTASeq(fasta_file,id_list,output_file):
    logging.info("Searching " + fasta_file + "...")
    record_dict = SeqIO.index(fasta_file, "fasta")
    records = []
    for readid in id_list:
        if readid in record_dict:
            seq_record = record_dict[readid]
            records.append(seq_record)
    count = SeqIO.write(records, output_file, "fasta")
    logging.info("Saved " + str(count) + " records from " + fasta_file + " to " + output_file)

"""
Function to run make search and extract FASTA file in parallel. 
"""
def RunExtractParallel(dbFileList,id_list,outFileList):
    numOfprocess = len(dbFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(ExtractFASTASeq, zip(dbFileList, repeat(id_list), outFileList))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Check if a given sample has any matched read at all.
"""
def sampleHasMatch(sample_list, sampleStr):
    for sample in sample_list:
        if sample in sampleStr:
            return True
    return False

"""
Function to run BLAST against a directory. 
"""
def RunExtractDirectoryPar(readsDir, readIDFile, ouputDir, outputFasta, fasta_file_ext, ncpus):
    try:
        df_reads = pd.read_csv(readIDFile, sep='\t')
        id_list = list(set(df_reads.readID.values.tolist()))
        sample_list = []
        for sample in list(set(df_reads.Sample.values.tolist())):
            sample_list.append(sample.split('-')[0])

        filePathDict={}
        for subdir, dirs, files in os.walk(readsDir):
            for file in files:
                filePath = os.path.join(subdir, file)
                match_str = r".*\."+ fasta_file_ext +"$"
                if re.match(match_str, file) and os.path.getsize(filePath) > 0:
                    filePathDict[filePath] = os.path.getsize(filePath)

        logging.info("Found " + str(len(filePathDict)) + " read files from which to extract.")

        dbFileList = []
        outFileList = []
        setCtr=0
        for filePath in sorted(filePathDict, key=filePathDict.get):
            file = os.path.basename(filePath)
            sampleStr = os.path.splitext(file)[0]
            if sampleHasMatch(sample_list, sampleStr):
                outputFileName = sampleStr + "." + fasta_file_ext
                outputFilePath = os.path.join(ouputDir, outputFileName)
                dbFileList.append(filePath)
                outFileList.append(outputFilePath)
            if len(dbFileList) >= ncpus:
                setCtr = setCtr + 1
                logging.info("Running " + str(setCtr) + " of " + str(int(len(filePathDict)/ncpus)) + " extract rounds.")
                RunExtractParallel(dbFileList, id_list, outFileList)
                dbFileList = []
                outFileList = []
        if len(dbFileList) > 0:
            logging.info("Running remaining " + str(len(dbFileList)) + " extracts.")
            RunExtractParallel(dbFileList, id_list, outFileList)

        with open(outputFasta, 'w') as outfile:
            for filename in os.listdir(ouputDir):
                if filename.endswith("."+fasta_file_ext):
                    filePath = os.path.join(ouputDir, filename)
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)
                else:
                    continue
    except:
        print("Metabgc-identify read extraction has failed. Please check your inputs and contact support on : https://github.com/donia-lab/MetaBGC")
        exit()

def RunExtractDescription(inputFasta, fasta_file_type):
    logging.info("Processing " + inputFasta + "...")
    record_dict={}
    for record in SeqIO.parse(inputFasta, fasta_file_type):
        record_dict[record.id] = record.description
    return record_dict
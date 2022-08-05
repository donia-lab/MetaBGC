#!/usr/bin/env python

#####################################################################################
# @author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

from Bio import SeqIO
import os
import re
import pandas as pd
import logging
from metabgc.src.producer_consumer import *

class ExtractTask:
    def __init__(self, fastaFile, id_list, output_file, exact_match=True):
        self.fastaFile = fastaFile
        self.id_list = id_list
        self.exact_match = exact_match
        self.output_file = output_file

"""
Function to parse fasta file based on text file with fasta header ids
"""
def ExtractFASTAConsumer(queue, lock):
    # Synchronize access to the console
    with lock:
        logging.info('Starting consumer => {}'.format(os.getpid()))
    # Run indefinitely
    while True:
        time.sleep(random.randint(0, 10))
        # If the queue is empty, queue.get() will block until the queue has data
        extract_task = queue.get()
        with lock:
            logging.info('{} got {}'.format(os.getpid(), extract_task.fastaFile))
        try:
            record_dict = SeqIO.index(extract_task.fastaFile, "fasta")
            records = []
            for readid in extract_task.id_list:
                if readid in record_dict and extract_task.exact_match:
                    seq_record = record_dict[readid]
                    records.append(seq_record)
                elif not extract_task.exact_match:
                    found_ids = [s for s in record_dict.keys() if s.find(readid) >= 0]
                    for ids in found_ids:
                        seq_record = record_dict[ids]
                        records.append(seq_record)
            count = SeqIO.write(records, extract_task.output_file, "fasta")
            record_dict.close()
            with lock:
                logging.info("Saved " + str(count) + " records from " + extract_task.fastaFile + " to " + extract_task.output_file)
        except Exception as e:
            with lock:
                logging.info("Failed to extract records from " + extract_task.fastaFile)
        queue.task_done()

"""
Check if a given sample has any matched read at all.
"""
def sampleHasMatch(sample_list, sampleStr):
    for sample in sample_list:
        if sample in sampleStr:
            return True
    return False

"""
Function to extract identified sequences from the sample FASTA files. 
"""
def RunExtractDirectoryPar(readsDir,
                           readIDFile,
                           ouputDir,
                           outputFasta,
                           fasta_file_ext,
                           match_exact=True,
                           check_sample_match=True,
                           ncpus=1):
    try:
        df_reads = pd.read_csv(readIDFile, sep='\t')
        id_list = list(set(df_reads.readID.values.tolist()))
        sample_list = []
        for sample in list(set(df_reads.Sample.values.tolist())):
            sample_list.append(sample.split('-')[0])

        filePathDict = {}
        for subdir, dirs, files in os.walk(readsDir):
            for file in files:
                filePath = os.path.join(subdir, file)
                match_str = r".*\." + fasta_file_ext + "$"
                if re.match(match_str, file) and os.path.getsize(filePath) > 0:
                    filePathDict[filePath] = os.path.getsize(filePath)

        logging.info("Found " + str(len(filePathDict)) + " read files from which to extract.")

        # Check and remove files that do not contribute a read to reduce the number of tasks
        extract_task_list = []
        for filePath in sorted(filePathDict, key=filePathDict.get):
            file = os.path.basename(filePath)
            sampleStr = os.path.splitext(file)[0]
            if sampleHasMatch(sample_list, sampleStr) and check_sample_match:
                outputFileName = sampleStr + "." + fasta_file_ext
                outputFilePath = os.path.join(ouputDir, outputFileName)
                extract_task = ExtractTask(filePath, id_list, outputFilePath, match_exact)
                extract_task_list.append(extract_task)
            else:
                outputFileName = sampleStr + "." + fasta_file_ext
                outputFilePath = os.path.join(ouputDir, outputFileName)
                extract_task = ExtractTask(filePath, id_list, outputFilePath, match_exact)
                extract_task_list.append(extract_task)

        logging.info('Number of pool processes:{0}.'.format(ncpus))
        print('Starting extract in: ' + str(len(extract_task_list)) + ' sample files.')
        # Create the Queue object
        queue = JoinableQueue()
        # Create a lock object to synchronize resource access
        lock = Lock()
        consumers = []

        print('Setting up consumers.')
        # Create consumer processes
        for i in range(ncpus):
            p = Process(target=ExtractFASTAConsumer, args=(queue, lock))
            consumers.append(p)
        # Start the producers and consumer
        # The Python VM will launch new independent processes for each Process object
        print('Starting up consumer sub-processes.')
        for c in consumers:
            c.start()

        print('Putting jobs in queue.')
        # Create our producer processes by passing the producer function and it's arguments
        for task_obj in extract_task_list:
            queue.put(task_obj)

        # join() method that synchronizes our program
        queue.join()
        print('Stopping sub-processes and terminating.')
        for c in consumers:
            c.terminate()
        print('Extraction of sequences completed. Starting concatenation...')
        logging.info('Extraction of sequences completed. Starting concatenation...')

        # Combining extracted reads
        with open(outputFasta, 'w') as outfile:
            for filename in os.listdir(ouputDir):
                if filename.endswith("." + fasta_file_ext):
                    filePath = os.path.join(ouputDir, filename)
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)
                else:
                    continue
    except:
        print(
            "Metabgc-identify read extraction has failed. Please check your inputs and contact support on : https://github.com/donia-lab/MetaBGC")
        exit()


def RunExtractDescription(inputFasta, fasta_file_type):
    logging.info("Processing " + inputFasta + "...")
    record_dict = {}
    for record in SeqIO.parse(inputFasta, fasta_file_type):
        record_dict[record.id] = record.description
    return record_dict

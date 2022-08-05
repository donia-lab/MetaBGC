#!/usr/bin/env python

#####################################################################################
# @author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
from metabgc.src.extractfastaseq import RunExtractDirectoryPar
from metabgc.src.hmmerrunlib import *
import os
import pandas as pd
from pathlib import Path
import logging


# Filter HMM results using predetermined spHMM models score cutoffs
def filter_spHMM_data(spHMM_df, cutoff_df):
    filter_spHMM_df = pd.DataFrame()
    for index, row in cutoff_df.iterrows():
        model_interval = cutoff_df.at[index, 'interval']
        cutoff = cutoff_df.at[index, 'cutoff']
        # filter spHMM dataframe with designated cutoffs
        filter_df = spHMM_df.query("HMMScore >= @cutoff & interval == @model_interval").reset_index(drop=True)
        filter_spHMM_df = filter_spHMM_df.append(filter_df)
    return (filter_spHMM_df)


def runidentify(hmm_file, cutoff_file, filteredTableFile, identifyReadIdFile):
    spHMM_df = pd.read_csv(hmm_file, sep="\t",
                           names=["readID", "sampleType", "Sample", "protType", "HMMScore", "window", "interval"])

    # Reformat IDs to remove frame identifier from transeq
    for i, row in spHMM_df.iterrows():
        translated_read_id = spHMM_df.at[i, 'readID']
        if translated_read_id.rfind('_') != -1:
            nucl_read_id = translated_read_id[:translated_read_id.rfind('_')]
            spHMM_df.at[i, 'readID'] = nucl_read_id

    cutoff_df = pd.read_csv(cutoff_file, sep="\t", header=0)
    spHMM_df_filtered = filter_spHMM_data(spHMM_df, cutoff_df)

    # Filter duplicate readIDs and keep the highest HMM Score
    spHMM_df_filtered_uniq = spHMM_df_filtered.groupby(['readID', 'Sample', 'sampleType', 'protType'])[
        'HMMScore'].max().reset_index()
    spHMM_df_filtered_uniq_interval = pd.merge(spHMM_df_filtered, spHMM_df_filtered_uniq, how='inner')

    spHMM_df_filtered_uniq_interval.to_csv(filteredTableFile, sep="\t", index=False)
    identifyReadIdList = list(set(spHMM_df_filtered_uniq_interval.readID.values.tolist()))

    with open(identifyReadIdFile, 'w') as outfile:
        for readID in identifyReadIdList:
            outfile.write(readID + '\n')


def mbgcidentify(sphmm_directory, cohort_name, nucl_seq_directory, prot_seq_directory,
                 seq_fmt, pair_fmt, r1_file_suffix, r2_file_suffix,
                 prot_family_name, hmm_search_output_directory, output_directory, cpu):
    try:
        if cpu is not None:
            CPU_THREADS = int(cpu)

        if r1_file_suffix is None:
            r1_file_suffix = ""
        if r2_file_suffix is None:
            r2_file_suffix = ""
        if prot_seq_directory is None:
            prot_seq_directory = ""

        if hmm_search_output_directory is None:
            hmm_search_output_directory = os.path.join(output_directory, 'hmm_identify_search')
        identify_directory = os.path.join(output_directory, 'identify_result')
        os.makedirs(identify_directory, 0o777, True)

        fasta_seq_dir = os.path.join(output_directory, 'fasta_seq_result')
        identifyReadIds = identify_directory + os.sep + "CombinedReadIds.txt"
        filteredHMMResult = identify_directory + os.sep + "spHMM-filtered-results.txt"
        allHMMResult = identify_directory + os.sep + "CombinedHmmSearch.txt"
        multiFastaFile = identify_directory + os.sep + "identified-biosynthetic-reads.fasta"

        nucl_seq_directory = PreProcessReadsPar(nucl_seq_directory,
                                                seq_fmt, pair_fmt,
                                                r1_file_suffix.strip(),
                                                r2_file_suffix.strip(),
                                                output_directory,
                                                CPU_THREADS)

        # Translate nucleotide seq
        if not os.path.isdir(prot_seq_directory):
            prot_seq_directory = TranseqReadsDir(output_directory, nucl_seq_directory, CPU_THREADS)

        # HMMER search
        if not os.path.exists(multiFastaFile):
            if not os.path.exists(allHMMResult):
                os.makedirs(hmm_search_output_directory, 0o777, True)
                for filename in os.listdir(sphmm_directory):
                    fileBase = Path(filename).resolve().stem
                    if filename.endswith(".hmm"):
                        hmmInterval = fileBase.split("__")[2]
                        hmmfilename = os.path.join(sphmm_directory, filename)
                        RunPCHMMDirectoryParallel(prot_seq_directory, hmmfilename, cohort_name, prot_family_name, "30_10",
                                                hmmInterval,
                                                hmm_search_output_directory, CPU_THREADS)
                found_hit_ctr = 0
                with open(allHMMResult, 'w') as outfile:
                    for subdir, dirs, files in os.walk(hmm_search_output_directory):
                        for file in files:
                            filePath = os.path.join(subdir, file)
                            if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                                with open(filePath) as infile:
                                    for line in infile:
                                        outfile.write(line)
                                        found_hit_ctr = found_hit_ctr + 1
                if found_hit_ctr == 0:
                    logging.info("Metabgc-identify has failed to identify any reads for this protein family model. "
                                 "Please try with a different metagenome.")
                    print("Metabgc-identify has failed to identify any reads for this protein family model. "
                        "Please try with a different metagenome.")
                    exit()
                print("Metabgc-identify HMMER search is complete.")
            else:
                print("Metabgc-identify is using the existing HMMER search result found.")

            cutoff_file = os.path.join(sphmm_directory, prot_family_name + "_F1_Cutoff.tsv")
            if not os.path.exists(cutoff_file):
                print("Metabgc-identify did not find the spHMM cutoff threshold file at the expected location.:" + cutoff_file)
                raise

            ##Run identify thresholding
            runidentify(allHMMResult, cutoff_file, filteredHMMResult, identifyReadIds)

            os.makedirs(fasta_seq_dir, 0o777, True)
            RunExtractDirectoryPar(nucl_seq_directory, filteredHMMResult,
                                   fasta_seq_dir, multiFastaFile, "fasta", True, True,
                                   CPU_THREADS)
        else:
            print("Metabgc-identify is returning the existing identified reads found.")
        return multiFastaFile
    except:
        print("Metabgc-identify has failed. Please check your inputs and contact support on : https://github.com/donia-lab/MetaBGC")
        exit()

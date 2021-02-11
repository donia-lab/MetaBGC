
import pandas as pd
from time import time
from Bio import SeqIO
import numpy as np
import json
import re
import os


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
### Remove the reading frame information from the readID.
# Function to aggregate identitical sample reads located at different reading frames
# and take the frame with the highest HMM score
def formatHMM(hmm_df):
    hmm_df[['readIDOnly','F_R_read_frame']] = hmm_df.readID.str.split('/',n=2,expand=True)
    hmm_df[['F_R','frameNumb']] = hmm_df.F_R_read_frame.str.split('_',n=2,expand=True)
    hmm_df['readID'] = hmm_df['readIDOnly'] + '/' + hmm_df['F_R']
    hmm_df = hmm_df.groupby(["readID", "sampleID", "sampleType", "protType", "window", "interval"]).agg({'HMMScore': ['max']}).reset_index()
    hmm_df.columns=["readID","Sample", "sampleType", "protType", "window", "interval", "HMMScore"]
    return(hmm_df)

### Positional read analysis in respect to location mapped to siderophore domain
# Keep reads that map to the interval of a given model and covers the model 90%
def filter_blast(all_blast_df, gene_positions):
    df_op = pd.DataFrame()
    for index, gene_data in gene_positions.iterrows():
        print("Completed " + gene_data['gene_name'] + " " + gene_data['start'])
        #filters datatframe for edges and internal reads compared to model interval

        interval_df_1 = all_blast_df[all_blast_df['qseqid'] == gene_data['gene_name']]
        interval_df_1 = interval_df_1[(interval_df_1['qstart']>=gene_data['start']) & (interval_df_1['qstart']<=gene_data['end'])]
        interval_df_1 = interval_df_1[(interval_df_1['qend']>=gene_data['start']) & (interval_df_1['qend']<=gene_data['end'])]

        interval_df_2 = all_blast_df[all_blast_df['qseqid'] == gene_data['gene_name']]
        interval_df_2 = interval_df_2[interval_df_2['qstart'] < gene_data['start']]
        interval_df_2 = interval_df_2[interval_df_2['qend'] > gene_data['end']]

        interval_df = pd.concat([interval_df_1,interval_df_2])

        if len(interval_df) > 0:
            gene_len = abs(gene_data['end'] - gene_data['start']) + 1
            interval_df['model_cov'] = (getOverlap([interval_df['qstart'], interval_df['qend']],[gene_data['start'], gene_data['end']]) / gene_len) * 100.0
            interval_df['interval'] = gene_data['interval']
            interval_df <- interval_df[interval_df['model_cov'] >= 90]
            df_op <- pd.concat([df_op, interval_df])
    return df_op

### Determine cutoffs for each interval within the model
def compare_reads(hmm_df, blast_df):
    blast_df.rename(columns={"sseqid": "readID"}, inplace=True)
    # remove columns to compare the two dataframe
    blastDF = blast_df.drop(['model_cov', 'interval', 'Sample', 'sampleType'])
    common_reads = hmm_df.merge(blastDF)
    common_reads['readCheck'] = "common-read"
    hmm_unique_reads = hmm_df.merge(blastDF, how='outer', indicator=True)
    hmm_unique_reads = hmm_unique_reads[hmm_unique_reads['_merge'] == 'left_only']
    hmm_unique_reads['readCheck'] = "hmm-unique-read"
    compared_data = pd.concat([common_reads, hmm_unique_reads])
    return(compared_data)

def evaluate_sphmms(HMMRunFile,
                    BLAST_TP_NoCov_File,
                    GeneIntervalPosFile,
                    HMM_Model_Name,
                    F1_Threshold,
                    HMMOutDir,
                    HMMHighPerfOutDir):

    ### Load segmented profiled HMMs for synthetic genomes
    #load HMM data and recode sampleType for the complexity of the synthetic sample (#of genomes in samples)
    hmm_df = pd.read_csv(HMMRunFile, delimiter="\t", header=None)
    hmm_df.columns=["readID", "sampleType", "sampleID", "protType", "HMMScore", "window","interval"]
    #Keep duplicated reads if they are in different reads
    hmm_df_recoded = formatHMM(hmm_df)

    ### Load the BLAST data for genes against our synthetic dataset.
    # BLAST unfiltered reads at 95% pident no readCoverage filter
    all_blast_df = pd.read_csv(BLAST_TP_NoCov_File, delimiter="\t")

    ### Positional information about domains and their locations in respect to the 30_10 spHMM models
    gene_positions = pd.read_csv(GeneIntervalPosFile, delimiter="\t")

    #Filter BLAST reads that are within the genes intervals and cover 90% of the model interval
    start_time = time()
    blast_intervals = filter_blast(all_blast_df, gene_positions)
    end_time = time()
    print(f'Genes interval finding took {round(end_time - start_time, 4)} seconds.')

    blast_bin = compare_reads(hmm_df_recoded, blast_intervals)

    df_hmm_cutoff_scores = blast_bin[blast_bin['readCheck'] == "common-read"]
    df_hmm_cutoff_scores = df_hmm_cutoff_scores.groupby(['interval'])['HMMScore'].median().reset_index(name='medianScore')
    df_hmm_cutoff_scores['medianScore'] = round(df_hmm_cutoff_scores['medianScore'])
    df_hmm_cutoff_scores.drop_duplicates(inplace=True)
    df_hmm_cutoff_scores.columns = ["interval", "read_check", "cutoff"]

    #Filter data with cutoffs to compare to BLAST interval reads
    filtered_median = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"])
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_median = pd.concat(filtered_median,df_found)

    filtered_subfive = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"]) - 5.0
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_subfive = pd.concat(filtered_subfive,df_found)

    filtered_plusfive = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"]) + 5.0
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_plusfive = pd.concat(filtered_plusfive,df_found)



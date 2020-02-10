#!/usr/bin/env python

#####################################################################################
#@author: shuowang
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

import pandas as pd
from sklearn.cluster import DBSCAN
from collections import Counter
from statistics import mean
from Bio import SeqIO
import numpy as np
import json
import re
import os

def PrintBinSeqs(binIds,df_read_labels,identifiedReadFile,outDir):
    seq_dict = SeqIO.index(identifiedReadFile, "fasta")
    single_records = []
    fastaDirGT10 = outDir + "/bin_fasta/gt10"
    fastaDirRem = outDir + "/bin_fasta/rem"
    os.makedirs(fastaDirGT10, 0o777, True)
    os.makedirs(fastaDirRem, 0o777, True)
    for bin in binIds:
        df_bin_reads = df_read_labels[df_read_labels['bin']==bin]
        records = []
        count_row = df_bin_reads.shape[0]
        if count_row >= 10:
            for index, row in df_bin_reads.iterrows():
                seq_record = seq_dict[row['qseqid']]
                records.append(seq_record)
            output_file = os.path.join(fastaDirGT10,str(bin)+".fasta")
            SeqIO.write(records, output_file, "fasta")
        elif count_row >= 2:
            for index, row in df_bin_reads.iterrows():
                seq_record = seq_dict[row['qseqid']]
                records.append(seq_record)
            output_file = os.path.join(fastaDirRem,str(bin)+".fasta")
            SeqIO.write(records, output_file, "fasta")
        else:
            seq_record = seq_dict[df_bin_reads['qseqid'].iloc[0]]
            single_records.append(seq_record)
    output_file = os.path.join(fastaDirRem, "singleton.fasta")
    SeqIO.write(single_records, output_file, "fasta")

def mbgccluster(table, tableAbundance, identifiedReadFile,max_dist,min_samples, readThresh, abundThresh, cpu):
    df = pd.read_csv(table, delimiter="\t")
    df_abundance = pd.read_csv(tableAbundance, delimiter="\t")
    mat = df.iloc[:,1:df.shape[1]].values
    read_names = df.iloc[:,0].values
    cl = DBSCAN(eps=max_dist, min_samples=min_samples, metric="correlation",
                algorithm="brute", n_jobs=cpu).fit_predict(mat)

    read_labels = dict(zip(read_names, cl.tolist()))
    out_file_json = re.sub("(.*)\\..*", r"\1_DBSCAN.json", table)
    dir_path = os.path.dirname(os.path.realpath(table))
    out_file_abund = os.path.join(dir_path,"ReadLevelAbundance.tsv")
    out_file_summary = os.path.join(dir_path, "BinSummary.txt")
    out_file_abund_sample = os.path.join(dir_path, "SampleAbundanceMatrix.tsv")
    with open(out_file_json, "w") as h:
        json.dump(read_labels, h, indent=4)

    clusterLabels = cl.tolist()
    noiseCluster = 0
    if -1 in clusterLabels:
        noiseCluster = 1
    clusterIDs = Counter(clusterLabels).keys()
    clusterCtr = len(Counter(clusterLabels).keys()) - noiseCluster

    clusterFreq = Counter(clusterLabels).values()
    clusterFreqGTThresh = [i for i in clusterFreq if i >= readThresh]
    clusterIdsGTThresh = [i for i,j in zip(clusterIDs,clusterFreq) if j >= readThresh]

    outF = open(out_file_summary, "w")
    outF.write("Number of Reads Quantified: {0}\n".format(len(clusterLabels)))
    outF.write("Number of Samples Quantified: {0}\n".format(len(Counter(df_abundance['Sample'].tolist()).keys())))
    outF.write("Number of Bins: {0}\n".format(clusterCtr))
    outF.write("Number of Bins with >= {0} Reads: {1}\n".format(readThresh,len(clusterFreqGTThresh)))
    outF.write("Number of Bins with >= 5 Reads: {0}\n".format(len([i for i in clusterFreq if i >= 5])))
    outF.write("Number of Bins with 2-4 Reads: {0}\n".format(len([i for i in clusterFreq if i >= 2 and i<=4])))
    outF.write("Number of Singleton Read Bins: {0}\n".format(len([i for i in clusterFreq if i == 1])))

    df_read_labels = pd.DataFrame(list(read_labels.items()),columns=['qseqid', 'bin'])
    df_abundance = pd.merge(df_abundance, df_read_labels, how='inner')
    PrintBinSeqs(clusterIDs,df_read_labels,identifiedReadFile,dir_path)

    df_abundance = df_abundance[df_abundance['bin'].isin(clusterIdsGTThresh)]
    df_abundance.rename(columns={'count': 'ReadAbundance'}, inplace=True)

    df_abundance_sample = df_abundance.groupby(['Sample', 'bin'])['ReadAbundance'].sum().reset_index()
    df_abundance_sample.rename(columns={'ReadAbundance': 'BinAbundance'}, inplace=True)

    df_bins = df_abundance_sample.groupby(['bin'])['BinAbundance'].sum().reset_index()

    outF.write("Average Bin Abundance of Bins with >= {0} Reads: {1}\n".format(readThresh,round(mean(df_bins['BinAbundance'].tolist()),2)))
    outF.write("Maximum Bin Abundance of Bins with >= {0} Reads: {1}\n".format(readThresh,round(max(df_bins['BinAbundance'].tolist()),2)))
    outF.write("Minimum Bin Abundance of Bins with >= {0} Reads: {1}\n".format(readThresh,round(min(df_bins['BinAbundance'].tolist()),2)))

    df_abundance_sample = df_abundance_sample[df_abundance_sample['BinAbundance']>=abundThresh]
    df_abundance_sample_pivot = df_abundance_sample.pivot_table(index='Sample', columns='bin',values='BinAbundance', fill_value=0)
    df_abundance_sample_pivot.to_csv(out_file_abund_sample,sep='\t')

    df_abundance = pd.merge(df_abundance, df_abundance_sample, how='inner')
    df_abundance['log10BinAdbundance'] = np.log10(df_abundance['BinAbundance'])
    df_abundance.to_csv(out_file_abund, index=False,sep='\t')

    outF.write("Average Number of Reads in Bins with >= {0} Reads: {1}\n".format(readThresh,round(mean(clusterFreqGTThresh),2)))
    outF.write("Average Abundance of Bins with >= {0} Reads: {1}\n".format(readThresh,round(mean(clusterFreqGTThresh),2)))
    outF.write("Number Samples in containing Bins with >= {0} Reads at Bin Abundance: {1}\n".format(readThresh, len(Counter(df_abundance_sample['Sample'].tolist()).keys())))
    outF.close()

    return out_file_json

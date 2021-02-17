import shutil
import pandas as pd
from time import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def getOverlap(row, b):
    return max(0, min(row['qend'], b[1]) - max(row['qstart'], b[0])) + 1
### Remove the reading frame information from the readID.
# Function to aggregate identitical sample reads located at different reading frames
# and take the frame with the highest HMM score
def formatHMM(hmm_df):
    hmm_df[['readIDOnly','F_R_read_frame']] = hmm_df.readID.str.split('/',n=2,expand=True)
    hmm_df[['F_R','frameNumb']] = hmm_df.F_R_read_frame.str.split('_',n=2,expand=True)
    hmm_df.loc[:,'readID'] = hmm_df['readIDOnly'] + '/' + hmm_df['F_R']
    hmm_df = hmm_df.groupby(["readID", "sampleID", "sampleType", "protType", "window", "interval"]).agg({'HMMScore': ['max']}).reset_index()
    hmm_df.columns=["readID","Sample", "sampleType", "protType", "window", "interval", "HMMScore"]
    return(hmm_df)

### Positional read analysis in respect to location mapped to siderophore domain
# Keep reads that map to the interval of a given model and covers the model 90%
def filter_blast(all_blast_df, gene_positions):
    df_op = pd.DataFrame()
    for index, gene_data in gene_positions.iterrows():
        print("Completed " + gene_data['gene_name'] + " " + str(gene_data['start']))
        #filters datatframe for edges and internal reads compared to model interval

        interval_df_1 = all_blast_df[all_blast_df['qseqid'] == gene_data['gene_name']]
        interval_df_1 = interval_df_1[((interval_df_1['qstart']>=gene_data['start']) & (interval_df_1['qstart']<=gene_data['end'])) |
                                      ((interval_df_1['qend']>=gene_data['start']) & (interval_df_1['qend']<=gene_data['end']))]

        interval_df_2 = all_blast_df[all_blast_df['qseqid'] == gene_data['gene_name']]
        interval_df_2 = interval_df_2[interval_df_2['qstart'] < gene_data['start']]
        interval_df_2 = interval_df_2[interval_df_2['qend'] > gene_data['end']]

        interval_df = pd.concat([interval_df_1,interval_df_2])

        if len(interval_df) > 0:
            gene_len = abs(gene_data['end'] - gene_data['start']) + 1
            interval_df['model_cov'] = interval_df.apply (lambda row: getOverlap(row, [gene_data['start'], gene_data['end']]), axis=1)
            interval_df['model_cov'] = (interval_df['model_cov'] / gene_len) * 100.0
            interval_df['interval'] = gene_data['interval']
            interval_df = interval_df[interval_df['model_cov'] >= 90]
            df_op = pd.concat([df_op, interval_df])
    return df_op

### Determine cutoffs for each interval within the model
def compare_reads(hmm_df, blast_df):
    blast_df_temp = blast_df.rename(columns={"sseqid": "readID"})
    # remove columns to compare the two dataframe
    blastDF = blast_df_temp.drop(columns=['model_cov', 'interval', 'Sample', 'sampleType'])
    common_reads = hmm_df.merge(blastDF)
    in_both = hmm_df['readID'].isin(common_reads['readID'])
    common_reads = hmm_df[in_both]
    common_reads.loc[:,'readCheck'] = "common-read"
    hmm_unique_reads = hmm_df.merge(blastDF, how='outer', indicator=True)
    hmm_unique_reads = hmm_unique_reads[hmm_unique_reads['_merge'] == 'left_only']
    hmm_unique_reads = hmm_unique_reads[hmm_df.columns]
    hmm_unique_reads.loc[:,'readCheck'] = "hmm-unique-read"
    compared_data = pd.concat([common_reads, hmm_unique_reads])
    return compared_data

#return the hmm-unique reads
def compare_hmm_unique(hmm_df, blast_df, pos_df):
    blast_df_renamed = blast_df.rename(columns={"sseqid": "readID"})
    hmm_unique_df = hmm_df.merge(blast_df_renamed, on = ["readID", "Sample", "sampleType", "protType"])
    intervals = hmm_unique_df.interval.unique()
    results = pd.DataFrame()
    #check that reads are in the same interval to throw out
    for interval in intervals:
        gene_interval_data = pos_df[pos_df['interval'] == interval]
        if len(gene_interval_data) > 0:
            for index, gene_data in gene_interval_data.iterrows():
                interval_df = hmm_unique_df[hmm_unique_df['qseqid'] == gene_data['gene_name']]
                interval_df = interval_df[((interval_df['qstart']>=gene_data['start']) & (interval_df['qstart']<=gene_data['end']))|
                                          ((interval_df['qend']>=gene_data['start']) & (interval_df['qend']<=gene_data['end']))]
                if len(interval_df) > 0:
                    results = pd.concat([results,interval_df])
    return results.drop_duplicates()

def return_hmm_unique(hmm_df, blast_df):
    temp_hmm_df = hmm_df.drop(columns=['window'])
    temp_hmm_df.loc[:,'interval'] = temp_hmm_df['interval'].apply(str)
    blast_df_temp = blast_df.rename(columns={"sseqid": "readID"})
    temp_hmm_unique = temp_hmm_df.merge(blast_df_temp, how='outer', on=["readID", "Sample", "sampleType", "interval"], indicator=True)
    temp_hmm_unique = temp_hmm_unique[temp_hmm_unique['_merge'] == 'left_only']
    temp_hmm_unique = temp_hmm_unique[temp_hmm_df.columns]
    return temp_hmm_unique


def calculate_F1(hmm_df, hmm_fp, blast_df, intervals):
    #added this because factor vector
    hmm_df.loc[:,'interval'] = hmm_df['interval'].apply(str)
    temp_hmm_df = hmm_df.drop(columns=['window'])
    blast_df_renamed = blast_df.rename(columns={"sseqid": "readID"})
    results = pd.DataFrame({'interval': pd.Series([], dtype='str'),
                            'F1': pd.Series([], dtype='float'),
                            'TP': pd.Series([], dtype='float'),
                            'FP': pd.Series([], dtype='float'),
                            'FN': pd.Series([], dtype='float')})
    for interval in intervals:
        model_hmm = hmm_df[hmm_df['interval'] == interval] # filter out model HMM results
        model_blast = blast_df_renamed[blast_df_renamed['interval'] == interval]
        model_hmm_fp = hmm_fp[hmm_fp['interval'] == interval]
        TP = len(model_hmm.merge(model_blast))
        FP = len(model_hmm_fp)
        temp_fn = model_blast.merge(model_hmm, how='outer', indicator=True)
        temp_fn = temp_fn[temp_fn['_merge'] == 'left_only']
        temp_fn = temp_fn[model_blast.columns]
        FN = len(temp_fn)
        F1_metric = (2*TP)/((2*TP)+FP+FN)
        row = [interval, F1_metric, TP, FP, FN]
        results.loc[len(results)] = row
    return results


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
    df_hmm_cutoff_scores = df_hmm_cutoff_scores.groupby(['interval', 'readCheck'])['HMMScore'].median().reset_index(name='medianScore')
    df_hmm_cutoff_scores['medianScore'] = round(df_hmm_cutoff_scores['medianScore'])
    df_hmm_cutoff_scores.drop_duplicates(inplace=True)
    df_hmm_cutoff_scores.rename(columns={"readCheck": "read_check", "medianScore": "cutoff"}, inplace=True)

    #Filter data with cutoffs to compare to BLAST interval reads
    filtered_median = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"])
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_median = pd.concat([filtered_median, df_found])

    filtered_subfive = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"]) - 5.0
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_subfive = pd.concat([filtered_subfive, df_found])

    filtered_plusfive = pd.DataFrame(columns=hmm_df_recoded.columns)
    for index, row in df_hmm_cutoff_scores.iterrows():
        intervalStr = str(row["interval"])
        cutoffScore = float(row["cutoff"]) + 5.0
        df_found = hmm_df_recoded[(hmm_df_recoded['interval'] == intervalStr) & (hmm_df_recoded['HMMScore'] >= cutoffScore)]
        filtered_plusfive = pd.concat([filtered_plusfive, df_found])


    gene_positions.loc[:,'protType']  = HMM_Model_Name
    all_blast_df.loc[:,'protType'] = HMM_Model_Name

    median_hmmunique = return_hmm_unique(filtered_median, blast_intervals)
    median_hmmunique_less_model_cov = compare_hmm_unique(median_hmmunique, all_blast_df, gene_positions)
    median_remaining_hmm = median_hmmunique.merge(median_hmmunique_less_model_cov, how='outer', indicator=True)
    median_remaining_hmm = median_remaining_hmm[median_remaining_hmm['_merge'] == 'left_only']
    median_remaining_hmm = median_remaining_hmm[median_hmmunique.columns]

    subfive_hmmunique = return_hmm_unique(filtered_subfive, blast_intervals)
    subfive_hmmunique_less_model_cov = compare_hmm_unique(subfive_hmmunique, all_blast_df, gene_positions)
    subfive_remaining_hmm = subfive_hmmunique.merge(subfive_hmmunique_less_model_cov, how='outer', indicator=True)
    subfive_remaining_hmm = subfive_remaining_hmm[subfive_remaining_hmm['_merge'] == 'left_only']
    subfive_remaining_hmm = subfive_remaining_hmm[subfive_hmmunique.columns]

    plusfive_hmmunique = return_hmm_unique(filtered_plusfive, blast_intervals)
    plusfive_hmmunique_less_model_cov = compare_hmm_unique(plusfive_hmmunique, all_blast_df, gene_positions)
    plusfive_remaining_hmm = plusfive_hmmunique.merge(plusfive_hmmunique_less_model_cov, how='outer', indicator=True)
    plusfive_remaining_hmm = plusfive_remaining_hmm[plusfive_remaining_hmm['_merge'] == 'left_only']
    plusfive_remaining_hmm = plusfive_remaining_hmm[plusfive_hmmunique.columns]

    f1_cutoff_median = calculate_F1(filtered_median, median_remaining_hmm, blast_intervals, hmm_df_recoded.interval.unique())
    f1_cutoff_median.loc[:,'cutoff'] = "median"
    FP_Reads_Fname = HMM_Model_Name + "-FP_Reads.csv"
    median_remaining_hmm.to_csv(os.path.join(HMMHighPerfOutDir,FP_Reads_Fname), sep='\t',index=False)

    f1_cutoff_plusfive = calculate_F1(filtered_plusfive, plusfive_remaining_hmm, blast_intervals, hmm_df_recoded.interval.unique())
    f1_cutoff_plusfive.loc[:,'cutoff'] = "+5"

    f1_cutoff_subfive = calculate_F1(filtered_subfive, subfive_remaining_hmm, blast_intervals, hmm_df_recoded.interval.unique())
    f1_cutoff_subfive.loc[:,'cutoff'] = "-5"

    f1_df = pd.concat([f1_cutoff_median, f1_cutoff_plusfive, f1_cutoff_subfive])
    f1_cutoff_df = f1_df.groupby(["interval"])['F1'].max().reset_index(name='F1_Max')
    f1_cutoff_df = f1_cutoff_df[f1_cutoff_df['F1_Max']>= F1_Threshold]
    f1_cutoff_df = f1_cutoff_df.rename(columns={"F1_Max": "F1"})
    f1_cutoff_df = f1_cutoff_df.merge(f1_df)

    f1_cutoff_df.rename(columns={"cutoff": "cutoff_diff"}, inplace=True)

    f1_cutoff_df.loc[f1_cutoff_df['cutoff_diff'] == 'median', 'cutoff_diff'] = 0
    f1_cutoff_df  = f1_cutoff_df.merge(df_hmm_cutoff_scores,on=["interval"])
    f1_cutoff_df = f1_cutoff_df.rename(columns={"cutoff_diff": "FinalCutoff"})
    cutoffFileName = HMM_Model_Name + "_F1_Cutoff.tsv"
    scoringFileName = HMM_Model_Name + "_Scores.tsv"
    f1_cutoff_df.to_csv(os.path.join(HMMHighPerfOutDir,cutoffFileName), sep='\t',index=False)
    f1_df.to_csv(os.path.join(HMMHighPerfOutDir,scoringFileName), sep='\t',index=False)

    for index, row in f1_cutoff_df.iterrows():
        spHMMInterval = row['interval']
        spHMMFileName = HMM_Model_Name + "__30_10__" + spHMMInterval + ".hmm"
        try:
            shutil.copy(os.path.join(HMMOutDir,spHMMFileName), HMMHighPerfOutDir)
        except IOError as e:
            print("Unable to copy file. %s" % e)
        except:
            print("Unexpected error:", sys.exc_info())

    # Setup dataframes for plotting.
    f1_cutoff_median[['interval_start','interval_end']] = f1_cutoff_median.interval.str.split("_",expand=True)
    f1_cutoff_median = f1_cutoff_median.astype({'interval_start': 'int32', 'interval_end': 'int32'})
    f1_cutoff_median.sort_values(by=['interval_start', 'interval_end'], inplace=True)

    f1_cutoff_plusfive[['interval_start','interval_end']] = f1_cutoff_plusfive.interval.str.split("_",expand=True)
    f1_cutoff_plusfive = f1_cutoff_plusfive.astype({'interval_start': 'int32', 'interval_end': 'int32'})
    f1_cutoff_plusfive.sort_values(by=['interval_start', 'interval_end'], inplace=True)

    f1_cutoff_subfive[['interval_start','interval_end']] = f1_cutoff_subfive.interval.str.split("_",expand=True)
    f1_cutoff_subfive = f1_cutoff_subfive.astype({'interval_start': 'int32', 'interval_end': 'int32'})
    f1_cutoff_subfive.sort_values(by=['interval_start', 'interval_end'], inplace=True)

    fig = plt.figure(figsize=(25, 6), dpi = 600.0)
    plt.title('HMM Score Cutoff')
    plt.xlabel('Interval')
    plt.ylabel('F1')
    plt.title('')
    line1 = plt.plot(f1_cutoff_median['interval'], f1_cutoff_median['F1'], marker='o', linewidth=1, markersize=2, color='seagreen', label='Median')
    line2 = plt.plot(f1_cutoff_plusfive['interval'], f1_cutoff_plusfive['F1'], marker='o', linewidth=1, markersize=2, color='cyan', label='+5')
    line3 = plt.plot(f1_cutoff_subfive['interval'], f1_cutoff_subfive['F1'], marker='o', linewidth=1, markersize=2, color='red', label='-5')
    plt.hlines(float(F1_Threshold), 0, len(f1_cutoff_plusfive), linestyles='dashed', color='lime')
    plt.legend(loc="upper right")
    F1_Plot_Fname  = HMM_Model_Name + "-F1_Plot.png"
    plt.savefig(os.path.join(HMMHighPerfOutDir,F1_Plot_Fname))

    fig = plt.figure(figsize=(25, 6), dpi = 600.0)
    plt.title('HMM Score Cutoff')
    plt.xlabel('Interval')
    plt.ylabel('TP')
    plt.title('')
    line1 = plt.plot(f1_cutoff_median['interval'], f1_cutoff_median['TP'], marker='o', linewidth=1, markersize=2, color='seagreen', label='Median')
    line2 = plt.plot(f1_cutoff_plusfive['interval'], f1_cutoff_plusfive['TP'], marker='o', linewidth=1, markersize=2, color='cyan', label='+5')
    line3 = plt.plot(f1_cutoff_subfive['interval'], f1_cutoff_subfive['TP'], marker='o', linewidth=1, markersize=2, color='red', label='-5')
    plt.legend(loc="upper right")
    F1_Plot_Fname  = HMM_Model_Name + "-TP_Plot.png"
    plt.savefig(os.path.join(HMMHighPerfOutDir,F1_Plot_Fname))

    fig = plt.figure(figsize=(25, 6), dpi = 600.0)
    plt.title('HMM Score Cutoff')
    plt.xlabel('Interval')
    plt.ylabel('FP')
    plt.title('')
    line1 = plt.plot(f1_cutoff_median['interval'], f1_cutoff_median['FP'], marker='o', linewidth=1, markersize=2, color='seagreen', label='Median')
    line2 = plt.plot(f1_cutoff_plusfive['interval'], f1_cutoff_plusfive['FP'], marker='o', linewidth=1, markersize=2, color='cyan', label='+5')
    line3 = plt.plot(f1_cutoff_subfive['interval'], f1_cutoff_subfive['FP'], marker='o', linewidth=1, markersize=2, color='red', label='-5')
    plt.legend(loc="upper right")
    F1_Plot_Fname  = HMM_Model_Name + "-FP_Plot.png"
    plt.savefig(os.path.join(HMMHighPerfOutDir,F1_Plot_Fname))






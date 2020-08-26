import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import math
from metabgc.src.utils import *
from metabgc.src.blastrunlib import *

def ReadLevelBinAnalytics(readTableAbundance,
                          sampleTableAbundance,
                          spHMMFile,
                          blastResultFile,
                          fastaFile,
                          cohortMetadataFile,
                          assembly_metadata_file,
                          output_dir):
    df_read_abundance = pd.read_csv(readTableAbundance, delimiter="\t")
    df_sample_abundance = pd.read_csv(sampleTableAbundance, delimiter="\t")
    df_spHMM = pd.read_csv(spHMMFile, delimiter="\t")

    df_metadata = pd.read_csv(cohortMetadataFile, delimiter=",", keep_default_na=False)
    df_metadata = df_metadata.drop(columns=['Subject','BodyAggSite','Subject_Status','Visits'])
    df_metadata.rename(columns={'Sample': 'derived_sample', 'Cohort': 'derived_cohort','Bodysite': 'derived_bodysite'}, inplace=True)

    df_BLAST = pd.read_csv(blastResultFile, delimiter="\t",header=None,names=['qseqid','sseqid','identity','e-value','staxids','genome','scomnames','sskingdoms','stitle'])
    df_BLAST = df_BLAST.drop(columns=['sseqid','staxids','scomnames','sskingdoms'])

    record_dict = SeqIO.index(fastaFile, "fasta")
    df_fasta = pd.DataFrame(columns=['qseqid', 'read_seq'])
    for read_id in record_dict:
        seq_str = str(record_dict[read_id].seq)
        df_fasta = df_fasta.append({'qseqid': read_id,
                                    'read_seq': seq_str}, ignore_index=True)

    # Read stats
    df_abundance_max = df_read_abundance.groupby(['qseqid','bin'])['ReadAbundance'].max().reset_index(name='max_sample_abundance_read')
    df_abundance_total = df_read_abundance.groupby(['qseqid'])['ReadAbundance'].sum().reset_index(name='total_abundance_read')
    df_abundance_median = df_read_abundance.groupby(['qseqid'])['ReadAbundance'].median().reset_index(name='median_abundance_read')
    df_abundance_count = df_read_abundance.groupby(['qseqid'])['Sample'].count().reset_index(name='n_sample_read')

    df_read_abundance_count =  df_read_abundance.groupby(['Sample','qseqid'])['ReadAbundance'].max().reset_index()
    df_read_abundance_count.rename(columns={'ReadAbundance': 'max_sample_abundance_read'}, inplace=True)

    df_abundance_max_median = pd.merge(df_abundance_max,df_abundance_median,on=['qseqid'],how='inner')
    df_read_max_sample = pd.merge(df_abundance_max_median,df_read_abundance_count,on=['qseqid','max_sample_abundance_read'],how='inner')
    df_read_max_sample.rename(columns={'Sample': 'max_sample_read'}, inplace=True)
    df_read_max_sample = pd.merge(df_read_max_sample, df_abundance_total, on=['qseqid'], how='inner')
    df_read = pd.merge(df_read_max_sample, df_abundance_count, on=['qseqid'], how='inner')

    # Bin stats
    sampleList = []
    df_sample = pd.DataFrame(columns=['bin', 'total_abundance_bin','n_sample_bin', 'max_sample_id_bin', 'median_sample_abundance_bin', 'max_sample_abundance_bin'])
    for (columnName, columnData) in df_sample_abundance.iteritems():
        if columnName == 'Sample':
            sampleList = columnData.values
            continue
        else:
            binID = columnName
            total_abundance_bin = sum(columnData.values)
            values_nz = np.argwhere(columnData.values)
            n_sample_bin= len(values_nz)
            median_sample_abundance_bin = np.median(values_nz)
            max_sample_abundance_bin = max(columnData.values)
            max_sample_id_bin = np.argmax(columnData.values)
            max_sample_bin = sampleList[max_sample_id_bin]
            df_sample = df_sample.append({'bin': int(binID),
                              'total_abundance_bin': total_abundance_bin,
                              'n_sample_bin': n_sample_bin,
                              'max_sample_bin': max_sample_bin,
                              'median_sample_abundance_bin': median_sample_abundance_bin,
                              'max_sample_abundance_bin': max_sample_abundance_bin}, ignore_index=True)
    # Read Bin Stats
    df_read_bin = pd.merge(df_read, df_sample,on=['bin'], how='inner')
    df_spHMM = df_spHMM[['readID','sampleType','Sample','HMMScore','interval']]
    df_spHMM = df_spHMM.drop(columns={'sampleType'})
    df_spHMM.rename(columns={'readID':'qseqid','Sample': 'derived_sample','HMMScore':'hmm_score'},inplace=True)

    df_summary = pd.merge(df_read_bin, df_spHMM,on=['qseqid'], how='inner')

    df_summary = pd.merge(df_summary, df_metadata, on=['derived_sample'], how='left')
    df_summary = pd.merge(df_summary, df_BLAST, on=['qseqid'], how='left')
    df_summary = pd.merge(df_summary, df_fasta, on=['qseqid'], how='left')

    df_summary = df_summary[['qseqid',
                            'derived_sample',
                            'derived_cohort',
                            'derived_bodysite',
                            'total_abundance_read',
                            'n_sample_read',
                            'max_sample_read',
                            'max_sample_abundance_read',
                            'median_abundance_read',
                            'bin',
                            'total_abundance_bin',
                            'n_sample_bin',
                            'max_sample_bin',
                            'max_sample_abundance_bin',
                            'median_sample_abundance_bin',
                            'hmm_score',
                            'interval',
                            'genome',
                            'identity',
                            'e-value',
                            'stitle',
                            'read_seq']]

    out_file_abund = os.path.join(output_dir, 'ReadLevelBinStats.tsv')
    df_summary.to_csv(out_file_abund, index=False, sep='\t')

    out_bin_max_sample_file = os.path.join(output_dir, 'MaxBinSample.tsv')
    df_sample.to_csv(out_bin_max_sample_file, index=False, sep='\t')


    ## For implementation of scaffold extraction
    df_assembly_metadata = pd.read_csv(assembly_metadata_file, delimiter=",", header=None, names=['Sample','Sample_Path'])
    df_assembly_metadata.set_index('Sample',inplace=True)
    f_asm = open(os.path.join(output_dir,"max_sample_assemblies.txt"), "w")
    f_bin = open(os.path.join(output_dir,"bin_paths.txt"), "w")
    for index, row in df_sample.iterrows():
        bin_id = row[0]
        max_sample_id = row[6]
        if max_sample_id in df_assembly_metadata.index:
            bin_path = os.path.join(output_dir,'bin_fasta/gt10',str(bin_id)+'.fasta')
            f_bin.write(bin_path + '\n')
            f_asm.write(df_assembly_metadata.loc[max_sample_id].at['Sample_Path'] + '\n')
    f_bin.close()
    f_asm.close()

def BinSubjectCount(subjectTableAbundance, cohortMetadataFile, output_dir):
    df_subject_abundance = pd.read_csv(subjectTableAbundance, delimiter="\t",keep_default_na=False)
    df_cohort = pd.read_csv(cohortMetadataFile,keep_default_na=False)
    df_cohort = df_cohort.drop(columns=['Visits'])
    #Add metadata
    df_subject_abundance = pd.merge(df_cohort,df_subject_abundance,on=['Sample'], how='inner')

    #Print bodysite counts
    df_subject_bodysite = df_cohort.drop(columns=['Sample'])
    df_subject_bodysite.drop_duplicates(keep='first', inplace=True)
    df_subject_count_bodysite = df_subject_bodysite.groupby(['Cohort','Bodysite','BodyAggSite','Subject_Status'])['Subject'].count().reset_index(name='Subject_Count')

    df_subject_bin_counts = df_subject_abundance.drop(columns=['Sample'])
    df_subject_bin_counts.rename(columns={'Subject': 'Subject_Count'}, inplace=True)
    df_subject_bin_counts = df_subject_bin_counts[0:0]
    bin_cols = df_subject_bin_counts.columns
    bin_cols = bin_cols[5:]
    for row_index, row in df_subject_count_bodysite.iterrows():
        df_subject_bin_counts = df_subject_bin_counts.append({'Cohort': row[0],
                                                              'Bodysite': row[1],
                                                              'BodyAggSite': row[2],
                                                              'Subject_Status':row[3],
                                                              'Subject_Count':row[4]}, ignore_index=True)
        for bin in bin_cols:
            bin_ctr = df_subject_abundance.loc[(df_subject_abundance['Cohort'] == row[0]) &
                                                 (df_subject_abundance['Bodysite'] == row[1]) &
                                                 (df_subject_abundance['BodyAggSite'] == row[2]) &
                                                 (df_subject_abundance['Subject_Status'] == row[3]) &
                                                 (df_subject_abundance[bin] > 0)]
            bin_ctr.drop_duplicates(subset=['Subject','Cohort','Bodysite','BodyAggSite','Subject_Status'], keep='first',
                                    inplace=True)
            df_subject_bin_counts.loc[(df_subject_bin_counts['Cohort'] == row[0]) &
                                                 (df_subject_bin_counts['Bodysite'] == row[1]) &
                                                 (df_subject_bin_counts['BodyAggSite'] == row[2]) &
                                                 (df_subject_bin_counts['Subject_Status'] == row[3]),bin] = len(bin_ctr)

    cols = list(df_subject_bin_counts.columns)
    cols.remove('Subject_Count')
    cols.insert(4,'Subject_Count')
    df_subject_bin_counts = df_subject_bin_counts[cols]
    out_subject_bin_count = os.path.join(output_dir, 'BinSubjectCounts_Bodysite.tsv')
    df_subject_bin_counts.to_csv(out_subject_bin_count, index=False, sep='\t')


    # Print bodysite aggregate counts
    df_subject_bodysite_agg = df_cohort.drop(columns=['Sample','Bodysite'])
    df_subject_bodysite_agg.drop_duplicates(keep='first', inplace=True)
    df_subject_count_bodysite_agg = df_subject_bodysite_agg.groupby(['Cohort','BodyAggSite','Subject_Status'])['Subject'].count().reset_index(name='Subject_Count')
    #out_subject_bodysite_agg_count = "C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/combined_2/SubjectTotals_Bodysite_Agg.tsv"
    #df_subject_count_bodysite_agg.to_csv(out_subject_bodysite_agg_count, index=False, sep='\t')

    df_agg_subject_bin_counts = df_subject_abundance.drop(columns=['Sample','Bodysite'])
    df_agg_subject_bin_counts.rename(columns={'Subject': 'Subject_Count'}, inplace=True)
    df_agg_subject_bin_counts = df_agg_subject_bin_counts[0:0]
    bin_cols = df_agg_subject_bin_counts.columns
    bin_cols = bin_cols[4:]
    for row_index, row in df_subject_count_bodysite_agg.iterrows():
        df_agg_subject_bin_counts = df_agg_subject_bin_counts.append({'Cohort': row[0],
                                                              'BodyAggSite': row[1],
                                                              'Subject_Status':row[2],
                                                              'Subject_Count':row[3]}, ignore_index=True)
        for bin in bin_cols:
            bin_ctr = df_subject_abundance.loc[(df_subject_abundance['Cohort'] == row[0]) &
                                                 (df_subject_abundance['BodyAggSite'] == row[1]) &
                                                 (df_subject_abundance['Subject_Status'] == row[2]) &
                                                 (df_subject_abundance[bin] > 0)]
            bin_ctr.drop_duplicates(subset=['Subject','Cohort','BodyAggSite','Subject_Status'],keep='first', inplace=True)
            df_agg_subject_bin_counts.loc[(df_agg_subject_bin_counts['Cohort'] == row[0]) &
                                                 (df_agg_subject_bin_counts['BodyAggSite'] == row[1]) &
                                                 (df_agg_subject_bin_counts['Subject_Status'] == row[2]),bin] = len(bin_ctr)

    cols = list(df_agg_subject_bin_counts.columns)
    cols.remove('Subject_Count')
    cols.insert(3,'Subject_Count')
    df_agg_subject_bin_counts = df_agg_subject_bin_counts[cols]
    out_agg_subject_bin_count = os.path.join(output_dir, 'BinSubjectCounts_BodysiteAgg.tsv')
    df_agg_subject_bin_counts.to_csv(out_agg_subject_bin_count, index=False, sep='\t')

def BinSampleCount(subjectTableAbundance, cohortMetadataFile, output_dir):

    df_subject_abundance = pd.read_csv(subjectTableAbundance, delimiter="\t",keep_default_na=False)
    df_cohort = pd.read_csv(cohortMetadataFile,keep_default_na=False)
    df_cohort = df_cohort.drop(columns=['Visits'])
    #Add metadata
    df_subject_abundance = pd.merge(df_cohort,df_subject_abundance,on=['Sample'], how='inner')

    #Print bodysite counts
    df_subject_bodysite = df_cohort.drop(columns=['Subject'])
    df_subject_count_bodysite = df_subject_bodysite.groupby(['Cohort','Bodysite','BodyAggSite','Subject_Status'])['Sample'].count().reset_index(name='Sample_Count')

    df_subject_bin_counts = df_subject_abundance.drop(columns=['Subject'])
    df_subject_bin_counts.rename(columns={'Sample': 'Sample_Count'}, inplace=True)
    df_subject_bin_counts = df_subject_bin_counts[0:0]
    bin_cols = df_subject_bin_counts.columns
    bin_cols = bin_cols[5:]
    for row_index, row in df_subject_count_bodysite.iterrows():
        df_subject_bin_counts = df_subject_bin_counts.append({'Cohort': row[0],
                                                              'Bodysite': row[1],
                                                              'BodyAggSite': row[2],
                                                              'Subject_Status':row[3],
                                                              'Sample_Count':row[4]}, ignore_index=True)
        for bin in bin_cols:
            bin_ctr = df_subject_abundance.loc[(df_subject_abundance['Cohort'] == row[0]) &
                                                 (df_subject_abundance['Bodysite'] == row[1]) &
                                                 (df_subject_abundance['BodyAggSite'] == row[2]) &
                                                 (df_subject_abundance['Subject_Status'] == row[3]) &
                                                 (df_subject_abundance[bin] > 0)]
            df_subject_bin_counts.loc[(df_subject_bin_counts['Cohort'] == row[0]) &
                                                 (df_subject_bin_counts['Bodysite'] == row[1]) &
                                                 (df_subject_bin_counts['BodyAggSite'] == row[2]) &
                                                 (df_subject_bin_counts['Subject_Status'] == row[3]),bin] = len(bin_ctr)

    cols = list(df_subject_bin_counts.columns)
    cols.remove('Sample_Count')
    cols.insert(4,'Sample_Count')
    df_subject_bin_counts = df_subject_bin_counts[cols]

    out_subject_bin_count = os.path.join(output_dir, 'BinSampleCounts_Bodysite.tsv')
    df_subject_bin_counts.to_csv(out_subject_bin_count, index=False, sep='\t')


    # Print bodysite aggregate counts
    df_subject_bodysite_agg = df_cohort.drop(columns=['Subject','Bodysite'])
    df_subject_count_bodysite_agg = df_subject_bodysite_agg.groupby(['Cohort','BodyAggSite','Subject_Status'])['Sample'].count().reset_index(name='Sample_Count')

    df_agg_subject_bin_counts = df_subject_abundance.drop(columns=['Subject','Bodysite'])
    df_agg_subject_bin_counts.rename(columns={'Sample': 'Sample_Count'}, inplace=True)
    df_agg_subject_bin_counts = df_agg_subject_bin_counts[0:0]
    bin_cols = df_agg_subject_bin_counts.columns
    bin_cols = bin_cols[4:]
    for row_index, row in df_subject_count_bodysite_agg.iterrows():
        df_agg_subject_bin_counts = df_agg_subject_bin_counts.append({'Cohort': row[0],
                                                              'BodyAggSite': row[1],
                                                              'Subject_Status':row[2],
                                                              'Sample_Count':row[3]}, ignore_index=True)
        for bin in bin_cols:
            bin_ctr = df_subject_abundance.loc[(df_subject_abundance['Cohort'] == row[0]) &
                                                 (df_subject_abundance['BodyAggSite'] == row[1]) &
                                                 (df_subject_abundance['Subject_Status'] == row[2]) &
                                                 (df_subject_abundance[bin] > 0)]
            df_agg_subject_bin_counts.loc[(df_agg_subject_bin_counts['Cohort'] == row[0]) &
                                                 (df_agg_subject_bin_counts['BodyAggSite'] == row[1]) &
                                                 (df_agg_subject_bin_counts['Subject_Status'] == row[2]),bin] = len(bin_ctr)

    cols = list(df_agg_subject_bin_counts.columns)
    cols.remove('Sample_Count')
    cols.insert(3,'Sample_Count')
    df_agg_subject_bin_counts = df_agg_subject_bin_counts[cols]
    out_agg_subject_bin_count = os.path.join(output_dir, 'BinSampleCounts_BodysiteAgg.tsv')
    df_agg_subject_bin_counts.to_csv(out_agg_subject_bin_count, index=False, sep='\t')

def HMP_BinConsistancy(sampleTableAbundance,cohortMetadataFile,output_dir):
    df_sample_abundance = pd.read_csv(sampleTableAbundance, delimiter="\t")
    df_bodysite = pd.read_csv(cohortMetadataFile, delimiter=",",keep_default_na=False)

    df_subject_visit_abundance = pd.merge(df_bodysite, df_sample_abundance, on=['Sample'], how='inner')

    df_subject_visit_abundance = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Cohort'] == 'HMP') |
                                                                (df_subject_visit_abundance['Cohort'] == 'HMP-1_2') |
                                                                (df_subject_visit_abundance['Cohort'] == 'MetaHit')]
    df_subject_visit_abundance["Visits"] = pd.to_numeric(df_subject_visit_abundance["Visits"])
    df_hmp_multivisit = df_subject_visit_abundance.groupby(['Subject','Cohort','Bodysite','BodyAggSite','Subject_Status'])['Subject'].count().reset_index(name='Subject_Count')
    df_bodysite_sample = pd.DataFrame(columns=['Subject', 'Cohort', 'Bodysite', 'BodyAggSite', 'Bin', 'Visit 1 Sample', 'Visit 1 Bin Abundance',
                 'Visit 2 Sample', 'Visit 2 Bin Abundance',
                 'Visit 3 Sample', 'Visit 3 Bin Abundance'])
    bin_cols = df_subject_visit_abundance.columns
    bin_cols = bin_cols[7:]
    for row_index, row in df_hmp_multivisit.iterrows():
        for bin in bin_cols:
            min_visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
                                                       (df_subject_visit_abundance['Cohort'] == row[1]) &
                                                       (df_subject_visit_abundance['Bodysite'] == row[2]) &
                                                       (df_subject_visit_abundance['BodyAggSite'] == row[3]) &
                                                       (df_subject_visit_abundance['Subject_Status'] == row[4]) &
                                                        (df_subject_visit_abundance[bin] > 0)]
            if len(min_visit_ctr) == 0:
                continue
            visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
                                                       (df_subject_visit_abundance['Cohort'] == row[1]) &
                                                       (df_subject_visit_abundance['Bodysite'] == row[2]) &
                                                       (df_subject_visit_abundance['BodyAggSite'] == row[3]) &
                                                       (df_subject_visit_abundance['Subject_Status'] == row[4])]

            if len(visit_ctr.loc[visit_ctr['Visits'] == 1]) == 1 and len(visit_ctr.loc[visit_ctr['Visits'] == 2]) == 1 and len(visit_ctr.loc[visit_ctr['Visits'] == 3]) == 1:

                # Update visit 1
                visit_1 = visit_ctr.loc[visit_ctr['Visits'] == 1]
                df_bodysite_sample = df_bodysite_sample.append({'Subject': row[0],
                                                                'Cohort': row[1],
                                                                'Bodysite': row[2],
                                                                'BodyAggSite': row[3],
                                                                'Bin': bin,
                                                                'Visit 1 Sample': visit_1['Sample'].item(),
                                                                'Visit 1 Bin Abundance': visit_1[bin].item()}, ignore_index=True)
                # Update visit 2
                visit_2 = visit_ctr.loc[visit_ctr['Visits'] == 2]
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 2 Sample'] = visit_2['Sample'].item()
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 2 Bin Abundance'] = visit_2[bin].item()

                # Update visit 3
                visit_3 = visit_ctr.loc[visit_ctr['Visits'] == 3]
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 3 Sample'] = visit_3['Sample'].item()
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 3 Bin Abundance'] = visit_3[bin].item()

            elif len(visit_ctr.loc[visit_ctr['Visits'] == 1]) == 1 and len(visit_ctr.loc[visit_ctr['Visits'] == 2]) == 1:
                # Update visit 1
                visit_1 = visit_ctr.loc[visit_ctr['Visits'] == 1]
                df_bodysite_sample = df_bodysite_sample.append({'Subject': row[0],
                                                                'Cohort': row[1],
                                                                'Bodysite': row[2],
                                                                'BodyAggSite': row[3],
                                                                'Bin': bin,
                                                                'Visit 1 Sample': visit_1['Sample'].item(),
                                                                'Visit 1 Bin Abundance': visit_1[bin].item()}, ignore_index=True)
                # Update visit 2
                visit_2 = visit_ctr.loc[visit_ctr['Visits'] == 2]
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 2 Sample'] = visit_2['Sample'].item()
                df_bodysite_sample.loc[(df_bodysite_sample['Subject'] == row[0]) &
                                               (df_bodysite_sample['Cohort'] == row[1]) &
                                               (df_bodysite_sample['Bodysite'] == row[2]) &
                                               (df_bodysite_sample['BodyAggSite'] == row[3]) &
                                               (df_bodysite_sample['Bin'] == bin),'Visit 2 Bin Abundance'] = visit_2[bin].item()
    df_bodysite_sample = df_bodysite_sample.sort_values(by=['Visit 3 Bin Abundance', 'Visit 2 Bin Abundance', 'Visit 1 Bin Abundance'],
                                                      ascending=[False, False, False])

    out_bin_visits = os.path.join(output_dir, "BodysiteBinVisitsAbundances.tsv")
    df_bodysite_sample.to_csv(out_bin_visits, index=False, sep='\t')

def HMP_BinSampleAbund(sampleTableAbundance, cohortMetadataFile, output_dir):
    df_sample_abundance = pd.read_csv(sampleTableAbundance, delimiter="\t")
    df_bodysite = pd.read_csv(cohortMetadataFile, delimiter=",",keep_default_na=False)

    df_subject_visit_abundance = pd.merge(df_bodysite, df_sample_abundance, on=['Sample'], how='inner')
    df_subject_visit_abundance = df_subject_visit_abundance.drop(columns=['Bodysite', 'Subject_Status'])
    df_subject_visit_abundance = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Cohort'] == 'HMP') |
                                                                (df_subject_visit_abundance['Cohort'] == 'HMP-1_2') |
                                                                (df_subject_visit_abundance['Cohort'] == 'MetaHit')]

    df_subject_visit_abundance["Visits"] = pd.to_numeric(df_subject_visit_abundance["Visits"])
    df_hmp_multivisit = df_subject_visit_abundance.groupby(['Subject','Cohort','BodyAggSite'])['Subject'].count().reset_index(name='Subject_Count')
    df_bodysite_sample = pd.DataFrame(columns=['Subject', 'Cohort', 'BodyAggSite', 'Bin', 'Total Visit 1 Samples', 'Visit 1 Samples',
                                               'Total Visit 2 Samples', 'Visit 2 Samples',
                                               'Total Visit 3 Samples', 'Visit 3 Samples'])
    bin_cols = df_subject_visit_abundance.columns
    bin_cols = bin_cols[5:]
    for row_index, row in df_hmp_multivisit.iterrows():
        visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
                                                   (df_subject_visit_abundance['Cohort'] == row[1]) &
                                                   (df_subject_visit_abundance['BodyAggSite'] == row[2])]
        all_visit_1 = visit_ctr.loc[visit_ctr['Visits'] == 1]
        all_visit_2 = visit_ctr.loc[visit_ctr['Visits'] == 2]
        all_visit_3 = visit_ctr.loc[visit_ctr['Visits'] == 3]

        for bin in bin_cols:
            min_visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
                                                       (df_subject_visit_abundance['Cohort'] == row[1]) &
                                                       (df_subject_visit_abundance['BodyAggSite'] == row[2]) &
                                                        (df_subject_visit_abundance[bin] > 0)]
            if len(min_visit_ctr) == 0:
                continue


            visit_1 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 1]
            visit_2 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 2]
            visit_3 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 3]
            df_bodysite_sample = df_bodysite_sample.append({'Subject': row[0],
                                                            'Cohort': row[1],
                                                            'BodyAggSite': row[2],
                                                            'Bin': bin,
                                                            'Total Visit 1 Samples': len(all_visit_1),
                                                            'Visit 1 Samples': len(visit_1),
                                                            'Total Visit 2 Samples': len(all_visit_2),
                                                            'Visit 2 Samples': len(visit_2),
                                                            'Total Visit 3 Samples': len(all_visit_3),
                                                            'Visit 3 Samples': len(visit_3)},
                                                           ignore_index=True)
    df_bodysite_sample = df_bodysite_sample.sort_values(by=['Visit 3 Samples', 'Visit 2 Samples', 'Visit 1 Samples'],
                                                      ascending=[False, False, False])
    out_bin_visits = os.path.join(output_dir, "SubjectBinVisitCombined.tsv")
    df_bodysite_sample.to_csv(out_bin_visits, index=False, sep='\t')

def SampleStacked(sampleTableAbundance, cohortMetadataFile, output_dir):
    df_sample_abundance = pd.read_csv(sampleTableAbundance, delimiter="\t")
    df_bodysite = pd.read_csv(cohortMetadataFile, delimiter=",",keep_default_na=False)

    df_subject_abundance = pd.merge(df_bodysite, df_sample_abundance, on=['Sample'], how='inner')
    df_subject_abundance = df_subject_abundance.drop(columns=['Visits'])

    col_list = ['Sample','Subject','Cohort','Bodysite','BodyAggSite','Subject_Status']
    for col_name in reversed(col_list):
        temp_col = df_subject_abundance.pop(col_name)
        df_subject_abundance.insert(0,col_name,temp_col)

    out_sample_abund = os.path.join(output_dir, "SampleAbundanceBodysite.tsv")
    df_subject_abundance.to_csv(out_sample_abund, index=False, sep='\t')

    df_subject_abundance_stacked = pd.melt(df_subject_abundance,id_vars=col_list)
    df_subject_abundance_stacked.rename(columns={'variable': 'Bin',
                                                 'value': 'read_abundance'}, inplace=True)

    out_stacked = os.path.join(output_dir, "SampleAbundanceBodysite_Stacked.tsv")
    df_subject_abundance_stacked.to_csv(out_stacked, index=False, sep='\t')

def mbgcanalytics(metabgc_op_dir,cohort_metadata_file,assembly_metadata_file,output_dir,cpu):
    try:
        CPU_THREADS = 4
        if cpu is not None:
            CPU_THREADS = int(cpu)
        readTableAbundance = os.path.join(metabgc_op_dir,'ReadLevelAbundance.tsv')
        sampleTableAbundance = os.path.join(metabgc_op_dir,'SampleAbundanceMatrix.tsv')
        spHMMFile = os.path.join(metabgc_op_dir, "identify_result", 'spHMM-filtered-results.txt')
        blastResultFile = os.path.join(metabgc_op_dir, "identify_result",'identified-biosynthetic-reads-blast.txt')
        fastaFile = os.path.join(metabgc_op_dir, "identify_result",'identified-biosynthetic-reads.fasta')

        if not os.path.exists(readTableAbundance) and \
            not os.path.exists(sampleTableAbundance) and \
            not os.path.exists(spHMMFile) and \
            not os.path.exists(readTableAbundance):
            print("Metabgc-analytics has failed because metabgc_output_dir is not pointing to a successful run of metabgc search")

        if not os.path.exists(blastResultFile):
            blastn_search_directory = os.path.join(output_dir, 'blastx_result')
            os.makedirs(blastn_search_directory, 0o777, True)
            records = list(SeqIO.parse(fastaFile, "fasta"))
            write_rec = []
            output_file_list=[]
            file_count = 0
            seq_count = math.ceil(len(records)/cpu) + 1
            for seq_rec in records:
                if len(write_rec) < seq_count:
                    write_rec.append(seq_rec)
                else:
                    output_file = os.path.join(blastn_search_directory, str(file_count)+"_seg.fasta")
                    count = SeqIO.write(write_rec, output_file, "fasta")
                    output_file_list.append(output_file)
                    file_count = file_count + 1
                    write_rec.clear()
            if len(write_rec)>0:
                output_file = os.path.join(blastn_search_directory, str(file_count) + "_seg.fasta")
                count = SeqIO.write(write_rec, output_file, "fasta")
                output_file_list.append(output_file)
                file_count = file_count + 1
                write_rec.clear()

            RunPCBlastSearch("nr", output_file_list, "blastx", "-max_target_seqs 1 -outfmt \"6 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle\" ", blastn_search_directory, CPU_THREADS)
            found_hit_ctr=0
            with open(blastResultFile, 'w') as outfile:
                for subdir, dirs, files in os.walk(blastn_search_directory):
                    for file in files:
                        filePath = os.path.join(subdir, file)
                        if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                            with open(filePath) as infile:
                                for line in infile:
                                    outfile.write(line)
                                    found_hit_ctr = found_hit_ctr + 1
            print("Metabgc-analytics BLAST search is complete.")

        #Generate read level analytics table
        ReadLevelBinAnalytics(readTableAbundance,
                              sampleTableAbundance,
                              spHMMFile,
                              blastResultFile,
                              fastaFile,
                              cohort_metadata_file,
                              assembly_metadata_file,
                              output_dir)
        #Generate sample level analytics table
        BinSubjectCount(sampleTableAbundance, cohort_metadata_file, output_dir)

        #Generate bodysite aggregate file
        BinSampleCount(sampleTableAbundance, cohort_metadata_file, output_dir)

        #HMP bin abundance for each subject for each visit grouped by cohort and bodysite
        HMP_BinConsistancy(sampleTableAbundance, cohort_metadata_file, output_dir)

        # HMP abundance for each bin for each subject across visits
        HMP_BinSampleAbund(sampleTableAbundance, cohort_metadata_file, output_dir)

        #Sample abundance pivot table
        SampleStacked(sampleTableAbundance, cohort_metadata_file, output_dir)
    except:
        print("Metabgc-analytics has failed. Please check your paths are correct and contact support on : https://github.com/donia-lab/MetaBGC")
        exit()
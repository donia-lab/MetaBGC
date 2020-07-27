import pytest
from metabgc.src.metabgcanalytics import mbgcanalytics
import os
import pandas as pd

def test1():
    data_dir = "AbcK/output/analytics_test"
    op_dir = "AbcK/output/analytics_test/output"
    cohort_metadata_file = "../metadata/combined_cohort_v2.csv"
    assembly_metadata_file = "../metadata/DoniaSampleAssemblyFileList.csv"
    BLASTDB = "nr"
    mbgcanalytics(data_dir, cohort_metadata_file, assembly_metadata_file, op_dir, 1)

# def test3():
#     df1 = pd.DataFrame({'subject': ['V1.CD1', 'V1.CD10', 'V1.CD12', 'V1.CD15', 'V1.CD15', 'V1.CD17', 'V1.CD17',
#                                     'V1.CD18', 'V1.CD18', 'V1.CD20', 'V1.CD21', 'V1.CD21', 'V1.CD25', 'V1.CD25',
#                                     'V1.CD35', 'V1.CD35', 'V1.CD41', 'V1.CD6', 'V1.CD6', 'V1.CD7', 'V1.CD7', 'O2.UC1',
#                                     'O2.UC11', 'O2.UC11', 'O2.UC1', 'O2.UC12', 'O2.UC12', 'O2.UC13', 'O2.UC13',
#                                     'O2.UC14', 'O2.UC14', 'O2.UC16', 'O2.UC16', 'O2.UC17', 'O2.UC17', 'O2.UC18',
#                                     'O2.UC18', 'O2.UC19', 'O2.UC19', 'O2.UC2', 'O2.UC20', 'O2.UC20', 'O2.UC21',
#                                     'O2.UC21', 'O2.UC2', 'O2.UC22', 'O2.UC22', 'O2.UC23', 'O2.UC23', 'O2.UC24',
#                                     'O2.UC24', 'O2.UC26', 'O2.UC26', 'O2.UC27', 'O2.UC27', 'O2.UC29', 'O2.UC29',
#                                     'O2.UC3', 'O2.UC30', 'O2.UC30', 'O2.UC31', 'O2.UC31', 'O2.UC3', 'O2.UC32',
#                                     'O2.UC32', 'O2.UC34', 'O2.UC34', 'O2.UC35', 'O2.UC35', 'O2.UC36', 'O2.UC37',
#                                     'O2.UC37', 'O2.UC38', 'O2.UC38', 'O2.UC4', 'O2.UC40', 'O2.UC40', 'O2.UC41',
#                                     'O2.UC41', 'O2.UC4', 'O2.UC43', 'O2.UC43', 'O2.UC44', 'O2.UC44', 'O2.UC46',
#                                     'O2.UC46', 'O2.UC47', 'O2.UC47', 'O2.UC48', 'O2.UC48', 'O2.UC5', 'O2.UC5', 'O2.UC9',
#                                     'O2.UC9', 'V1.UC1', 'V1.UC10', 'V1.UC10', 'V1.UC11', 'V1.UC11', 'V1.UC12',
#                                     'V1.UC12', 'V1.UC13', 'V1.UC13', 'V1.UC14', 'V1.UC14', 'V1.UC15', 'V1.UC15',
#                                     'V1.UC17', 'V1.UC17', 'V1.UC2', 'V1.UC21', 'V1.UC21', 'V1.UC23', 'V1.UC23',
#                                     'V1.UC2', 'V1.UC25', 'V1.UC25', 'V1.UC26', 'V1.UC26', 'V1.UC3', 'V1.UC31',
#                                     'V1.UC31', 'V1.UC3', 'V1.UC35', 'V1.UC35', 'V1.UC38', 'V1.UC39', 'V1.UC4',
#                                     'V1.UC40', 'V1.UC4', 'V1.UC45', 'V1.UC47', 'V1.UC49', 'V1.UC49', 'V1.UC5',
#                                     'V1.UC50', 'V1.UC51', 'V1.UC51', 'V1.UC52', 'V1.UC52', 'V1.UC5', 'V1.UC53',
#                                     'V1.UC53', 'V1.UC54', 'V1.UC55', 'V1.UC55', 'V1.UC56', 'V1.UC58', 'O2.UC49',
#                                     'O2.UC49', 'O2.UC50', 'O2.UC50', 'O2.UC51', 'O2.UC51', 'O2.UC52', 'O2.UC52',
#                                     'O2.UC53', 'O2.UC53', 'O2.UC54', 'O2.UC54', 'O2.UC55', 'O2.UC55', 'O2.UC56',
#                                     'O2.UC56', 'O2.UC57', 'O2.UC57', 'O2.UC58', 'O2.UC58', 'O2.UC59', 'O2.UC59',
#                                     'O2.UC60', 'O2.UC60', 'V1.CD11', 'V1.CD13', 'V1.CD14', 'V1.CD16', 'V1.CD19',
#                                     'V1.CD2', 'V1.CD22', 'V1.CD24', 'V1.CD27', 'V1.CD28', 'V1.CD29', 'V1.CD3',
#                                     'V1.CD30', 'V1.CD31', 'V1.CD34', 'V1.CD36', 'V1.CD38', 'V1.CD4', 'V1.CD42',
#                                     'V1.CD8', 'V1.CD9', 'V1.UC16', 'V1.UC18', 'V1.UC19', 'V1.UC24', 'V1.UC27',
#                                     'V1.UC28', 'V1.UC29', 'V1.UC30', 'V1.UC32', 'V1.UC33', 'V1.UC34', 'V1.UC36',
#                                     'V1.UC37', 'V1.UC41', 'V1.UC42', 'V1.UC43', 'V1.UC44', 'V1.UC46', 'V1.UC48',
#                                     'V1.UC57', 'V1.UC6', 'V1.UC61', 'V1.UC62', 'V1.UC7', 'V1.UC8', 'V1.UC9'],
#                         'visits': [0, 0, 0, 0, 3, 0, 4, 0, 3, 4, 0, 4, 0, 4, 0, 1, 0, 0, 4, 0, 4, 0, 0, 2, 2, 0, 2, 0,
#                                    2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 0, 2, 0, 2, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0,
#                                    2, 0, 0, 2, 0, 2, 2, 0, 2, 0, 2, 0, 2, 2, 0, 2, 0, 2, 0, 0, 2, 0, 2, 2, 0, 2, 0, 2,
#                                    0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 0, 2, 0, 5, 0, 4, 0, 3, 0, 1, 0, 3, 0, 2, 0, 0, 4,
#                                    0, 1, 4, 0, 1, 0, 4, 0, 0, 4, 2, 0, 4, 0, 0, 0, 0, 5, 0, 0, 0, 1, 0, 1, 0, 4, 0, 1,
#                                    3, 0, 4, 0, 0, 4, 0, 0, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
#                                    0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]})
#
#     df_out = pd.DataFrame(columns=['Sample', 'Visit'])
#
#     processed_subjects = []
#
#     for index, row in df1.iterrows():
#         if row[0] not in processed_subjects:
#             df_found = df1.loc[df1['subject'] == row[0]]
#             df_found.sort_values(by=['visits'], inplace=True)
#             visit_ctr = 1
#             for index_f, row_f in df_found.iterrows():
#                 sample_id = row_f[0] + "." + str(row_f[1])
#                 df_out.loc[len(df_out)] = [sample_id, visit_ctr]
#                 visit_ctr = visit_ctr + 1
#             processed_subjects.append(row[0])
#     df_out.to_csv('C:/Users/ab50/Documents/git/MetaBGC/MetaBGC-Development/metabgc/metadata/out.csv', index=False, sep=',')

# def test2():
#     data_dir = "C:/Users/ab50/Documents/data/MetaBGCRuns"
#     cohort_metadata_file = "../metadata/combined_cohort.csv"
#     sample_file = os.path.join(data_dir,"Cluster-Classification.csv")
#
#     df_sample_abundance = pd.read_csv(sample_file, delimiter=",")
#     df_bodysite = pd.read_csv(cohort_metadata_file, delimiter=",",keep_default_na=False)
#
#     df_sample_abundance = df_sample_abundance.pivot(index='Sample',columns='ClusterClassification', values='Abund')
#     df_sample_abundance.fillna(0, inplace=True)
#
#     df_subject_visit_abundance = pd.merge(df_bodysite, df_sample_abundance, on=['Sample'], how='left')
#     df_subject_visit_abundance = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Bodysite'] == 'supragingival_plaque') |
#                                                                 (df_subject_visit_abundance['Bodysite'] == 'Supragingival_plaque')]
#
#     df_subject_visit_abundance = df_subject_visit_abundance.drop(columns=['Bodysite', 'Subject_Status'])
#     df_subject_visit_abundance = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Cohort'] == 'HMP') |
#                                                                 (df_subject_visit_abundance['Cohort'] == 'HMP-1_2')]
#     df_subject_visit_abundance.fillna(0, inplace=True)
#     df_subject_visit_abundance["Visits"] = pd.to_numeric(df_subject_visit_abundance["Visits"])
#     df_hmp_multivisit = df_subject_visit_abundance.groupby(['Subject','Cohort','BodyAggSite'])['Subject'].count().reset_index(name='Subject_Count')
#     df_bodysite_sample = pd.DataFrame(columns=['Subject', 'Cohort', 'BodyAggSite', 'Bin', 'Total Visit 1 Samples', 'Visit 1 Samples',
#                                                'Total Visit 2 Samples', 'Visit 2 Samples',
#                                                'Total Visit 3 Samples', 'Visit 3 Samples'])
#     bin_cols = df_subject_visit_abundance.columns
#     bin_cols = bin_cols[5:]
#
#     for row_index, row in df_hmp_multivisit.iterrows():
#         visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
#                                                    (df_subject_visit_abundance['Cohort'] == row[1]) &
#                                                    (df_subject_visit_abundance['BodyAggSite'] == row[2])]
#         all_visit_1 = visit_ctr.loc[visit_ctr['Visits'] == 1]
#         all_visit_2 = visit_ctr.loc[visit_ctr['Visits'] == 2]
#         all_visit_3 = visit_ctr.loc[visit_ctr['Visits'] == 3]
#
#         for bin in bin_cols:
#             min_visit_ctr = df_subject_visit_abundance.loc[(df_subject_visit_abundance['Subject'] == row[0]) &
#                                                        (df_subject_visit_abundance['Cohort'] == row[1]) &
#                                                        (df_subject_visit_abundance['BodyAggSite'] == row[2]) &
#                                                         (df_subject_visit_abundance[bin] > 0)]
#             #if len(min_visit_ctr) == 0:
#             #    continue
#
#
#             visit_1 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 1]
#             visit_2 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 2]
#             visit_3 = min_visit_ctr.loc[min_visit_ctr['Visits'] == 3]
#             df_bodysite_sample = df_bodysite_sample.append({'Subject': row[0],
#                                                             'Cohort': row[1],
#                                                             'BodyAggSite': row[2],
#                                                             'Bin': bin,
#                                                             'Total Visit 1 Samples': len(all_visit_1),
#                                                             'Visit 1 Samples': len(visit_1),
#                                                             'Total Visit 2 Samples': len(all_visit_2),
#                                                             'Visit 2 Samples': len(visit_2),
#                                                             'Total Visit 3 Samples': len(all_visit_3),
#                                                             'Visit 3 Samples': len(visit_3)},
#                                                            ignore_index=True)
#
#     df_bodysite_sample_pos = df_bodysite_sample.loc[(df_bodysite_sample['Visit 1 Samples'] > 0)|
#                                                     (df_bodysite_sample['Visit 2 Samples'] > 0)|
#                                                     (df_bodysite_sample['Visit 3 Samples'] > 0)]
#     df_bodysite_sample_counts = df_bodysite_sample_pos.groupby(['Subject','Cohort','BodyAggSite'])['Subject'].count().reset_index(name='Subject_Count')
#     df_bodysite_sample_pos_dis_counts = df_bodysite_sample_counts.loc[df_bodysite_sample_counts['Subject_Count'] > 1]
#
#     df_bodysite_sample_concord = df_bodysite_sample.loc[~df_bodysite_sample.index.isin(df_bodysite_sample.merge(df_bodysite_sample_pos_dis_counts.assign(a='key'), how='left').dropna().index)]
#     df_bodysite_sample_concord = df_bodysite_sample_concord.sort_values(by=['Visit 3 Samples', 'Visit 2 Samples', 'Visit 1 Samples'],
#                                                       ascending=[False, False, False])
#     out_bin_visits = os.path.join(data_dir, "SubjectBinVisitCombined.tsv")
#     df_bodysite_sample_concord.to_csv(out_bin_visits, index=False, sep='\t')
#
#
#     df_bodysite_sample_discord = df_bodysite_sample.loc[~df_bodysite_sample.index.isin(df_bodysite_sample.merge(df_bodysite_sample_concord.assign(a='key'), how='left').dropna().index)]
#     df_bodysite_sample_discord = df_bodysite_sample_discord.sort_values(by=['Visit 3 Samples', 'Visit 2 Samples', 'Visit 1 Samples'],
#                                                       ascending=[False, False, False])
#     out_bin_visits = os.path.join(data_dir, "SubjectBinVisitCombined_Pos_Discord.tsv")
#     df_bodysite_sample_discord.to_csv(out_bin_visits, index=False, sep='\t')
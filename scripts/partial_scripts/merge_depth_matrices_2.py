import pandas as pd
import os
import re
import sys


def extract_abund(depth_file, drop_var = True):
    print("Processing depth: "+depth_file, flush=True)
    depth_df = pd.DataFrame()
    if os.path.exists(depth_file):
        depth_df = pd.read_csv(depth_file, sep='\t')
    col_name_list = depth_df.columns
    for col_name in col_name_list:
        if "-vs-" in col_name:
            new_col_name = re.sub(".*-vs-", '', col_name)
            depth_df.rename(columns = {col_name:new_col_name}, inplace = True)

    # Delete "-var" and other columns
    col_name_list = depth_df.columns
    drop_cols = []
    if drop_var:
        for col_name in col_name_list:
            if col_name.endswith(".bam-var"):
                drop_cols.append(col_name)
    drop_cols.append("totalAvgDepth")
    depth_df.drop(drop_cols, axis=1, inplace=True)
    return depth_df

if __name__ == '__main__':

    base_dir = sys.argv[2]
    sample_name = sys.argv[1]

    binning_dir = os.path.join(base_dir, sample_name,sample_name + '-binning')
    result = [os.path.join(dp, f)
              for dp, dn, filenames in os.walk(binning_dir)
              for f in filenames if (os.path.basename(f).startswith('depth') and
                                     os.path.splitext(f)[1] == '.txt')]
    if len(result) > 0 :
        merged_depth_df = extract_abund(result[0])
        for depth_file in result[1:]:
            depth_df = extract_abund(depth_file)
            depth_df.drop(['contigLen'], axis=1, inplace=True)
            merged_depth_df = pd.merge(merged_depth_df,depth_df,on=['contigName'])

        # Remove the '.bam' from column names
        col_name_list = merged_depth_df.columns
        for col_name in col_name_list:
            if col_name.endswith(".bam"):
                new_col_name = re.sub(".bam$", '', col_name)
                merged_depth_df.rename(columns = {col_name: new_col_name}, inplace = True)
        merged_depth_df.to_csv(os.path.join(binning_dir, 'merged_depth_abund_only.tsv'), sep='\t', index=False)
        merged_depth_df.drop(['contigLen'], axis=1, inplace=True)
        col_ctr = len(merged_depth_df.columns)
        merged_depth_df['totalAvgDepth'] = merged_depth_df.iloc[:, -(col_ctr - 1):-1].sum(axis=1)
        merged_totals_df = merged_depth_df[['contigName','totalAvgDepth']].copy()
        # Generate bigger report with var columns for the binning tools to use
        merged_depth_df = extract_abund(result[0], False)
        for depth_file in result[1:]:
            depth_df = extract_abund(depth_file, False)
            depth_df.drop(['contigLen'], axis=1, inplace=True)
            merged_depth_df = pd.merge(merged_depth_df,depth_df,on=['contigName'])
        merged_depth_df = pd.merge(merged_depth_df,merged_totals_df,on=['contigName'])
        totalAvgDepth_col = merged_depth_df.pop('totalAvgDepth')
        merged_depth_df.insert(2,'totalAvgDepth',totalAvgDepth_col)
        merged_depth_df.to_csv(os.path.join(binning_dir, 'merged_depth_all.tsv'), sep='\t', index=False)
    else:
        print(sample_name + ' no depth matrices found!')

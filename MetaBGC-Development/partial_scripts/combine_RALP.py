import os
import pandas as pd
import glob
import sys


if __name__ == '__main__':
    binning_tool='vamb'
    ralp_op_dir = 'C:/Users/ab50/Documents/data/binning/' + binning_tool + '_bin_assignment'
    all_summary = {}
    sample_order = [''] * 100
    for ralp_file in os.listdir(ralp_op_dir): # Loop over the reference alignment files
        if ralp_file.endswith('summary_RLAP.csv'):
            sample_name = ralp_file.replace('_summary_RLAP.csv','')
            sample_index = int(sample_name.replace('synthetic_1_algae_S',''))
            sample_order[sample_index] = sample_name
            all_summary[sample_name] = {}
            with open(os.path.join(ralp_op_dir,ralp_file),'r') as f_handle:
                # Skip header
                f_handle.readline()
                for line in f_handle:
                    line_tok = line.split(',')
                    ref_name = line_tok[0]
                    coverage_breath = line_tok[4]
                    all_summary[sample_name][ref_name] = coverage_breath
    op_dir = r'C:\Users\ab50\Documents\data\binning'
    df_summary = pd.DataFrame.from_dict(all_summary)
    df_summary = df_summary.fillna(0)
    df_summary = df_summary.reindex(columns=sample_order)
    df_summary = df_summary.sort_values(by=sample_order,ascending=False)
    df_summary.to_csv(os.path.join(op_dir,'synthetic_1_'+ binning_tool + '_reference_breath.csv'))


import os
from Bio import SeqIO
import sys
import pandas as pd

if __name__ == '__main__':
    ref_file = sys.argv[1]
    base_dir = sys.argv[2]

    ref_seq_list = []
    for record in SeqIO.parse(ref_file, "fasta"):
        ref_seq_list.append(str(record.id))

    df_breath_all = pd.DataFrame(ref_seq_list, columns = ['genes'])
    df_depth_all = pd.DataFrame(ref_seq_list, columns = ['genes'])

    rpkm_files = [os.path.join(dp, f)
              for dp, dn, filenames in os.walk(base_dir)
              for f in filenames if '-rpkm.tsv' in f]

    for file in rpkm_files:
        if os.path.getsize(file) == 0:
            continue
        df_file = pd.read_csv(file, sep='\t', header=None)
        sample_name = df_file.iloc[0, 3]
        # Read breath and make dataframe
        df_breath = df_file.iloc[:, [0, 1]]
        df_breath.columns = ['genes', sample_name]
        df_breath[['genes']] = df_breath[['genes']].astype(str)
        # Read depth and make dataframe
        df_depth = df_file.iloc[:, [0, 2]]
        df_depth.columns = ['genes', sample_name]
        df_depth[['genes']] = df_depth[['genes']].astype(str)
        df_breath_all = pd.merge(df_breath_all,df_breath,on=['genes'],how='left')
        df_depth_all = pd.merge(df_depth_all,df_depth,on=['genes'],how='left')
    df_breath_all = df_breath_all.T
    df_depth_all = df_depth_all.T
    df_breath_all.to_csv('quantified_breath.csv', index=False, header=True)
    df_depth_all.to_csv('quantified_depth.csv', index=False, header=True)
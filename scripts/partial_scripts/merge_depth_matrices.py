import pandas as pd
import os
from Bio import SeqIO
import re
def extract_abund(sample_dir, sample_name, contig_list):
    depth_file = os.path.join(sample_dir, sample_name+"-binning", "depth.txt")
    print("Processing depth: "+depth_file, flush=True)
    depth_df = pd.DataFrame()
    if os.path.exists(depth_file):
        depth_df = pd.read_csv(depth_file, sep='\t')
        depth_df = depth_df[depth_df['contigName'].isin(contig_list)]
        depth_df = depth_df.assign(sample=sample_name)
    col_name_list = depth_df.columns
    for col_name in col_name_list:
        if "-vs-" in col_name:
            new_col_name = re.sub(".*-vs-", '', col_name)
            depth_df.rename(columns = {col_name:new_col_name}, inplace = True)
    return depth_df


if __name__ == '__main__':

    sample_dir_file = "/tigress/DONIA/scripts/Algae_100_Metagenome_Binning/binning_run_23_29/all-sample-dirnames"
    sample_dir_dict={}
    with open(sample_dir_file) as sample_list:
        for line in sample_list:
            sample_dir = os.path.join("/tigress/DONIA/data/donia/", line.strip())
            sample_name = os.path.basename(sample_dir)
            sample_dir_dict[sample_name] = sample_dir

    #Load bin info.
    all_contig_abund = pd.DataFrame()
    drep_bin_dir = "/tigress/DONIA/data/donia/Algae_100_Metagenome_Analysis/binning_run_23_29/vamb_drep_output/dereplicated_genomes_decom"
    for filename in os.listdir(drep_bin_dir):
        if filename.endswith(".fna"):
            print("Processing: "+filename, flush=True)
            bin_filename = os.path.join(drep_bin_dir, filename)
            sample_filename_tok = os.path.splitext(filename)[0]
            sample_name = os.path.splitext(sample_filename_tok)[0]
            bin_id = os.path.splitext(sample_filename_tok)[1][1:]
            contig_list = []
            with open(bin_filename) as handle:
                for record in SeqIO.parse(handle, format="fasta"):
                    contig_list.append(record.id)
            bin_abund_df = extract_abund(sample_dir_dict[sample_name], sample_name, contig_list)
            bin_abund_df = bin_abund_df.assign(bin_id=bin_id)
            all_contig_abund = pd.concat([all_contig_abund, bin_abund_df], axis=0)
    all_contig_abund.to_csv("vamb_contig_abund.csv", index=False)

    # Delete "-var" and other columns
    col_name_list = all_contig_abund.columns
    drop_cols = []
    for col_name in col_name_list:
        if col_name.endswith(".bam-var"):
            drop_cols.append(col_name)
    drop_cols.append("contigLen")
    drop_cols.append("contigName")
    drop_cols.append("totalAvgDepth")
    all_contig_abund.drop(drop_cols, axis=1, inplace=True)

    # Trim out ".bam" from column names
    col_name_list = all_contig_abund.columns
    group_cols=[]
    metric_cols=[]
    for col_name in col_name_list:
        if col_name.endswith(".bam"):
            new_col_name = re.sub(".bam$", '', col_name)
            all_contig_abund.rename(columns = {col_name:new_col_name}, inplace = True)
            all_contig_abund[new_col_name] = pd.to_numeric(all_contig_abund[new_col_name])
            metric_cols.append(new_col_name)
        else:
            group_cols.append(col_name)
    aggs = all_contig_abund.groupby(group_cols)[metric_cols].mean()
    all_contig_abund.drop(metric_cols, axis=1, inplace=True)
    all_contig_abund.drop_duplicates(subset=group_cols, keep='last', inplace=True)
    all_contig_abund = all_contig_abund.merge(right=aggs, right_index=True, left_on=group_cols, how='right')

    # Load taxanomy file
    taxa_file = "/tigress/DONIA/data/donia/Algae_100_Metagenome_Analysis/binning_run_23_29/vamb_gtdbtk_taxa/classify/gtdbtk.bac120.summary.tsv"
    df_taxa = pd.DataFrame(columns=['bin_id', 'sample', 'domain', 'phylum','class','order','family','genus','species'])
    with open(taxa_file) as taxa_f:
        taxa_f.readline()
        for line in taxa_f:
            line_tok = line.split('\t')
            sample_filename_tok = os.path.splitext(line_tok[0])[0]
            sample_name =  os.path.splitext(sample_filename_tok)[0]
            bin_id =  os.path.splitext(sample_filename_tok)[1][1:]
            taxa_string = line_tok[1].split(sep=';')
            df_taxa = df_taxa.append({'bin_id': bin_id,
                                        'sample': sample_name,
                                        'domain': taxa_string[0][3:],
                                        'phylum': taxa_string[1][3:],
                                        'class': taxa_string[2][3:],
                                        'order': taxa_string[3][3:],
                                        'family': taxa_string[4][3:],
                                        'genus': taxa_string[5][3:],
                                        'species': taxa_string[6][3:]}, ignore_index=True)
    all_contig_abund = pd.merge(all_contig_abund, df_taxa, on=['bin_id','sample'])
    all_contig_abund.to_csv("vamb_bin_abund.csv", index=False)










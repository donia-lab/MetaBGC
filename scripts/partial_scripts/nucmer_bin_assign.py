import os
import pandas as pd
from Bio import SeqIO
import glob
import numpy as np
import scipy.optimize as sopt    # RLAP solver
import matplotlib.pyplot as plt  # visualizatiion
import seaborn as sns            # """
import sys

def parse_file_name(nucmer_align_file, binning_tool):
    if binning_tool == 'metabat':
        return nucmer_align_file.split('.')[1].split('_')[0]
    else:
        return nucmer_align_file.split('.')[0]


def read_nucmer_align(nucmer_align_file):
    with open(nucmer_align_file , 'r') as f:
        # Skip 3 lines on top
        f.readline()
        f.readline()
        f.readline()
        nucmer_align_data=[]
        row_count = 0
        for line in f:
            if not line.startswith('='): # Skip header separator line
                line_tok_l1 = line.split('|')
                line_tok_l2=[]
                if row_count == 0: # Header row
                    for tok in line_tok_l1:
                        tok = tok.replace(' ','')
                        tok = tok.replace('[','')
                        tok = tok.replace(']',' ')
                        line_tok_l2.extend(tok.split())
                else:
                    for tok in line_tok_l1:
                        line_tok_l2.extend(tok.split())
                nucmer_align_data.append(line_tok_l2)
                row_count = row_count + 1
        headers = nucmer_align_data.pop(0)
        headers.append('AssembledContig')
        df = pd.DataFrame(nucmer_align_data, columns=headers)
        #df['COVR'] = df['COVR'].astype(float)
        for i in headers[:-2]:
            df[i] = df[i].astype(float)
        return df

def compute_assignment(bin_coverages):
    coverages_array = np.array(bin_coverages)
    max_cost = np.amax(coverages_array)
    harvest_profit = max_cost - coverages_array

    row_ind, col_ind = sopt.linear_sum_assignment(harvest_profit)
    return row_ind, col_ind

if __name__ == '__main__':
    ref_genomes_dir = sys.argv[1]
    sample_name = sys.argv[2]
    base_dir = sys.argv[3]
    binning_tool = sys.argv[4]
    nucmer_align_dir= os.path.join(base_dir,sample_name,sample_name + '-binning', binning_tool + '_bins', 'nucmer_align')

    # Load read counts
    read_length = 100.0
    number_of_reads = 20000000
    contig_read_count_file = os.path.join(base_dir,sample_name,sample_name + '_abundance.csv')
    contig_read_counts = {}
    with open(contig_read_count_file , 'r') as f:
        f.readline()
        for line in f:
            line_tok = line.split(',')
            contig = line_tok[1].replace('.fasta', '')
            num_reads = int(line_tok[2])
            if contig not in contig_read_counts:
                contig_read_counts[contig] = num_reads
            else:
                contig_read_counts[contig] = contig_read_counts[contig] + num_reads

    # Load sequence length metadata
    ref_len_dict = {}
    ref_coverage_dict = {}
    ref_count_ratio = {}
    for ref_file in os.listdir(ref_genomes_dir):
        ref_file_path = os.path.join(ref_genomes_dir,ref_file)
        ref_total_len = 0
        ref_total_read_count = 0
        for record in SeqIO.parse(ref_file_path, "fasta"):
            ref_total_len = ref_total_len + len(record)
            if record.id in contig_read_counts:
                ref_total_read_count = ref_total_read_count + contig_read_counts[record.id]
        ref_count_ratio[os.path.splitext(ref_file)[0]] = ref_total_read_count / number_of_reads
        ref_coverage_dict[os.path.splitext(ref_file)[0]] = ref_total_read_count*read_length / ref_total_len
        ref_len_dict[os.path.splitext(ref_file)[0]] = ref_total_len


    # Loop over all the reference directories

    bin_coverages_dict = []
    bin_identities_dict = []
    references_list = []
    for ref_dir in os.listdir(nucmer_align_dir):
        ref_name = ref_dir
        ref_dir = os.path.join(nucmer_align_dir,ref_dir)
        if os.path.isdir(ref_dir):
            references_list.append(ref_name)
            bin_count = len(glob.glob1(ref_dir,"*.coords"))
            coverages = {}
            identities = {}
            for nucmer_align_file in os.listdir(ref_dir): # Loop over the reference alignment files
                if nucmer_align_file.endswith('.coords'):
                    bin_id = parse_file_name(nucmer_align_file, binning_tool)
                    nucmer_align_file = os.path.join(ref_dir, nucmer_align_file)

                    # Read the alignment file and return the alignment in a df
                    df_nucmer_align = read_nucmer_align(nucmer_align_file)

                    # Get the total coverage % of each contig in the reference
                    df_ref_coverage = df_nucmer_align.groupby(['TAGS'])['COVR'].sum().reset_index(name='TotalRefContigCov')
                    # Get the mean % identity of each contig in the reference
                    df_ref_identity = df_nucmer_align.groupby(['TAGS'])['%IDY'].mean().reset_index(name='AvgRefContigIdy')

                    # Get the size of each contig in the reference
                    df_ref_len = df_nucmer_align.groupby(['TAGS','LENR']).size().to_frame(name = 'COUNT').reset_index()

                    # Compute to weighted sum to get the total coverage of the reference
                    df_ref = pd.merge(df_ref_coverage,df_ref_len,on=['TAGS'], how='inner')
                    df_ref = pd.merge(df_ref,df_ref_identity,on=['TAGS'], how='inner')
                    df_ref['WtCoverage'] = df_ref['TotalRefContigCov'] * df_ref['LENR']
                    df_ref['WtCoverage'] = df_ref['WtCoverage'] /100.0
                    total_coverage = df_ref['WtCoverage'].sum()
                    coverages[bin_id] = round(total_coverage/ref_len_dict[ref_name] * 100.0,3)
                    identities[bin_id] = df_ref['AvgRefContigIdy'].mean()
            bin_coverages_dict.append(coverages)
            bin_identities_dict.append(identities)

    # Sort the bin_ids and create a 2D matrix of the coverages
    bin_coverages = []
    bin_ids_mapping = []
    for coverage_dict in bin_coverages_dict:
        coverage_list = []
        bin_id_list = []
        for bin in sorted (coverage_dict) :
            coverage_list.append(coverage_dict[bin])
            bin_id_list.append(bin)
        bin_coverages.append(coverage_list)
        bin_ids_mapping.append(bin_id_list)
    ref_index, bin_index = compute_assignment(bin_coverages)

    # Print coverage matrix
    f = open(os.path.join(base_dir, binning_tool + '_bin_assignment',sample_name + '_ref_coverage_matrix.csv'),'w')
    f.write('reference/bins,')
    for bin in sorted (bin_coverages_dict[0]):
        f.write(bin + ',')
    f.write('\n')
    for j in range(0,len(references_list)):
        f.write(references_list[j] + ',')
        for bin in sorted (bin_coverages_dict[j]):
            f.write(str(bin_coverages_dict[j][bin]) + ',')
        f.write('\n')
    f.close()

    # Print summary of the assignments
    f = open(os.path.join(base_dir, binning_tool + '_bin_assignment',sample_name + '_summary_RLAP.csv'),'w')
    f.write('reference_genome,reference_genome_size,synthetic_read_coverage_DEPTH,bin_id,bin_coverage_BREATH,bin_mean_identity\n')
    for index, ref_i in enumerate(ref_index):
        ref_name = references_list[ref_i]
        bin_i = bin_ids_mapping[ref_i][bin_index[index]]
        print('{0} genome is best covered by bin: {1} with coverage % of {2}'.format(ref_name,bin_i,bin_coverages_dict[ref_i][bin_i]))
        f.write('{0},{1},{2},{3},{4},{5}\n'.format(ref_name,ref_len_dict[ref_name],ref_coverage_dict[ref_name],bin_i,bin_coverages_dict[ref_i][bin_i],bin_identities_dict[ref_i][bin_i]))
    f.close()

# Script to generate a all-vs-all coverage matrix from mummer alignment of multiple
# references to multiple other references or samples
import os
import pandas as pd
from Bio import SeqIO
import glob
import sys


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

if __name__ == '__main__':
    data_dir = sys.argv[1]
    ref_genomes_dir = os.path.join(data_dir, 'fna')
    ref_mummer_dir = os.path.join(data_dir, 'mummer_align')
    output_dir = os.path.join(data_dir, 'output_dir')

    # Load genus taxa file
    taxa_file = os.path.join(data_dir, 'genus-bacterial-genome-taxonomic-groupings.csv')
    df_taxa = pd.read_csv(taxa_file)
    genus_group_dict = {k: v.drop(['group_id','classification'], axis=1)['sample_name'].to_list()
                        for k, v in df_taxa.groupby('group_id')}
    genus_dict = {k: v.drop(['group_id', 'sample_name'], axis=1)['classification'].to_list()
                  for k, v in df_taxa.groupby('group_id')}

    # Load family taxa file
    taxa_file = os.path.join(data_dir, 'family-bacterial-genome-taxonomic-groupings.csv')
    df_taxa = pd.read_csv(taxa_file)
    family_group_dict = {k: v.drop(['group_id', 'classification'], axis=1)['sample_name'].to_list()
                        for k, v in df_taxa.groupby('group_id')}
    family_dict = {k: v.drop(['group_id', 'sample_name'], axis=1)['classification'].to_list()
                  for k, v in df_taxa.groupby('group_id')}

    # Load sequence length metadata
    ref_len_dict = {}
    ref_coverage_dict = {}
    ref_count_ratio = {}
    for ref_file in os.listdir(ref_genomes_dir):
        ref_file_path = os.path.join(ref_genomes_dir, ref_file)
        ref_total_len = 0
        ref_total_read_count = 0
        for record in SeqIO.parse(ref_file_path, "fasta"):
            ref_total_len = ref_total_len + len(record)
        ref_len_dict[os.path.splitext(ref_file)[0]] = ref_total_len

    # Loop over all the reference directories

    bin_coverages_dict = {}
    bin_identities_dict = {}
    references_list = []
    for ref_dir in os.listdir(ref_mummer_dir):
        ref_name = ref_dir
        ref_dir = os.path.join(ref_mummer_dir,ref_dir)
        if os.path.isdir(ref_dir):
            references_list.append(ref_name)
            bin_count = len(glob.glob1(ref_dir,"*.coords"))
            coverages = {}
            identities = {}
            for nucmer_align_file in os.listdir(ref_dir): # Loop over the reference alignment files
                if nucmer_align_file.endswith('.coords'):
                    align_id = os.path.splitext(os.path.basename(nucmer_align_file))[0]
                    nucmer_align_file = os.path.join(ref_dir, nucmer_align_file)

                    # Read the alignment file and return the alignment in a df
                    df_nucmer_align = read_nucmer_align(nucmer_align_file)

                    # Get the total coverage % of each contig in the reference
                    #df_ref_coverage = df_nucmer_align.groupby(['TAGS','AssembledContig'])['COVR'].max().reset_index()
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
                    coverages[align_id] = round(total_coverage/ref_len_dict[ref_name] * 100.0,3)
                    identities[align_id] = df_ref['AvgRefContigIdy'].mean()
            bin_coverages_dict[ref_name] = coverages
            bin_identities_dict[ref_name] = identities

    # Print coverage matrix
    f = open(os.path.join(output_dir, 'nucmer_ref_coverage_matrix.csv'), 'w')
    f.write('reference,')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
    f.write('\n')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
        for idx_2, ref_name_2 in enumerate(references_list):
            f.write(str(bin_coverages_dict[ref_name][ref_name_2]) + ',')
        f.write('\n')
    f.close()

    # Print identity matrix
    f = open(os.path.join(output_dir, 'nucmer_ref_identity_matrix.csv'),'w')
    f.write('reference,')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
    f.write('\n')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
        for idx_2, ref_name_2 in enumerate(references_list):
            f.write(str(bin_identities_dict[ref_name][ref_name_2]) + ',')
        f.write('\n')
    f.close()

    # Group by genus

    genus_coverages_dict = {}
    for k_1, v_1 in genus_group_dict.items():
        genus_coverages_dict[k_1] = {}
        for k_2, v_2 in genus_group_dict.items():
            cov_sum = 0
            for ref_1 in v_1:
                for ref_2 in v_2:
                    cov_sum = cov_sum + bin_coverages_dict[ref_1][ref_2]
            genus_coverages_dict[k_1][k_2] = cov_sum / (len(v_1) * len(v_2))

    f = open(os.path.join(output_dir, 'nucmer_ref_coverage_matrix_genus.csv'), 'w')
    f.write('group_id:genus,')
    for k, v in genus_group_dict.items():
        f.write(str(k) + ':' + genus_dict[k][0] + ',')
    f.write('\n')
    for k_1, v_1 in genus_group_dict.items():
        f.write(str(k_1) + ':' + genus_dict[k_1][0] + ',')
        for k_2, v_2 in genus_group_dict.items():
            f.write(str(genus_coverages_dict[k_1][k_2]) + ',')
        f.write('\n')
    f.close()

    # Group by family

    family_coverages_dict = {}
    for k_1, v_1 in family_group_dict.items():
        family_coverages_dict[k_1] = {}
        for k_2, v_2 in family_group_dict.items():
            cov_sum = 0
            for ref_1 in v_1:
                for ref_2 in v_2:
                    cov_sum = cov_sum + bin_coverages_dict[ref_1][ref_2]
            family_coverages_dict[k_1][k_2] = cov_sum / (len(v_1) * len(v_2))

    f = open(os.path.join(output_dir, 'nucmer_ref_coverage_matrix_family.csv'), 'w')
    f.write('group_id:family,')
    for k, v in family_group_dict.items():
        f.write(str(k) + ':' + family_dict[k][0] + ',')
    f.write('\n')
    for k_1, v_1 in family_group_dict.items():
        f.write(str(k_1) + ':' + family_dict[k_1][0] + ',')
        for k_2, v_2 in family_group_dict.items():
            f.write(str(family_coverages_dict[k_1][k_2]) + ',')
        f.write('\n')
    f.close()
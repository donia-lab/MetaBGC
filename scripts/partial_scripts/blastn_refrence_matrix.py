
# Script to generate a all-vs-all coverage matrix from mummer alignment of multiple
# references to multiple other references or samples
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import sys

if __name__ == '__main__':
    data_dir = sys.argv[1] # r'C:\Users\ab50\Documents\data\binning\algae_genomes\'
    ref_genomes_dir = os.path.join(data_dir, 'fna')
    ref_blastn_dir = os.path.join(data_dir, 'blastn_align')
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
    references_list = []
    for ref_dir in os.listdir(ref_blastn_dir):
        ref_name = ref_dir
        ref_dir = os.path.join(ref_blastn_dir,ref_dir)
        if os.path.isdir(ref_dir):
            references_list.append(ref_name)
            print('Processing:', ref_name)
            coverages = {}
            identities = {}
            for blast_align_file in os.listdir(ref_dir): # Loop over the reference alignment files
                if blast_align_file.endswith('.out'):
                    align_id = blast_align_file.split('_+_')[0]
                    blast_align_file = os.path.join(ref_dir, blast_align_file)

                    # Get the len bitmap of the scaffolds in the reference to number of bases covered
                    contig_coverage_dict = {}
                    ref_file_path = os.path.join(ref_genomes_dir, ref_name+'.fna')
                    for record in SeqIO.parse(ref_file_path, "fasta",):
                        coverage_arr = np.zeros(len(record))
                        contig_coverage_dict[record.id] = coverage_arr

                    # Read the alignment records and compute the bases covered by the alignment
                    # Check and read the blast alignment file
                    # -outfmt "6 sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue"
                    aligned_seq_records = []
                    if os.path.getsize(blast_align_file) > 0:
                        with open(blast_align_file) as f:
                            for line in f:
                                subject_scaf_name = line.split(sep='\t')[0]
                                start_index = int(line.split(sep='\t')[2]) - 1
                                end_index = int(line.split(sep='\t')[3])
                                ident_pct = float(line.split(sep='\t')[9])
                                if ident_pct > 70:
                                    contig_coverage_dict[subject_scaf_name][start_index:end_index] = 1
                        ref_total_len = ref_len_dict[ref_name]
                        total_coverage = sum(np.sum(v) for k, v in contig_coverage_dict.items())
                        coverages[align_id] = round(total_coverage/ref_total_len*100.0, 2)
            bin_coverages_dict[ref_name] = coverages

    # Print coverage matrix
    f = open(os.path.join(output_dir, 'blastn_ref_coverage_matrix.csv'),'w')
    f.write('reference,')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
    f.write('\n')
    for idx, ref_name in enumerate(references_list):
        f.write(ref_name + ',')
        for idx_2, ref_name_2 in enumerate(references_list):
            if ref_name_2 in bin_coverages_dict[ref_name]:
                f.write(str(bin_coverages_dict[ref_name][ref_name_2]) + ',')
            else:
                f.write('0,')
        f.write('\n')
    f.close()

    # Group by genus

    genus_coverages_dict = {}
    for k_1, v_1 in genus_group_dict.items():
        genus_coverages_dict[k_1] = {}
        for k_2, v_2 in genus_group_dict.items():
            cov_sum = 0
            for ref_1 in v_1:
                if ref_1 in bin_coverages_dict:
                    for ref_2 in v_2:
                        if ref_2 in bin_coverages_dict[ref_1]:
                            cov_sum = cov_sum + bin_coverages_dict[ref_1][ref_2]
            genus_coverages_dict[k_1][k_2] = cov_sum / (len(v_1) * len(v_2))

    f = open(os.path.join(output_dir, 'blastn_ref_coverage_matrix_genus.csv'), 'w')
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
                if ref_1 in bin_coverages_dict:
                    for ref_2 in v_2:
                        if ref_2 in bin_coverages_dict[ref_1]:
                            cov_sum = cov_sum + bin_coverages_dict[ref_1][ref_2]
            family_coverages_dict[k_1][k_2] = cov_sum / (len(v_1) * len(v_2))

    f = open(os.path.join(output_dir, 'blastn_ref_coverage_matrix_family.csv'), 'w')
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
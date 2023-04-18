
# Script to generate a all-vs-all coverage matrix from mummer alignment of multiple
# references to multiple other references or samples
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import sys

if __name__ == '__main__':
    ref_genomes_dir = sys.argv[1] # r'C:\Users\ab50\Documents\data\binning\algae_genomes\fna'
    ref_blastn_dir = sys.argv[2] #r'C:\Users\ab50\Documents\data\binning\algae_genomes\blastn_align'

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
    f = open(os.path.join('ref_coverage_matrix.csv'),'w')
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
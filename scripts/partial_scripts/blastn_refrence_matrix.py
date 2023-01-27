import os
import pandas as pd
from Bio import SeqIO
import glob
import sys
import shutil
from statistics import mean 
import numpy as np

if __name__ == '__main__':
    ref_genome = sys.argv[1]
    sample_list_file = sys.argv[2]
    output_dir = sys.argv[3]
    ext_str = '.fna'
    # Load sequence length metadata
    ref_total_len = 0
    ref_total_read_count = 0
    for record in SeqIO.parse(ref_genome, "fasta",):
        ref_total_len = ref_total_len + len(record)

    # Loop over all the sample alignment directories
    bin_coverages_dict = {}
    bin_identities_dict = {}
    bin_aligned_dict = {}
    with open(sample_list_file, 'r') as s_file:
        for line in s_file:
            sample_dir = line.strip().split(sep=',')[0]
            sample_name = line.strip().split(sep=',')[1]
            print('Processing: ' + sample_name)
            bin_dir = os.path.join(sample_dir, sample_name + '-binning/vamb_bins/bins')
            align_dir = os.path.join(bin_dir, '../blastn_align/2716884990/')
            if os.path.isdir(align_dir):
                bin_count = len(glob.glob1(align_dir,"*.out"))
                coverage_pct = {}
                aligned_pct = {}
                identities = {}
                for align_file in os.listdir(align_dir): # Loop over the alignment files
                    if align_file.endswith('.out'):
                        #bin_id = os.path.basename(align_file).split(sep='.fna_')[0] + '.fna'
                        bin_id = os.path.basename(align_file).split(sep='_2716884990')[0]
                        align_file = os.path.join(align_dir, align_file)
                        sequence_file = os.path.join(bin_dir, bin_id + ext_str)

                        bin_total_len = 0
                        bin_seqs_records = {}
                        for record in SeqIO.parse(sequence_file, "fasta",):
                            bin_total_len = bin_total_len + len(record)
                            bin_seqs_records[record.id] = record

                        # Check and read the blast alignment file
                        # -outfmt "6 sseqid slen sstart send qseqid qlen qstart qend pident evalue"
                        aligned_seq_records = []
                        if os.path.getsize(align_file) > 0:
                            coverage_arr = np.zeros(ref_total_len)
                            ident_arr = []
                            bin_aligned_len = 0
                            with open(align_file) as f:
                                for line in f:
                                    start_index = int(line.split(sep='\t')[2])
                                    end_index = int(line.split(sep='\t')[3])
                                    ident_pct = float(line.split(sep='\t')[8])
                                    coverage_arr[start_index:end_index] = 1
                                    ident_arr.append(ident_pct)

                                    q_seq_id = line.split(sep='\t')[4]
                                    q_start_index = int(line.split(sep='\t')[6])
                                    q_end_index = int(line.split(sep='\t')[7])
                                    q_len = int(line.split(sep='\t')[5])
                                    bin_aligned_len = bin_aligned_len + (q_end_index - q_start_index + 1)
                                    if bin_aligned_len > (q_len * 0.5) and q_seq_id in bin_seqs_records:
                                        aligned_seq_records.append(bin_seqs_records[q_seq_id])
                            total_coverage = np.sum(coverage_arr)
                            coverage_pct[bin_id] = round(total_coverage/ref_total_len*100.0, 6)
                            aligned_pct[bin_id] = round(bin_aligned_len/bin_total_len*100.0, 6)
                            identities[bin_id] = round(mean(ident_arr), 6)
                            if coverage_pct[bin_id] > 5.0:
                                bin_file = os.path.join(bin_dir, bin_id + ext_str)
                                shutil.copy2(bin_file, output_dir)
                                extract_bin_file = os.path.join(output_dir, 'cEK_Aligned_' + bin_id + ext_str)
                                SeqIO.write(aligned_seq_records, extract_bin_file, "fasta")
                bin_coverages_dict[sample_name] = coverage_pct
                bin_identities_dict[sample_name] = identities
                bin_aligned_dict[sample_name] = aligned_pct

    # Print coverage matrix
    f = open(os.path.join('ref_coverage_matrix.csv'),'w')
    f.write('sample_name,bin_name,ref_coverage_pct,alignment_identity,bin_aligned_pct')
    f.write('\n')
    for sample_name, sample_coverages in bin_coverages_dict.items():
        for bin, cov in sample_coverages.items():
            if cov > 1.0:
                f.write(sample_name + ',' + bin + ',' + str(cov) + ',' + str(bin_identities_dict[sample_name][bin]) +
                        ',' + str(bin_aligned_dict[sample_name][bin]))
                f.write('\n')
    f.close()
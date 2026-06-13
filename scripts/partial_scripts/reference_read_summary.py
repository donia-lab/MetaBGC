import os
import pandas as pd
from Bio import SeqIO
import glob
import sys

if __name__ == '__main__':
    ref_genomes_dir = sys.argv[1]
    base_dir = sys.argv[2]

    # Load read counts
    read_length = 100.0

    # Open file and write header
    f_ratio = open(os.path.join('synthetic_read_ratio.csv'), 'w')
    f_count = open(os.path.join('synthetic_read_count.csv'), 'w')
    f_cov = open(os.path.join('synthetic_read_coverage.csv'),'w')
    f_ratio.write('sample_name,')
    f_count.write('sample_name,')
    f_cov.write('sample_name,')
    for ref_file in os.listdir(ref_genomes_dir):
        f_ratio.write(os.path.splitext(ref_file)[0] + ',')
        f_count.write(os.path.splitext(ref_file)[0] + ',')
        f_cov.write(os.path.splitext(ref_file)[0]+',')
    f_ratio.write('\n')
    f_count.write('\n')
    f_cov.write('\n')
    for i in range(0,100):
        sample_name = 'synthetic_1_algae_S' + str(i)
        contig_read_count_file = os.path.join(base_dir,sample_name,sample_name + '_abundance.csv')
        contig_read_counts = {}
        total_sample_reads = 0
        with open(contig_read_count_file , 'r') as f:
            for line in f:
                line_tok = line.split(',')
                contig = line_tok[1].replace('.fasta', '')
                num_reads = int(line_tok[2])
                total_sample_reads = total_sample_reads + num_reads
                if contig not in contig_read_counts:
                    contig_read_counts[contig] = num_reads
                else:
                    contig_read_counts[contig] = contig_read_counts[contig] + num_reads
        print('Processing:' + sample_name + ', Size:' + str(sum(contig_read_counts.values())))
        # Load sequence length metadata
        f_ratio.write(sample_name + ',')
        f_count.write(sample_name + ',')
        f_cov.write(sample_name+',')
        dup_list = []
        for ref_file in os.listdir(ref_genomes_dir):
            ref_file_path = os.path.join(ref_genomes_dir,ref_file)
            ref_total_len = 0
            ref_total_read_count = 0
            for record in SeqIO.parse(ref_file_path, "fasta"):
                ref_total_len = ref_total_len + len(record)
                if record.id in contig_read_counts and record.id not in dup_list:
                    ref_total_read_count = ref_total_read_count + contig_read_counts[record.id]
                    dup_list.append(record.id)
                elif record.id in dup_list:
                    print('Dup Found:' + record.id)
            f_ratio.write(str(ref_total_read_count / total_sample_reads) + ',')
            f_count.write(str(ref_total_read_count) + ',')
            f_cov.write(str(ref_total_read_count*read_length / ref_total_len) + ',')
        f_ratio.write('\n')
        f_count.write('\n')
        f_cov.write('\n')
    f_ratio.close()
    f_count.close()
    f_cov.close()

    # Print reference lengths
    with open('ref_lengths.csv','w') as f_ref:
        f_ref.write('Reference,Reference Length\n')
        for ref_file in os.listdir(ref_genomes_dir):
            ref_file_path = os.path.join(ref_genomes_dir,ref_file)
            ref_total_len = 0
            ref_total_read_count = 0
            for record in SeqIO.parse(ref_file_path, "fasta"):
                ref_total_len = ref_total_len + len(record)
            f_ref.write(os.path.splitext(ref_file)[0]+','+str(ref_total_len)+'\n')


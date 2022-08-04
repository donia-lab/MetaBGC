import os
from Bio import SeqIO
import sys


def parse_file_name(rpkm_file):
    return rpkm_file.replace('-rpkm.tsv','')

def gather_reference_breath(ref_genomes_dir, sample_name, base_dir):
    print('Processing: ' + sample_name,flush=True)
    quant_dir= os.path.join(base_dir,sample_name,sample_name + '-ref_align', sample_name + '-ref_align')
    # Load read counts
    read_length = 100.0
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
    for ref_file in os.listdir(ref_genomes_dir):
        ref_name = os.path.splitext(ref_file)[0]
        ref_file_path = os.path.join(ref_genomes_dir,ref_file)
        ref_total_len = 0
        ref_total_read_count = 0
        ref_len_dict[ref_name] = {}
        for record in SeqIO.parse(ref_file_path, "fasta"):
            ref_total_len = ref_total_len + len(record)
            ref_len_dict[ref_name][record.id] = len(record)
            if record.id in contig_read_counts:
                ref_total_read_count = ref_total_read_count + contig_read_counts[record.id]
        ref_coverage_dict[ref_name] = ref_total_read_count*read_length / ref_total_len


    # Loop over all the -rpkm.tsv directories
    ref_breath_dict={}
    for rpkm_file in os.listdir(quant_dir): # Loop over the reference alignment files
        if rpkm_file.endswith('-rpkm.tsv'):
            ref_name = parse_file_name(rpkm_file)
            ref_bp_covered = 0
            ref_len = 0
            with open(os.path.join(quant_dir, rpkm_file) , 'r') as f:
                for line in f:
                    line_tok_l1 = line.split('\t')
                    ref_contig = line_tok_l1[0]
                    ref_contig_breath = float(line_tok_l1[1])
                    ref_contig_len = ref_len_dict[ref_name][ref_contig]
                    ref_bp_covered = ref_bp_covered +(ref_contig_breath*ref_contig_len)
                    ref_len = ref_len + ref_contig_len

            ref_bp_covered = ref_bp_covered / ref_len
            ref_breath_dict[ref_name] = ref_bp_covered
    return ref_breath_dict, ref_coverage_dict


if __name__ == '__main__':
    ref_genomes_dir = sys.argv[1]
    base_dir = sys.argv[2]

    # Get list of references
    references=[]
    for ref_file in os.listdir(ref_genomes_dir):
        references.append(os.path.splitext(ref_file)[0])

    ref_sample_breath = {}
    ref_sample_depth = {}
    for sample_dir in os.listdir(base_dir):
        sample_name = sample_dir
        sample_dir = os.path.join(base_dir,sample_name)
        if sample_name.startswith('synthetic_1_') and os.path.isdir(sample_dir):
            ref_breath_dict, ref_coverage_dict = gather_reference_breath(ref_genomes_dir, sample_name, base_dir)
            ref_sample_breath[sample_name] = ref_breath_dict
            ref_sample_depth[sample_name] = ref_coverage_dict


    # Print summary of the coverage breath
    f = open(os.path.join(base_dir, 'sample_read_breath.csv'),'w')
    f.write('reference_genome,')
    for ref in references:
        f.write(ref+',')
    f.write('\n')
    for sample_name, sample_refs in ref_sample_breath.items():
        f.write(sample_name + ',')
        for ref in references:
            if ref in sample_refs:
                f.write(str(sample_refs[ref])+',')
            else:
                f.write('0,')
        f.write('\n')
    f.close()

    # Print summary of the coverage depth
    f = open(os.path.join(base_dir, 'sample_read_depth.csv'),'w')
    f.write('reference_genome,')
    for ref in references:
        f.write(ref+',')
    f.write('\n')
    for sample_name, sample_refs in ref_sample_depth.items():
        f.write(sample_name + ',')
        for ref in references:
            if ref in sample_refs:
                f.write(str(sample_refs[ref])+',')
            else:
                f.write('0,')
        f.write('\n')
    f.close()
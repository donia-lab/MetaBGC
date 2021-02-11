#!/usr/bin/env python3


from Bio import SeqIO
from glob import glob
import numpy as np
import os
import re
import json
import subprocess
from multiprocessing import Pool
from functools import partial
import config
import argparse


def get_records(outdir, organism):
    """
    Split an input fasta into separate fasta files
    
    :param organism: the fasta file of an organism
    :returns: a tuple of fasta names and their lengths
    """
    record_names = []
    record_len = []
    with open(organism, "r") as h:
        for record in SeqIO.parse(h, "fasta"):
            f_name = "{}/{}.fasta".format(outdir, 
                     re.sub(r".*\|(.*)\|$", r"\1", record.name))
            with open(f_name, "w") as f:
                SeqIO.write(record, f, "fasta")
            record_names.append(f_name)
            record_len.append(len(record))
    return record_names, record_len


def assign_reads(seq_records_background, seq_records_PKS,
                 args, sample_no):
    """
    Assign number of reads per organism

    :param seq_records: a list of records returned by get_records
    :param sample_no: the number for the synthetic sample
    :returns: a dictionary of form fasta_file:number of reads
    """
    np.random.seed(args.seed + sample_no)

    # select background genomes
    weight = [1] * int(len(seq_records_background) * 0.4) + \
             [0.7] * int(len(seq_records_background) * 0.2) + \
             [0.5] * int(len(seq_records_background) * 0.2) + \
             [0.3] * int(len(seq_records_background) * 0.2) + \
             [0.1] * (len(seq_records_background) - \
                      int(len(seq_records_background) * 0.4) \
             - 3 * int(len(seq_records_background) * 0.2))
    selected_idx = np.random.choice(range(len(seq_records_background)),
                                    int(args.prop * \
                                        len(seq_records_background)),
                                    replace=False,
                                    p=np.divide(weight, sum(weight)))
    selected_records = [seq_records_background[i] for i in selected_idx]

    # select PKS-containing genomes
    selected_size = np.random.choice([0, 1, 2, 3], 1)
    if selected_size > 0:
        # select genomes
        selected_idx = np.random.choice(range(len(seq_records_PKS)),
                                        selected_size,
                                        replace=False)
        selected_records_PKS = [seq_records_PKS[i] for i in selected_idx]
        selected_records = selected_records + selected_records_PKS

    # sample abundance for each selected genome with log-normal distribution
    abund_raw = np.random.lognormal(mean=1, sigma=2,
                                    size=len(selected_records))
    # number of reads per organism
    num_reads = abund_raw / abund_raw.sum() * args.num_reads

    # generate number of reads per fasta file
    num_dict = {}
    for idx1, num_read in enumerate(num_reads):
        arr = np.array(selected_records[idx1][1])
        totlen = arr.sum()
        # number of reads per scaffold
        for idx2, record_name in enumerate(selected_records[idx1][0]):
            num_dict[record_name] = int(arr[idx2] / totlen * num_read)

    # save simulation profile
    sample_dir = "{}/{}{}".format(args.outdir,
                                  args.base_name,
                                  sample_no)
    try:
        os.mkdir(sample_dir)
    except:
        pass

    return num_dict 

   
def call_ART(args, num_dict, ART_path, sample_no):
    """
    Call ART program
 
    :param num_dict: dictionary generated from assign_reads
    :returns: None
    """
    for key in num_dict:
        with open(os.devnull, "w") as f:
            subprocess.call("{} -ss {} -i {} -p -q -l {} " \
                            "-c {} -m {} -s {} -o {} -na"\
                            .format(ART_path, args.system, key,
                                    args.length, 
                                    num_dict[key], args.mflen,
                                    args.mflensd, "{}_{}{}-".format(\
                                    re.sub(r"(.*)\..*", r"\1", key),
                                    args.base_name, sample_no)), 
                            shell=True, stdout=f)
    
    read1_to_cat = ["{}_{}{}-1.fq".format(re.sub(r"(.*)\..*", r"\1", 
                                                 key),
                                          args.base_name, sample_no) \
                    for key in num_dict]
    read2_to_cat = ["{}_{}{}-2.fq".format(re.sub(r"(.*)\..*", r"\1",
                                                 key),
                                          args.base_name, sample_no) \
                    for key in num_dict]


    sample_dir = "{}/{}{}".format(args.outdir, args.base_name, sample_no)

    try:        
        os.mkdir("{}/{}{}-raw-reads-fastq".format(sample_dir, 
                                                  args.base_name, 
                                                  sample_no))

    except:
        pass

    finally:
        for idx in range(len(read1_to_cat)):
            subprocess.call("cat {0} >> {1}/{2}{3}-raw-reads-fastq/" \
                            "{2}{3}-1.fastq".format(read1_to_cat[idx], 
                                                    sample_dir, 
                                                    args.base_name,
                                                    sample_no),
                            shell=True)
            subprocess.call("cat {0} >> {1}/{2}{3}-raw-reads-fastq/" \
                            "{2}{3}-2.fastq".format(read2_to_cat[idx],
                                                    sample_dir,
                                                    args.base_name,
                                                    sample_no),
                            shell=True)
            os.remove(read1_to_cat[idx])
            os.remove(read2_to_cat[idx])


def generate_metagenomic_sample(sample_no,seq_records_background,
                                seq_records_PKS, args, ART_path):
    """
    Generate a metagenomic sample

    """
    num_dict = assign_reads(seq_records_background, seq_records_PKS,
                            args, sample_no)

    sample_dir = "{}/{}{}".format(args.outdir, args.base_name, 
                                  sample_no)

    with open("{}/{}{}_abundance.json".format(sample_dir, 
                                              args.base_name,
                                              sample_no), "w") as h:
        json.dump(num_dict, h, indent=4)

    call_ART(args, num_dict, ART_path, sample_no)


def main(args):
    tmp_dir1 = args.outdir + "/tmp1"
    tmp_dir2 = args.outdir + "/tmp2"
    try:
        os.mkdir(args.outdir)
    except:
        pass

    try:
        os.mkdir(tmp_dir1)
    except:
        pass

    try:
        os.mkdir(tmp_dir2)
    except:
        pass

    ART_path = config.ART_path

    background_fasta = glob(args.indir1 + "/*.fa") + \
                      glob(args.indir1 + "/*.fna") + \
                      glob(args.indir1 + "/*.fasta") + \
                      glob(args.indir1 + "/*.fsa") + \
                      glob(args.indir1 + "/*.fas")

    PKS_fasta = glob(args.indir2 + "/*.fa") + \
                glob(args.indir2 + "/*.fna") + \
                glob(args.indir2 + "/*.fasta") + \
                glob(args.indir2 + "/*.fsa") + \
                glob(args.indir2 + "/*.fas")

    seq_records_background = [get_records(tmp_dir1, i) for i in background_fasta]
    seq_records_PKS = [get_records(tmp_dir2, i) for i in PKS_fasta]
    
    with Pool(processes=args.processes) as pool:
        comp_list = pool.map(partial(generate_metagenomic_sample, 
                                     seq_records_background=\
                                     seq_records_background,
                                     seq_records_PKS=seq_records_PKS,
                                     args=args,
                                     ART_path=ART_path), 
                             range(args.samples), chunksize=1)
    print("Synthesis complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=\
        "synthesize metagenomic samples with Illumina sequencing")
    parser.add_argument("-i1", "--indir1", type=str, required=True, 
                        help="input directory of background fasta " \
                             "files for simulation")
    parser.add_argument("-i2", "--indir2", type=str, required=True,
                        help="input directory of PKS fasta files " \
                             "for simulation")
    parser.add_argument("-o", "--outdir", type=str, required=True, 
                        help="output directory for simulated samples")
    parser.add_argument("-ss", "--system", type=str, required=True, 
                        help="Illumina sequencing system " \
                             "(HS10, HS25, MSv3, etc.). Same " \
                             "options as -ss in art_illumina")
    parser.add_argument("-l", "--length", type=int, required=True, 
                        help="read length in bp")
    parser.add_argument("-fl", "--mflen", type=int, required=True, 
                        help="mean fragment size in bp")
    parser.add_argument("-sd", "--mflensd", type=int, required=True, 
                        help="standard dev of fragment size in bp")
    parser.add_argument("-nr", "--num_reads", type=int, required=True, 
                        help="total number of read pairs") 
    parser.add_argument("-ns", "--samples", type=int, required=True, 
                        help="number of samples to generate")
    parser.add_argument("-p", "--prop", type=float, required=True, 
                        help="proportion of organisms to draw for " \
                             "each metagenomic sample. Should be " \
                             "between 0 and 1")
    parser.add_argument("-b", "--base_name", type=str, default="S",
                        help="prefix of sample name. Default is 'S'")
    parser.add_argument("-t", "--processes", type=int, default=1, 
                        help="number of processes to run. Default is " \
                             "1. Should be smaller than --samples")
    parser.add_argument("-rs", "--seed", type=int, default=915,
                        help="random seed")
    args = parser.parse_args()
    main(args) 

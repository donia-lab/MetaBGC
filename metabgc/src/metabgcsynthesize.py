from Bio import SeqIO
from glob import glob
import numpy as np
import os
import re
import json
import subprocess
from multiprocessing import Pool
from functools import partial

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


def assign_reads(seq_records_background, seq_records_PKS, sample_no, seed, prop, num_reads, outdir, base_name):
    """
    Assign number of reads per organism

    :param seq_records: a list of records returned by get_records
    :param sample_no: the number for the synthetic sample
    :returns: a dictionary of form fasta_file:number of reads
    """
    np.random.seed(seed + sample_no)

    # select background genomes
    weight = [1] * int(len(seq_records_background) * 0.4) + \
             [0.7] * int(len(seq_records_background) * 0.2) + \
             [0.5] * int(len(seq_records_background) * 0.2) + \
             [0.3] * int(len(seq_records_background) * 0.2) + \
             [0.1] * (len(seq_records_background) - int(len(seq_records_background) * 0.4) - 3 * int(len(seq_records_background) * 0.2))
    selected_idx = np.random.choice(range(len(seq_records_background)),
                                    int(prop * len(seq_records_background)),
                                    replace=False,
                                    p =np.divide(weight, sum(weight)))
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
    num_reads = abund_raw / abund_raw.sum() * num_reads

    # generate number of reads per fasta file
    num_dict = {}
    for idx1, num_read in enumerate(num_reads):
        arr = np.array(selected_records[idx1][1])
        totlen = arr.sum()
        # number of reads per scaffold
        for idx2, record_name in enumerate(selected_records[idx1][0]):
            num_dict[record_name] = int(arr[idx2] / totlen * num_read)

    # save simulation profile
    sample_dir = "{}/{}{}".format(outdir, base_name, sample_no)
    try:
        os.mkdir(sample_dir)
    except:
        pass
    return num_dict

def call_ART(num_dict, sample_no, system, length, mflen, mflensd, base_name, outdir):
    """
    Call ART program
 
    :param num_dict: dictionary generated from assign_reads
    :returns: None
    """
    for key in num_dict:
        with open(os.devnull, "w") as f:
            subprocess.call("{} -ss {} -i {} -p -q -l {} " \
                            "-c {} -m {} -s {} -o {} -na"\
                            .format('art_illumina', system, key,
                                    length, num_dict[key], mflen,
                                    mflensd, "{}_{}{}-".format(\
                                    re.sub(r"(.*)\..*", r"\1", key),
                                    base_name, sample_no)),
                            shell=True, stdout=f)
    
    read1_to_cat = ["{}_{}{}-1.fq".format(re.sub(r"(.*)\..*", r"\1", 
                                                 key),
                                          base_name, sample_no) \
                    for key in num_dict]
    read2_to_cat = ["{}_{}{}-2.fq".format(re.sub(r"(.*)\..*", r"\1",
                                                 key),
                                          base_name, sample_no) \
                    for key in num_dict]


    sample_dir = "{}/{}{}".format(outdir, base_name, sample_no)

    try:        
        os.mkdir("{}/{}{}-raw-reads-fastq".format(sample_dir, 
                                                  base_name,
                                                  sample_no))

    except:
        pass

    finally:
        for idx in range(len(read1_to_cat)):
            subprocess.call("cat {0} >> {1}/{2}{3}-raw-reads-fastq/" \
                            "{2}{3}-1.fastq".format(read1_to_cat[idx], 
                                                    sample_dir, 
                                                    base_name,
                                                    sample_no),
                            shell=True)
            subprocess.call("cat {0} >> {1}/{2}{3}-raw-reads-fastq/" \
                            "{2}{3}-2.fastq".format(read2_to_cat[idx],
                                                    sample_dir,
                                                    base_name,
                                                    sample_no),
                            shell=True)
            os.remove(read1_to_cat[idx])
            os.remove(read2_to_cat[idx])


def generate_metagenomic_sample(sample_no,seq_records_background,
                                seq_records_PKS, seed, prop, num_reads, system, length, mflen, mflensd, base_name, outdir):
    """
    Generate a metagenomic sample

    """
    num_dict = assign_reads(seq_records_background, seq_records_PKS, sample_no, seed, prop, num_reads, outdir, base_name)

    sample_dir = "{}/{}{}".format(outdir, base_name, sample_no)

    with open("{}/{}{}_abundance.json".format(sample_dir, base_name, sample_no), "w") as h:
        json.dump(num_dict, h, indent=4)

    with open("{}/{}{}_abundance.csv".format(sample_dir, base_name, sample_no), "w") as h:
        for key in num_dict:
            h.write('{},{},{},{}\n'.format(sample_no, os.path.basename(key), str(num_dict[key]), str(num_dict[key]/num_reads)))
    call_ART(num_dict, sample_no, system, length, mflen, mflensd, base_name, outdir)

def mbgcsynthesize(indir1, indir2, system, length, mflen, mflensd, num_reads,
               samples, prop, output_directory, base_name, cpu, seed):
    tmp_dir1 = output_directory + "/tmp1"
    tmp_dir2 = output_directory + "/tmp2"
    try:
        os.mkdir(output_directory)
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

    background_fasta = glob(indir1 + "/*.fa") + \
                      glob(indir1 + "/*.fna") + \
                      glob(indir1 + "/*.fasta") + \
                      glob(indir1 + "/*.fsa") + \
                      glob(indir1 + "/*.fas")

    PKS_fasta = glob(indir2 + "/*.fa") + \
                glob(indir2 + "/*.fna") + \
                glob(indir2 + "/*.fasta") + \
                glob(indir2 + "/*.fsa") + \
                glob(indir2 + "/*.fas")

    seq_records_background = [get_records(tmp_dir1, i) for i in background_fasta]
    seq_records_PKS = [get_records(tmp_dir2, i) for i in PKS_fasta]
    
    with Pool(processes=cpu) as pool:
        comp_list = pool.map(partial(generate_metagenomic_sample, 
                                     seq_records_background=seq_records_background,
                                     seq_records_PKS=seq_records_PKS,
                                     seed=seed, prop=prop, num_reads=num_reads,
                                     system=system, length=length, mflen=mflen,
                                     mflensd=mflensd, base_name=base_name, outdir=output_directory),
                             range(samples), chunksize=1)

    f_agg = open(output_directory + os.sep +'sample_abund.csv', 'w')
    f_agg.write('Sample,genome,num_reads,abund_ratio\n')
    for sample_no in range(0, samples):
        sample_dir = "{}/{}{}".format(output_directory, base_name, sample_no)
        with open("{}/{}{}_abundance.csv".format(sample_dir, base_name, sample_no)) as infile:
            for line in infile:
                f_agg.write(line)
    f_agg.close()
    print("Synthesis complete!")
#!/usr/bin/python

from __future__ import division
from Bio import AlignIO
import os 
from CreateSpHMMs import 

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--aln_file', type= str,required=True) 
    parser.add_argument('--window_len',type= int, required=False, default=10)
    parser.add_argument('--kmer_len', type= int, required=True ) 
    parser.add_argument('--outdir', type= str, required=True)
    parser.add_argument('--hmmName', type= str, required=True)
    parser.add_argument('--start', type = int, required=False)
    parser.add_argument('--end', type = int, required=False)        
    args = parser.parse_args()
    GenerateSpHMM(args.aln_file, args.window_len, args.kmer_len, args.outdir, args.hmmName, args.start, args.end )



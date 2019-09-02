import os, re, time, sys, string, math
import argparse
import ntpath
from Bio import AlignIO
from Bio import SeqIO
from Utils.Utils import RunHMMDirectory
from Utils.Utils import RunBLASTNDirectory
from Utils.Utils import runMakeBLASTDB
from Utils.Utils import runBLASTN
from Utils.Utils import runMUSCLE
from Utils.Utils import parseHMM
from Utils.Utils import createPandaDF
from Utils.Utils import runTranSeq
from Utils.Utils import PreProcessReads
from CreateSpHMMs import GenerateSpHMM
from rpy2.robjects.packages import STAP
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from fuzzywuzzy import fuzz
from shutil import copyfile

CPU_THREADS = 4

if __name__ == '__main__':
    startTime = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--prot_alignment', required=True, help="Alignment of the protein homologs in FASTA format.")
    parser.add_argument('--prot_family_name', required=True, help="Name of the protein family.")
    parser.add_argument('--cohort_name', required=True, help="Name of the sample/cohort.")
    parser.add_argument('--nucl_seq_directory', required=True, help="Directory with synthetic read files of the cohort.")
    parser.add_argument('--seq_fmt', required=True, help="Sequence file format and extension.: {fasta,fastq}.")
    parser.add_argument('--pair_fmt', required=True, help="Sequence pair format: {single, split, interleaved}.")
    parser.add_argument('--R1_file_suffix', required=False, help="Suffix including extension of the file name specifying R1 reads. Not specified for single or interleaved reads.")
    parser.add_argument('--R2_file_suffix', required=False, help="Suffix including extension of the file name specifying R2 reads. Not specified for single or interleaved reads.")
    parser.add_argument('--tp_genes_nucl', required=True, help="Multi-FASTA with the nucleotide sequence of the true positive genes.")
    parser.add_argument('--F1_Thresh', required=True, help="F1 score threshold.")
    parser.add_argument('--output_directory', required=True, help="Directory to save results.")
    parser.add_argument('--cpu', required=False, help="Number of threads. Def.: 4")
    args = parser.parse_args()

    if args.cpu is not None:
        CPU_THREADS = args.cpu

    # Gen spHMMs and interval pos
    hmm_directory = os.path.join(args.output_directory, 'spHMMs')
    os.makedirs(hmm_directory,0o777,True)
    copyfile(args.prot_alignment, hmm_directory+os.sep+ntpath.basename(args.prot_alignment))
    prot_aln_file = os.path.join(hmm_directory,ntpath.basename(args.prot_alignment))
    alignment = AlignIO.read(args.prot_alignment, "fasta")
    hmmDict = GenerateSpHMM(prot_aln_file, 10, 30, hmm_directory, args.prot_family_name, 1, alignment.get_alignment_length()+1)

    tp_genes_prot = args.output_directory+os.sep+"TPGenes.faa"
    runTranSeq(args.tp_genes_nucl,"1",tp_genes_prot)
    joinFilenames = [tp_genes_prot, prot_aln_file]
    tmpFile = os.path.join(args.output_directory,"tmp.fa")
    with open(tmpFile, 'w') as outfile:
        for fname in joinFilenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    alnOutput = os.path.join(args.output_directory,"tmp.afa")
    runMUSCLE(tmpFile, alnOutput)
    #gene_pos_file = os.path.join(args.output_directory, 'Gene_Interval_Pos.txt')
    # outfile = open(gene_pos_file, 'w')
    # outfile.write("gene_name\tstart\tend\tinterval\tcyclase_type\n")
    # for record in SeqIO.parse(tp_genes_prot, "fasta"):
    #     for hmmInterval, hmmFile in hmmDict.items():
    #         seqFile = hmmFile.split('.hmm')[0] +".fas"
    #         hmmIntervalSeqs = list(SeqIO.parse(seqFile, "fasta"))
    #         tmpFile = os.path.join(args.output_directory,"tmp.fa")
    #         alnOutput = os.path.join(args.output_directory,"tmp.afa")
    #         alnRecord = [record]
    #         for hmmRec in hmmIntervalSeqs:
    #             alnRecord.append(hmmRec)
    #         SeqIO.write(alnRecord, tmpFile, "fasta")
    #         runMUSCLE(tmpFile,alnOutput)
    #         alnSeqs = list(SeqIO.parse(alnOutput, "fasta"))
    #         foundCtr=0
    #         startPos = alnLen = -1
    #         for hmmRec in hmmIntervalSeqs:
    #             for alnRec in alnSeqs:
    #                 gapStripSeq = alnRec.seq.strip("-")
    #                 fuzzyScore = fuzz.token_set_ratio(hmmRec.seq, gapStripSeq)
    #                 if fuzzyScore > 70 :
    #                     foundCtr += 1
    #                     startPos = alnRec.seq.find(gapStripSeq)
    #                     if alnLen < len(hmmRec):
    #                         alnLen = len(hmmRec)
    #                     break;
    #         if foundCtr == len(hmmIntervalSeqs):
    #             startPos = startPos*3
    #             endPos = startPos + (alnLen * 3)
    #             outfile.write(record.id+"\t"+str(startPos)+"\t"+str(endPos)+"\t"+hmmInterval+"\t"+args.prot_family_name+"\n")
    # outfile.close()

    #Preprocess synthetic reads
    nucl_seq_directory = PreProcessReads(args.nucl_seq_directory,args.seq_fmt,args.pair_fmt,args.R1_file_suffix.strip(),args.R2_file_suffix.strip(),args.output_directory)
    # Translate nucleotide seq
    prot_seq_directory = os.path.join(args.output_directory, 'prot_seq_dir')
    os.makedirs(prot_seq_directory, 0o777, True)
    for subdir, dirs, files in os.walk(nucl_seq_directory):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*\.fasta$", file) and os.path.getsize(filePath) > 0:
                prot_file = prot_seq_directory + os.sep + ntpath.basename(filePath)
                runTranSeq(filePath, "6", prot_file)

    # HMMER Search
    hmm_search_directory = os.path.join(args.output_directory, 'hmm_result')
    os.makedirs(hmm_search_directory,0o777,True)
    for hmmInterval, hmmFile in hmmDict.items():
        RunHMMDirectory(prot_seq_directory,hmmFile, args.cohort_name, args.prot_family_name, "30_10", hmmInterval, hmm_search_directory, CPU_THREADS)
    allHMMResult = hmm_search_directory + os.sep + "CombinedHmmSearch.txt"
    with open(allHMMResult, 'w') as outfile:
        for subdir, dirs, files in os.walk(hmm_search_directory):
            for file in files:
                filePath = os.path.join(subdir, file)
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)

    # BLAST Alignment
    blastn_search_directory = os.path.join(args.output_directory, 'blastn_result')
    os.makedirs(blastn_search_directory,0o777,True)
    RunBLASTNDirectory(nucl_seq_directory, args.tp_genes_nucl, blastn_search_directory,CPU_THREADS)
    allBLASTResult = blastn_search_directory + os.sep + "CombinedBLASTSearch.txt"
    with open(allBLASTResult, 'w') as outfile:
        outfile.write("sseqid\tslen\tsstart\tsend\tqseqid\tqlen\tqstart\tqend\tpident\tevalue\tSample\tsampleType\n")
        for subdir, dirs, files in os.walk(blastn_search_directory):
            for file in files:
                filePath = os.path.join(subdir, file)
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            sampleName = ntpath.basename(filePath).split(".txt")[0]
                            outfile.write(line.strip() + "\t" + sampleName + "\t" + args.cohort_name + "\n")

    # Eval spHMMs
    #rpackages.importr('base')
    #packageNames = ('tidyverse','ggsci','ggpubr')
    #utils = rpackages.importr('utils')
    #utils.chooseCRANmirror(ind=1)
    #packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
    #if len(packnames_to_install) > 0:
    #    utils.install_packages(StrVector(packnames_to_install))

    rpackages.importr('tidyverse')
    rpackages.importr('ggsci')
    rpackages.importr('ggpubr')

    hp_hmm_directory = os.path.join(args.output_directory, 'HiPer_spHMMs')
    os.makedirs(hp_hmm_directory,0o777,True)
    with open('EvaluateSpHMMs.R', 'r') as f:
        rStr = f.read()
    myfunc = STAP(rStr, "EvaluateSpHMM")
    myfunc.EvaluateSpHMM(allHMMResult, allBLASTResult, gene_pos_file, args.prot_family_name, float(args.F1_Thresh), hmm_directory, hp_hmm_directory)

    timeTaken = time.time() - startTime
    mins = int(timeTaken / 60)
    secs = int(timeTaken) % 60
    print("\nTotal time taken : " + str(mins) + " mins " + str(secs) + " seconds")

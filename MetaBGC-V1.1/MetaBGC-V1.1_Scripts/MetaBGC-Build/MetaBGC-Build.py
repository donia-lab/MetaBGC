import os, re, time, sys, string, math
import argparse
import ntpath
from Bio import AlignIO
from Bio import SeqIO
from Utils.Utils import RunHMMDirectory
from Utils.Utils import RunBLASTNDirectoryPar
from Utils.Utils import runMakeBLASTDB
from Utils.Utils import runBLASTN
from Utils.Utils import runMUSCLE
from Utils.Utils import parseHMM
from Utils.Utils import createPandaDF
from Utils.Utils import runTranSeq
from Utils.Utils import TranseqReadsDir
from Utils.Utils import PreProcessReadsPar
from CreateSpHMMs import GenerateSpHMM
from rpy2.robjects.packages import STAP
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
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
        CPU_THREADS = int(args.cpu)

    build_op_dir = args.output_directory + os.sep + "build"

    # Gen spHMMs and interval pos
    hmm_directory = os.path.join(build_op_dir, 'spHMMs')
    os.makedirs(hmm_directory,0o777,True)
    copyfile(args.prot_alignment, hmm_directory+os.sep+ntpath.basename(args.prot_alignment))
    prot_aln_file = os.path.join(hmm_directory,ntpath.basename(args.prot_alignment))
    alignment = AlignIO.read(args.prot_alignment, "fasta")
    hmmDict = GenerateSpHMM(prot_aln_file, 10, 30, hmm_directory, args.prot_family_name, 1, alignment.get_alignment_length()+1)

    tp_genes_prot = build_op_dir+os.sep+"TPGenes.faa"
    runTranSeq(args.tp_genes_nucl,"1",tp_genes_prot)
    tmpFile = os.path.join(build_op_dir,"tmp.fa")

    # Join true positives in the sample with the BGC proteins
    joinedSeqs = []
    tpGeneSeqs = list(SeqIO.parse(tp_genes_prot, "fasta"))
    for seq in tpGeneSeqs:
        seq.id = seq.id[:-2]
        seq.description = ""
        joinedSeqs.append(seq)
    protAlnSeqs = list(SeqIO.parse(prot_aln_file, "fasta"))
    for seq in protAlnSeqs:
        joinedSeqs.append(seq)
    SeqIO.write(joinedSeqs, tmpFile, "fasta")

    # MUSCLE align TP genes with markers
    alnOutput = os.path.join(build_op_dir,"tmp.afa")
    runMUSCLE(tmpFile, alnOutput)
    muscleAlnSeqs = list(SeqIO.parse(alnOutput, "fasta"))
    protPosList = []
    for i,pseq in enumerate(protAlnSeqs):
        for j,mseq in enumerate(muscleAlnSeqs):
            if(pseq.id==mseq.id):
                protPosList.append(j)

    # Extract spHMM coordinates from MUSCLE alignment
    gene_pos_file = os.path.join(build_op_dir, 'Gene_Interval_Pos.txt')
    outfile = open(gene_pos_file, 'w')
    outfile.write("gene_name\tstart\tend\tinterval\tcyclase_type\n")
    for i, mseq in enumerate(muscleAlnSeqs):
        if i not in protPosList:
            protPos = min(protPosList, key=lambda x: abs(x - i))
            protSeq = str(protAlnSeqs[protPosList.index(protPos)].seq)
            protSeqNoGap = protSeq.replace('-', '')
            protMuscleSeq = str(muscleAlnSeqs[protPos].seq)
            ungappedCoords = [0] * len(protMuscleSeq)
            aaIndex = 0
            for i,aa in enumerate(protMuscleSeq):
                if aa=='-':
                    ungappedCoords[i]=-1
                else:
                    ungappedCoords[i] = aaIndex
                    aaIndex=aaIndex+1
            for key in hmmDict:
                startPos = int(key.split('_')[0])
                endPos = int(key.split('_')[1])
                windowSeq = protSeq[startPos:endPos]
                windowSeqNoGap = windowSeq.replace('-', '')
                beginGaps = 0
                endGaps = 0
                isEndGap = 0
                for aaVal in windowSeq:
                    if isEndGap==1 and aaVal=='-':
                        endGaps= endGaps +1
                    elif isEndGap == 1 and aaVal != '-':
                        endGaps = 0
                    else:
                        if aaVal!='-':
                            isEndGap = 1
                        else:
                            beginGaps = beginGaps + 1
                startPosNoGap = protSeqNoGap.find(windowSeqNoGap)
                endPosNoGap = startPosNoGap + len(windowSeqNoGap)

                startPosGapped = ungappedCoords.index(startPosNoGap)
                endPosGapped = ungappedCoords.index(endPosNoGap)

                #Add begin end gaps
                startPosGapped = 0 if (startPosGapped-beginGaps)<0 else (startPosGapped-beginGaps)
                endPosGapped = len(protMuscleSeq) if (endPosGapped + endGaps) >= len(protMuscleSeq) else (endPosGapped + endGaps)

                startPosGapped = startPosGapped * 3
                endPosGapped = endPosGapped * 3
                outfile.write(mseq.id+"\t"+str(startPosGapped)+"\t"+str(endPosGapped)+"\t"+key+"\t"+args.prot_family_name+"\n")
    outfile.close()

    #Preprocess synthetic reads
    nucl_seq_directory = PreProcessReadsPar(args.nucl_seq_directory,
                                            args.seq_fmt,args.pair_fmt,
                                            args.R1_file_suffix.strip(),
                                            args.R2_file_suffix.strip(),
                                            build_op_dir,
                                            CPU_THREADS)
    # Translate nucleotide seq
    prot_seq_directory = TranseqReadsDir(build_op_dir, nucl_seq_directory, CPU_THREADS)

    # HMMER Search
    hmm_search_directory = os.path.join(build_op_dir, 'hmm_result')
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
    blastn_search_directory = os.path.join(build_op_dir, 'blastn_result')
    os.makedirs(blastn_search_directory,0o777,True)
    RunBLASTNDirectoryPar(nucl_seq_directory, args.tp_genes_nucl, blastn_search_directory,CPU_THREADS)
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
    rpackages.importr('base')
    #packageNames = ('tidyverse','ggsci','ggpubr')
    #utils = rpackages.importr('utils')
    #utils.chooseCRANmirror(ind=1)
    #packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
    #if len(packnames_to_install) > 0:
    #    utils.install_packages(StrVector(packnames_to_install))

    rpackages.importr('tidyverse')
    rpackages.importr('ggsci')
    rpackages.importr('ggpubr')

    hp_hmm_directory = os.path.join(build_op_dir, 'HiPer_spHMMs')
    os.makedirs(hp_hmm_directory,0o777,True)
    with open('EvaluateSpHMMs.R', 'r') as f:
        rStr = f.read()
    myfunc = STAP(rStr, "EvaluateSpHMM")
    myfunc.EvaluateSpHMM(allHMMResult, allBLASTResult, gene_pos_file, args.prot_family_name, float(args.F1_Thresh), hmm_directory, hp_hmm_directory)

    timeTaken = time.time() - startTime
    mins = int(timeTaken / 60)
    secs = int(timeTaken) % 60
    print("\nTotal time taken : " + str(mins) + " mins " + str(secs) + " seconds")

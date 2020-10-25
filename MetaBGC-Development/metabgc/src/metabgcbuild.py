import time, sys
from Bio import AlignIO
from Bio import SeqIO
from metabgc.src.utils import *
from metabgc.src.blastrunlib import *
import metabgc.src.createsphmms as createhmm
from rpy2.robjects.packages import STAP
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from shutil import copyfile
import os

def ungappedseqsearch(reference_str, query_str):
    for i in range(0, len(reference_str)):
        if reference_str[i]=='-':
            continue
        k = i
        j = 0
        while j < len(query_str) and k < len(reference_str):
            if query_str[j] == '-':
                j = j + 1
            elif reference_str[k] == '-':
                k = k + 1
            elif reference_str[k] != query_str[j]:
                break
            else:
                k = k + 1
                j = j + 1
        if j == len(query_str):
            return [i, k]
    return [-1, -1]

def gensphmmfiles(prot_family_name,prot_aln_file,tp_prot_file,hmm_directory,gene_pos_file,gene_pos_file_aa):
    alignment = AlignIO.read(prot_aln_file, "fasta")
    hmmDict = createhmm.GenerateSpHMM(prot_aln_file, tp_prot_file,10, 30, hmm_directory, prot_family_name, 1, alignment.get_alignment_length()+1, gene_pos_file,gene_pos_file_aa)
    return hmmDict

def gengeneposlist(prot_family_name,protAlnSeqs,hmmDict,alnOutput,gene_pos_file):
    muscleAlnSeqs = list(SeqIO.parse(alnOutput, "fasta"))
    protPosList = []
    for i,pseq in enumerate(protAlnSeqs):
        for j,mseq in enumerate(muscleAlnSeqs):
            if(pseq.id==mseq.id):
                protPosList.append(j)
    outfile = open(gene_pos_file, 'w')
    outfile.write("gene_name\tstart\tend\tinterval\tprot_type\n")
    for i, mseq in enumerate(muscleAlnSeqs):
        if i not in protPosList:
            protPos = min(protPosList, key=lambda x: abs(x - i))
            protSeq = str(protAlnSeqs[protPosList.index(protPos)].seq)
            protMuscleSeq = str(muscleAlnSeqs[protPos].seq)
            protTPSeq = str(mseq.seq).replace('-', '')
            for key in hmmDict:
                startPos = int(key.split('_')[0])
                endPos = int(key.split('_')[1])
                windowSeq = protSeq[startPos:endPos]
                coordMUSCLEList = ungappedseqsearch(protMuscleSeq,windowSeq)

                startTPPosGapped = coordMUSCLEList[0]
                endTPPosGapped = coordMUSCLEList[1]
                windowTPSeqGapped = mseq.seq[startTPPosGapped:endTPPosGapped]

                coordTPList = ungappedseqsearch(protTPSeq, windowTPSeqGapped)
                startPosGapped = coordTPList[0]
                endPosGapped = coordTPList[1]

                startPosGapped = startPosGapped * 3
                endPosGapped = endPosGapped * 3
                outfile.write(mseq.id+"\t"+str(startPosGapped)+"\t"+str(endPosGapped)+"\t"+
                              str(hmmDict[key].intervalStart)+"_"+str(hmmDict[key].intervalEnd)+
                              "\t"+prot_family_name+"\n")
    outfile.close()

def mbgcbuild(prot_alignment, prot_family_name, cohort_name,
              nucl_seq_directory, prot_seq_directory, seq_fmt, pair_fmt, r1_file_suffix,
              r2_file_suffix, tp_genes_nucl, blast_db_directory_map_file,
              blastn_search_directory, hmm_search_directory, f1_thresh, output_directory, cpu):
    try:
        CPU_THREADS = 4
        startTime = time.time()
        if cpu is not None:
            CPU_THREADS = int(cpu)

        # setup paths
        build_op_dir = output_directory + os.sep + "build"
        hmm_directory = os.path.join(build_op_dir, 'spHMMs')
        tp_genes_prot = build_op_dir + os.sep + "TPGenes.faa"
        alnOutput = os.path.join(build_op_dir,"TP_Homolog_Alignment.afa")
        gene_pos_file = os.path.join(build_op_dir, 'Gene_Interval_Pos.txt')
        gene_pos_file_aa = os.path.join(build_op_dir, 'Gene_Interval_Pos_AA.txt')
        if hmm_search_directory is None:
            hmm_search_directory = os.path.join(build_op_dir, 'hmm_result')
        allHMMResult = os.path.join(build_op_dir,"CombinedHmmSearch.txt")
        if blastn_search_directory is None:
            blastn_search_directory = os.path.join(build_op_dir, 'blastn_result')
        allBLASTResult = os.path.join(build_op_dir,"CombinedBLASTSearch.txt")
        if prot_seq_directory is None:
            prot_seq_directory = ""

        # Create OP dirs
        os.makedirs(hmm_directory, 0o777, True)

        joinedSeqs = []
        tmpFile = os.path.join(build_op_dir, "TP_Homolog.faa")
        # Translate protein sequence
        if os.path.getsize(tp_genes_nucl) > 0:
            runTranSeq(tp_genes_nucl,"1",tp_genes_prot)
            # Join true positives in the sample with the BGC proteins
            tpGeneSeqs = list(SeqIO.parse(tp_genes_prot, "fasta"))
            # Removing _1 added by TranSeq
            for seq in tpGeneSeqs:
                seq.id = seq.id[:-2]
                seq.description = ""
                joinedSeqs.append(seq)
            SeqIO.write(joinedSeqs,tp_genes_prot,"fasta")

        protAlnSeqs = list(SeqIO.parse(prot_alignment, "fasta"))
        for seq in protAlnSeqs:
            joinedSeqs.append(seq)
        SeqIO.write(joinedSeqs, tmpFile, "fasta")

        # MUSCLE align TP genes with markers
        runMUSCLE(tmpFile, alnOutput)

        # Gen spHMMs and interval pos
        # Extract spHMM coordinates from MUSCLE alignment
        hmmDict = gensphmmfiles(prot_family_name, alnOutput, tp_genes_prot,
                                hmm_directory, gene_pos_file, gene_pos_file_aa)

        if r1_file_suffix is None:
            r1_file_suffix = ""
        if r2_file_suffix is None:
            r2_file_suffix = ""

        # #Preprocess synthetic reads
        nucl_seq_directory = PreProcessReadsPar(nucl_seq_directory,
                                                seq_fmt,pair_fmt,
                                                r1_file_suffix.strip(),
                                                r2_file_suffix.strip(),
                                                build_op_dir,
                                                CPU_THREADS)

        #Check if BLAST DB directory mapping file is provided or not
        if blast_db_directory_map_file is None:
            blast_db_directory_map_file = ""

        # Translate nucleotide seq
        if not os.path.isdir(prot_seq_directory):
            prot_seq_directory = TranseqReadsDir(build_op_dir, nucl_seq_directory, CPU_THREADS)

        # HMMER Search
        if not os.path.exists(allHMMResult):
            os.makedirs(hmm_search_directory,0o777,True)
            for hmmSeqPosKey, hmmFileObj in hmmDict.items():
                hmmInterval = str(hmmDict[hmmSeqPosKey].intervalStart)+"_"+str(hmmDict[hmmSeqPosKey].intervalEnd)
                RunHMMDirectoryParallel(prot_seq_directory,hmmFileObj.hmmFile, cohort_name, prot_family_name, "30_10", hmmInterval, hmm_search_directory, CPU_THREADS)

            with open(allHMMResult, 'w') as outfile:
                for subdir, dirs, files in os.walk(hmm_search_directory):
                    for file in files:
                        filePath = os.path.join(subdir, file)
                        if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                            with open(filePath) as infile:
                                for line in infile:
                                    outfile.write(line)
        else:
            print("Using existing HMM search result file:" + allHMMResult)

        # BLAST Alignment
        if not os.path.exists(allBLASTResult):
            if not os.path.isdir(blastn_search_directory):
                print("Constructing BLAST Search Dir:" + blastn_search_directory)
                os.makedirs(blastn_search_directory,0o777,True)
                RunPCMakeDBandBlastN(nucl_seq_directory, blast_db_directory_map_file,
                                   tp_genes_nucl, "blastn", "-max_target_seqs 10000 -perc_identity 90.0 -outfmt \"6 sseqid slen sstart send qseqid qlen qstart qend pident evalue\" ",
                                   blastn_search_directory, CPU_THREADS)

            with open(allBLASTResult, 'w') as outfile:
                outfile.write("sseqid\tslen\tsstart\tsend\tqseqid\tqlen\tqstart\tqend\tpident\tevalue\tSample\tsampleType\n")
                for subdir, dirs, files in os.walk(blastn_search_directory):
                    for file in files:
                        filePath = os.path.join(subdir, file)
                        if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                            with open(filePath) as infile:
                                for line in infile:
                                    sampleName = os.path.basename(filePath).split(".txt")[0]
                                    outfile.write(line.strip() + "\t" + sampleName + "\t" + cohort_name + "\n")
        else:
            print("Using existing BLAST search result file:" + allBLASTResult)

        # Eval spHMMs
        rpackages.importr('base')
        utils = rpackages.importr('utils')
        packageNames = ('tidyverse','ggsci','ggpubr','dplyr','ggplot2')
        packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
        if len(packnames_to_install) > 0:
            utils.install_packages(StrVector(packnames_to_install))
        rpackages.importr('tidyverse')
        rpackages.importr('ggsci')
        rpackages.importr('ggpubr')
        rpackages.importr('dplyr')
        rpackages.importr('ggplot2')

        hp_hmm_directory = os.path.join(build_op_dir, 'HiPer_spHMMs')
        os.makedirs(hp_hmm_directory,0o777,True)
        module_dir = os.path.dirname(os.path.abspath(createhmm.__file__))
        print("\nR-script path : " + module_dir)
        r_script = os.path.join(module_dir,'EvaluateSpHMMs.R')

        with open(r_script, 'r') as f:
            rStr = f.read()
        myfunc = STAP(rStr, "EvaluateSpHMM")
        myfunc.EvaluateSpHMM(allHMMResult, allBLASTResult, gene_pos_file, prot_family_name, float(f1_thresh), hmm_directory, hp_hmm_directory)
        timeTaken = time.time() - startTime
        mins = int(timeTaken / 60)
        secs = int(timeTaken) % 60
        print("\nTotal time taken : " + str(mins) + " mins " + str(secs) + " seconds")
        return hp_hmm_directory
    except:
        print("Metabgc-build has failed. Please check your inputs and contact support on : https://github.com/donia-lab/MetaBGC")
        exit()
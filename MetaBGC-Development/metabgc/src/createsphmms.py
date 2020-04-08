"""
Created on Tues Oct 10 03:08:34 2017
Updated on Mon Feb 5 18:40:09 2018
@author: francinecamacho
"""
from __future__ import division
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
import os 
import subprocess
from metabgc.src.hmmrecord import HMMFile

"""
Script takes in an alignment file (fasta format), kmer len, and sliding window len. The script 
parses the alignment file with kmer len and sliding winder to create a new alignment file.
Furthermore, the script builds a HMM profile with parsed portion of the new alignment file.
"""

"""
Function builds HMM profile with the parsed portion of the alignment file. 
"""
def runHMMBuild(alnFile, modelName):
    print("Running HMM Build on:", alnFile)
    hmmFile = alnFile.split('.fas')[0] +".hmm"
    cmd = "hmmbuild -n " + modelName + " --amino "+ hmmFile + " "+ alnFile 
    print(cmd) 
    subprocess.call(cmd, shell=True)
    print("Done Running HMM Build on:",alnFile)
    return hmmFile

"""
Function parses alignment file into kmer parts of the alignment file,
with a sliding window method. 
"""
def getKmers(k, interval, outdir, msaFile, tp_prot_file, modelName, start, end, gene_pos_file,gene_pos_file_aa):
    pprot_TP_dict = {}
    for record in SeqIO.parse(tp_prot_file, "fasta"):
        pprot_TP_dict[record.id] = str(record.seq)

    alignment = AlignIO.read(msaFile, "fasta")
    if (start != None) & (end !=None):
        alignment = alignment[:,start-1:end-1]
    print ("Number of domains: %i" % len(alignment))
    print ("Alignment length: %i" % alignment.get_alignment_length())
    hmmDict = {}
    counter = int(((alignment.get_alignment_length()- k) /interval)+1)
    j = 0

    for i in range(alignment.get_alignment_length()):
        alnCol = alignment[:,i]
        if '-' in alnCol:
            j=j+1
        else:
            break;

    seqCtr = alignment.get_alignment_length()
    for i in range(alignment.get_alignment_length()-1, -1, -1):
        alnCol = alignment[:, i]
        if '-' in alnCol:
            seqCtr=seqCtr-1
        else:
            break;
    gene_pos_out_aa = open(gene_pos_file_aa, 'w')
    gene_pos_out_aa.write("gene_name\tstart\tend\tinterval\tprot_type\n")
    gene_pos_out = open(gene_pos_file, 'w')
    gene_pos_out.write("gene_name\tstart\tend\tinterval\tprot_type\n")
    for i in range(counter):
        startPos = j 
        endPos = j+k
        if endPos <= seqCtr:
            kmer = alignment[:,startPos:endPos] #[ rows (different domains),columns (Amino Acids)]
            spHMMAlign = MultipleSeqAlignment([],Gapped(IUPAC.extended_protein, "-"))

            if str(kmer[0].seq).count("-") > 15:
                # Remove the TP genes
                for align in kmer:
                    if align.id not in pprot_TP_dict:
                        spHMMAlign.append(align)
                    else:
                        prot_seq = str(align.seq)
                        prot_seq = prot_seq.replace("-", "")
                        start_tp_coord = pprot_TP_dict[align.id].find(prot_seq)
                        end_tp_coord = start_tp_coord + len(prot_seq)
                        gene_pos_out_aa.write(align.id + "\t" + str(start_tp_coord) + "\t" + str(end_tp_coord) + "\t" +
                                      str(i*interval) + "_" + str(i*interval+k) +
                                      "\t" + modelName + "\n")
                        gene_pos_out.write(align.id + "\t" + str(start_tp_coord*3) + "\t" + str(end_tp_coord*3) + "\t" +
                                      str(i*interval) + "_" + str(i*interval+k) +
                                      "\t" + modelName + "\n")

                outputFile = outdir + os.sep +modelName + "__" + str(k)+ "_"+ str(interval) + "__"+ str(i*interval)+ "_"+ str(i*interval+k) + ".fas"
                AlignIO.write(spHMMAlign, outputFile, "fasta")
                hmmFile = runHMMBuild(outputFile, modelName)
                hmmSegment = str(startPos)+ "_"+ str(endPos)
                hmmDict[hmmSegment] = HMMFile(i*interval,i*interval+k,hmmFile)

            j = j+interval
        else:
            break
    gene_pos_out_aa.close()
    gene_pos_out.close()
    return hmmDict

def GenerateSpHMM(aln_file, tp_prot_file, window_len, kmer_len, outdir, hmmName, start, end,gene_pos_file,gene_pos_file_aa):
    return getKmers(kmer_len, window_len, outdir, aln_file, tp_prot_file, hmmName, start, end,gene_pos_file,gene_pos_file_aa)


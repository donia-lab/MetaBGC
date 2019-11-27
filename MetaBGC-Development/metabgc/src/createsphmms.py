"""
Created on Tues Oct 10 03:08:34 2017
Updated on Mon Feb 5 18:40:09 2018
@author: francinecamacho
"""
from __future__ import division
from Bio import AlignIO
import os 
import subprocess
import ntpath

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
def getKmers(k, interval, outdir, msaFile, modelName, start, end):

    alignment = AlignIO.read(msaFile, "fasta")
    if (start != None) & (end !=None):
        alignment = alignment[:,start-1:end-1]
    print ("Number of domains: %i" % len(alignment))
    print ("Alignment length: %i" % alignment.get_alignment_length())
    hmmDict = {}
    counter = int(((alignment.get_alignment_length()- k) /interval)+1)
    j = 0
    seqCtr = len(alignment)
    for i in range(alignment.get_alignment_length()):
        alnCol = alignment[:,i]
        if '-' in alnCol:
            j=j+1
        else:
            break;
    for i in range(counter):
        startPos = j 
        endPos = j+k
        if endPos <= alignment.get_alignment_length():
            kmer = alignment[:,startPos:endPos] #[ rows (different domains),columns (Amino Acids)]
            outputFile = outdir + os.sep +modelName + "__" + str(k)+ "_"+ str(interval) + "__"+ str(i*interval)+ "_"+ str(i*interval+k) + ".fas"
            AlignIO.write(kmer, outputFile, "fasta")
            hmmFile = runHMMBuild(outputFile, modelName)
            hmmSegment = str(startPos)+ "_"+ str(endPos)
            hmmDict[hmmSegment] = hmmFile
            j = j+interval
        else:
            break
    return hmmDict

def GenerateSpHMM(aln_file, window_len, kmer_len, outdir, hmmName, start, end ):
    return getKmers(kmer_len, window_len, outdir, aln_file, hmmName, start, end)


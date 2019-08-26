from Bio import SearchIO
import os
import subprocess
from Utils.HMMRecord import HMMRecord
import pandas as pd

"""
Function searches FASTA file against HMM. 
"""
def runHMMSearch(fastaFile, hmmFile, outputFile):
    print('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
    cmd = "hmmsearch --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + outputFile + " " + hmmFile + " "+ fastaFile + " > /dev/null"
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running HMM Build with:",fastaFile)

"""
Function builds HMM profile with the parsed portion of the alignment file. 
"""
def runHMMBuild(alnFile, hmmFile, modelName):
    cmd = "hmmbuild -n " + modelName + " --amino "+ hmmFile + " "+ alnFile
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running HMM Build on:",alnFile)

"""
Function build a BLAST DB with a FASTA. 
"""
def runMakeBLASTDB(fastaFile, database, type):
    cmd = "makeblastdb -in " + fastaFile + " -dbtype nucl -parse_seqids -out " + database
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done building BLAST Build on:",fastaFile)

"""
Function BLAST search a FASTA. 
"""
def runBLASTN(fastaFile, database, outFile):
    cmd = "blastn -query " + fastaFile + " -db database -perc_identity 90.0 -max_target_seqs 10000 -outfmt \"6 sseqid qseqid slen qlen qcovs pident evalue qstart qend\" -out " + outFile
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done running BLAST Build on:",fastaFile)

"""
Function to parse HMM file into HMMRecord dict. 
"""
def parseHMM(hmmPathFile, sampleType, sampleID, cyclaseType, window, interval):
    with open(hmmPathFile, 'rU') as handle:
        results_dict = {}
        try:
            for record in SearchIO.parse(handle, 'hmmer3-text'):
                hits = record.hits
                num_hits = len(hits)  #calculate how many hits per query
                if num_hits > 0:  #extract hits data
                    for i in range(0, num_hits):
                        hmm_name = hits[i].id
                        hmm_bitscore = hits[i].bitscore
                        if hmm_name not in results_dict: # add hits to results dictionary
                            readID_Class = HMMRecord(hmm_name, sampleType, sampleID, cyclaseType, hmm_bitscore, window, interval)
                            results_dict[hmm_name] = readID_Class
                            # the bio module wont allow duplicated query IDs
            handle.close()
        except ValueError:
            print ("ERROR: duplicated queryIDs in hmmer result file:", hmmPathFile)
    return results_dict


def createPandaDF(cyclase_dict, outfile):
    outputDF = pd.DataFrame()
    if len(cyclase_dict) != 0:  # check if counter dict is not empty to add to df
        df = [(k, v.sampleType, v.sampleID, v.cyclaseType, v.bitscore, v.window, v.interval) for k, v in list(cyclase_dict.items())]  # convert dictionary to list
        outputDF = outputDF.append(df)  # append list to our initalized dataframe
        outputDF.columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMMScore", "window","interval"]  # rename column names
        sorted_outputDF = outputDF.sort_values(by=['HMMScore'],ascending=[False])  # sort in decending order by Hit.counts column
        sorted_outputDF.to_csv(outfile, index=False, sep='\t', header=False)
    else:
        outputDF_empty_columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMMScore", "window", "interval"]  # rename column names
        outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
        outputDF_empty.to_csv(outfile, index=False, sep='\t', header=False)
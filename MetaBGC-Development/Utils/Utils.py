from Bio import SearchIO
from Bio import SeqIO
import os
import subprocess
from Utils.HMMRecord import HMMRecord
import pandas as pd
import re
import itertools

"""
Function searches all FASTA file in a directory against a HMM. 
"""
def RunHMMDirectory(inputDir, hmmModel, sampleType, protType, window, interval, ouputDir, ncpus=4):
    for subdir, dirs, files in os.walk(inputDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.fasta$", file) and os.path.getsize(filePath) > 0:
                sampleStr = file.split(".")[0]
                hmmTblFileName = sampleStr +"_"+interval+".tbl"
                hmmTblFilePath = os.path.join(ouputDir, hmmTblFileName)
                runHMMSearch(filePath, hmmModel,hmmTblFilePath,ncpus)
                result_dict = parseHMM(hmmTblFilePath, sampleType, sampleStr, protType, window, interval)
                hmmSearchFileName = sampleStr +"_"+interval+".txt"
                hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
                createPandaDF(result_dict, hmmSearchFilePath)

"""
Function searches FASTA file against HMM. 
"""
def runHMMSearch(fastaFile, hmmFile, outputFile, ncpus=4):
    print('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
    cmd = "hmmsearch --cpu " + str(ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + outputFile + " " + hmmFile + " "+ fastaFile + " > /dev/null"
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
def runMakeBLASTDB(fastaFile, dbName, dbPath, type):
    dbOut = dbPath + os.sep + dbName
    cmd = "makeblastdb -in " + fastaFile + " -title " + dbName +" -dbtype nucl -parse_seqids -out " + dbOut
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done building BLAST Build on:",fastaFile)

"""
Function BLAST search a FASTA. 
"""
def runBLASTN(fastaFile, database, outFile, ncpus=4):
    cmd = "blastn -num_threads " + str(ncpus) +  " -query " + fastaFile + " -db " + database + " -perc_identity 90.0 -max_target_seqs 10000 -outfmt \"6 sseqid slen sstart send qseqid qlen qstart qend pident evalue\" -out " + outFile
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done running BLAST Build on:",fastaFile)

"""
Function to run BLAST against a directory. 
"""
def RunBLASTNDirectory(dbDir, queryFile, ouputDir,ncpus=4):
    for subdir, dirs, files in os.walk(dbDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*\.fasta$", file) and os.path.getsize(filePath) > 0:
                runMakeBLASTDB(filePath, "TMPDB", ouputDir, 'nucl')
                sampleStr = file.split(".")[0]
                outputFileName = sampleStr + ".txt"
                outputFilePath = os.path.join(ouputDir, outputFileName)
                runBLASTN(queryFile, ouputDir+os.sep+"TMPDB", outputFilePath, ncpus)

"""
Function to parse HMM file into HMMRecord dict. 
"""
def parseHMM(hmmPathFile, sampleType, sampleID, cyclaseType, window, interval):
    with open(hmmPathFile, 'rU') as handle:
        results_dict = {}
        try:
            for record in SearchIO.parse(handle, 'hmmer3-tab'):
                hits = record.hits
                num_hits = len(hits)  #calculate how many hits per query
                if num_hits > 0:  #extract hits data
                    for i in range(0, num_hits):
                        hmm_name = hits[i].id
                        hmm_bitscore = hits[i].bitscore
                        if hmm_name not in results_dict: # add hits to results dictionary
                            hmmRec = HMMRecord(hmm_name, sampleType, sampleID, cyclaseType, hmm_bitscore, window, interval)
                            results_dict[hmm_name] = hmmRec
                            print(hmmRec)
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


"""
Function runs MSA against FASTA file. 
"""
def runMUSCLE(fastaFile,outputFile):
    print('Running MUSCLE with {0}.'.format(fastaFile))
    cmd = "muscle -in " + fastaFile + " -out " + outputFile
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running MUSCLE with:",fastaFile)

"""
Function to generate a 6 frame translated amino acid sequence from the nucleotide sequence. 
"""
def runTranSeq(fastaFile,frameCode,outputFile):
    print('Running transeq with {0}.'.format(fastaFile))
    cmd = "transeq " + fastaFile + " " + outputFile + " -frame="+ frameCode +" -table=0 -sformat pearson"
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running transeq with:",fastaFile)


"""
Function runs cd-hit on a FASTA file. 
"""
def runCDHit(fastaFile,outputFile,ncpu):
    print('Running CD-Hit with {0}.'.format(fastaFile))
    opPrefix = outputFile.split(".fasta")[0]
    cmd = "cd-hit-est -i " + fastaFile + " -o " + opPrefix + " -c .95 -n 10 -d 0 -aS .95 -T "+ str(ncpu)
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running Cd-Hit with:",fastaFile)


"""
Interleave FR reads. 
"""
def interleave(iter1, iter2) :
    for (forward, reverse) in zip(iter1,iter2):
        yield forward
        yield reverse

"""

"""
def PreProcessReads(nucl_seq_directory,seq_fmt,pair_fmt,R1_file_suffix,R2_file_suffix,ouputDir):
    if seq_fmt.lower() == "fasta" and (pair_fmt.lower() == "interleaved" or pair_fmt.lower() == "single"):
        return nucl_seq_directory

    out_seq_directory = os.path.join(ouputDir, 'nucl_seq_dir')
    if os.path.isdir(out_seq_directory):
        return out_seq_directory
    elif seq_fmt.lower() == "fasta" and pair_fmt.lower() == "split":
        os.makedirs(out_seq_directory, 0o777, True)
        for subdir, dirs, files in os.walk(nucl_seq_directory):
            for file in files:
                filePath = os.path.join(subdir, file)
                regF = r".*" + R1_file_suffix + "$"
                regR = r".*" + R2_file_suffix + "$"
                if re.match(regF, file) and not re.match(regR, file) and os.path.getsize(filePath) > 0:
                    sampleName = file.split(R1_file_suffix)[0]
                    r2FilePath = os.path.join(subdir, sampleName + R2_file_suffix)
                    file_out = os.path.join(out_seq_directory, sampleName + ".fasta")
                    records_f = SeqIO.parse(open(filePath, "rU"), seq_fmt.lower())
                    records_r = SeqIO.parse(open(r2FilePath, "rU"), seq_fmt.lower())
                    handle = open(file_out, "w")
                    count = SeqIO.write(interleave(records_f, records_r), handle, "fasta")
                    handle.close()
        return out_seq_directory
    elif seq_fmt.lower() == "fastq":
        os.makedirs(out_seq_directory, 0o777, True)
        if pair_fmt.lower() == "interleaved" or pair_fmt.lower() == "single":
            for subdir, dirs, files in os.walk(nucl_seq_directory):
                for file in files:
                    filePath = os.path.join(subdir, file)
                    if re.match(r".*.fastq$", file) and os.path.getsize(filePath) > 0:
                        print("Pre-processing:", file)
                        out_seq_path = os.path.join(out_seq_directory, file)
                        count = SeqIO.convert(filePath, "fastq", out_seq_path, "fasta")
        elif pair_fmt.lower() == "split":
            for subdir, dirs, files in os.walk(nucl_seq_directory):
                for file in files:
                    filePath = os.path.join(subdir, file)
                    regF = r".*"+R1_file_suffix+"$"
                    regR = r".*"+R2_file_suffix+"$"
                    if re.match(regF, file) and not re.match(regR, file) and os.path.getsize(filePath) > 0:
                        print("Pre-processing:", file)
                        sampleName = file.split(R1_file_suffix)[0]
                        r2FilePath = os.path.join(subdir,sampleName + R2_file_suffix)
                        file_out = os.path.join(out_seq_directory,sampleName + ".fasta")
                        records_f = SeqIO.parse(open(filePath, "rU"), seq_fmt.lower())
                        records_r = SeqIO.parse(open(r2FilePath, "rU"), seq_fmt.lower())
                        handle = open(file_out, "w")
                        count = SeqIO.write(interleave(records_f, records_r), handle, "fasta")
                        handle.close()
        else:
            print("Invalid sequence pair format inputs.\n")
            exit(0)
        return out_seq_directory
    else:
        print("Invalid file format inputs.\n")
        exit(0)


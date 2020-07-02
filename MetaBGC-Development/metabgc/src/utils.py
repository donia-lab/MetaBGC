from Bio import SearchIO
from Bio import SeqIO
import os
import subprocess
import metabgc.src.hmmrecord as hmmrecord
import pandas as pd
import re
from multiprocessing import Pool, freeze_support
from itertools import repeat
import shutil
import csv
import logging

"""
Function searches all FASTA file in a directory against a HMM. 
"""
def RunHMMDirectory(inputDir, hmmModel, sampleType, protType, window, interval, ouputDir, ncpus=16):
    for subdir, dirs, files in os.walk(inputDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.fasta$", file) and os.path.getsize(filePath) > 0:
                sampleStr = os.path.splitext(file)[0]
                runHMMSearch(filePath, hmmModel, ouputDir, sampleType, sampleStr, protType, window, interval,ncpus)

"""
Function searches all FASTA file in a directory against a HMM in parallel. 
"""
def RunHMMDirectoryParallel(inputDir, hmmModel, sampleType, protType, window, interval, ouputDir, ncpus=4):
    logging.info('Number of pool processes:{0}.'.format(ncpus))
    for subdir, dirs, files in os.walk(inputDir):
        fastaFileList=[]
        sampleStrList = []
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.fasta$", file) and os.path.getsize(filePath) > 0:
                sampleStr = os.path.splitext(file)[0]
                fastaFileList.append(filePath)
                sampleStrList.append(sampleStr)
                if len(fastaFileList) >= ncpus:
                    HMMSearchParallel(fastaFileList, hmmModel, ouputDir, sampleType, sampleStrList, protType, window, interval)
                    fastaFileList = []
                    sampleStrList = []
        if len(fastaFileList) > 0:
            HMMSearchParallel(fastaFileList, hmmModel, ouputDir, sampleType, sampleStrList, protType, window, interval)

"""
Function to run make HMM search against FASTA files in parallel. 
"""
def HMMSearchParallel(fastaFileList, hmmFile, ouputDir, sampleType, sampleStrList, protType, window, interval):
    numOfprocess = len(fastaFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(runHMMSearch, zip(fastaFileList, repeat(hmmFile), repeat(ouputDir), repeat(sampleType),
                                   sampleStrList, repeat(protType),repeat(window),repeat(interval),repeat(1)))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Function searches FASTA file against HMM. 
"""
def runHMMSearch(fastaFile, hmmFile, ouputDir, sampleType, sampleStr, protType, window, interval, ncpus=16):
    hmmTblFileName = sampleStr + "__" + interval + ".tbl"
    hmmTblFilePath = os.path.join(ouputDir, hmmTblFileName)
    if not os.path.exists(hmmTblFilePath):
        logging.info('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
        cmd = "hmmsearch --cpu " + str(ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + hmmTblFilePath + " " + hmmFile + " "+ fastaFile + " > /dev/null"
        logging.info(cmd)
        subprocess.call(cmd, shell=True)
    else:
        logging.info("HMM search skipped... Using existing result for: " + fastaFile)
    result_dict = parseHMM(hmmTblFilePath, "hmmer3-tab", sampleType, sampleStr, protType, window, interval)
    hmmSearchFileName = sampleStr + "__" + interval + ".txt"
    hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
    createPandaDF(result_dict, hmmSearchFilePath)
    logging.info("Done Running HMM Build with:" + fastaFile)

"""
Function builds HMM profile with the parsed portion of the alignment file. 
"""
def runHMMBuild(alnFile, hmmFile, modelName):
    cmd = "hmmbuild -n " + modelName + " --amino "+ hmmFile + " "+ alnFile
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done Running HMM Build on:" + alnFile)

"""
Function build a BLAST DB with a FASTA. 
"""
def runMakeBLASTDB(fastaFile, dbName, dbOpPath, type):
    dbOut = dbOpPath + os.sep + dbName
    cmd = "makeblastdb -in " + fastaFile + " -title " + dbName +" -dbtype nucl -out " + dbOut
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done building BLAST build on:" + fastaFile)
    return dbOut

"""
Function BLAST search a FASTA. 
"""
def runBLASTN(fastaFile, database, blastCommand, blastParamStr, outFile, ncpus=4):
    cmd = blastCommand + " -num_threads " + str(ncpus) + " -query " + fastaFile + " -db " + database + " " + blastParamStr + " -out " + outFile
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done running BLAST search on:" + fastaFile)

"""
Function to run make BLAST db and search a FASTA file. 
"""
def MakeSearchBLASTN(dbInputFile, existingDbDirMapFile, dbOpPath, searchFile, blastCmdString, blastParamStr, outFile):
    sample_basename = os.path.basename(dbInputFile)
    dbOut = ""
    if existingDbDirMapFile:
        logging.info("Existing database path map provided:" + existingDbDirMapFile)
        map_dict = {}
        map_file = os.path.join(existingDbDirMapFile)
        with open(map_file) as f:
            reader = csv.reader(f, skipinitialspace=True)
            map_dict = dict(reader)
        dbOut = map_dict[sample_basename]

    if not os.path.isfile(dbOut):
        logging.info("Constructing BLAST DB for:" + dbInputFile)
        dbOpPath = dbOpPath + os.sep + sample_basename
        os.makedirs(dbOpPath, 0o777, True)
        dbName = os.path.splitext(sample_basename)[0]
        dbOut = runMakeBLASTDB(dbInputFile, dbName, dbOpPath, 'nucl')
        runBLASTN(searchFile, dbOut, blastCmdString, blastParamStr, outFile, 1)
        shutil.rmtree(dbOpPath)
    else:
        logging.info("Found existing database path:" + dbOut)
        runBLASTN(searchFile, dbOut, blastCmdString, blastParamStr, outFile, 1)
"""
Function to run make BLAST db and search a FASTA file. 
"""
def MakeSearchBLASTNParallel(dbFileList, existingDbDir, dbOpPath,
                             searchFileList, blastCmdString, blastParamStr,
                             outFileList):
    numOfprocess = len(dbFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(MakeSearchBLASTN, zip(dbFileList, repeat(existingDbDir),repeat(dbOpPath),
                                       searchFileList, repeat(blastCmdString), repeat(blastParamStr),
                                       outFileList))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Function to create blast databases from the fasta files in dbDir and running a query against the directory. 
"""
def RunMakeDBandBlastN(dbDir, existingDbDir, queryFile, blastCmdString, blastParamStr, ouputDir, ncpus=4):
    for subdir, dirs, files in os.walk(dbDir):
        dbFileList=[]
        searchFileList = []
        outFileList = []
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*\.fasta$", file) and os.path.getsize(filePath) > 0:
                sampleStr = os.path.splitext(file)[0]
                outputFileName = sampleStr + ".txt"
                outputFilePath = os.path.join(ouputDir, outputFileName)
                dbFileList.append(filePath)
                searchFileList.append(queryFile)
                outFileList.append(outputFilePath)
                if len(dbFileList)>=ncpus:
                    MakeSearchBLASTNParallel(dbFileList,existingDbDir,
                                             ouputDir,searchFileList,
                                             blastCmdString,blastParamStr,
                                             outFileList)
                    dbFileList = []
                    searchFileList = []
                    outFileList = []
        if len(dbFileList) > 0:
            MakeSearchBLASTNParallel(dbFileList,existingDbDir,
                                     ouputDir, searchFileList, blastCmdString,
                                     blastParamStr, outFileList)



"""
Function to run BLAST search against a FASTA file. 
"""
def BLASTParallel(blastdb, searchFileList, blastCmdString, blastParamStr,outFileList):
    numOfprocess = len(searchFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(runBLASTN, zip(searchFileList, repeat(blastdb),
                                repeat(blastCmdString), repeat(blastParamStr),
                                outFileList,repeat(1)))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Function to search a bunch of blast queries against the database provided. 
"""
def RunBlastSearch(blastdb, queryFileList, blastCmdString, blastParamStr, ouputDir, ncpus=4):
    searchFileList = []
    outFileList = []
    for file in queryFileList:
        if os.path.exists(file):
            searchFileList.append(file)
            rootFileName = os.path.splitext(os.path.basename(file))[0]
            outFileName = os.path.join(ouputDir,rootFileName + ".txt")
            outFileList.append(outFileName)
        if len(searchFileList) >= ncpus:
            BLASTParallel(blastdb, searchFileList,
                            blastCmdString, blastParamStr,outFileList)
            searchFileList = []
            outFileList = []
    if len(searchFileList) > 0:
        BLASTParallel(blastdb, searchFileList,
                      blastCmdString, blastParamStr, outFileList)

"""
Function to parse HMM file into HMMRecord dict. 
"""
def parseHMM(hmmPathFile, hmm_string_fmt,sampleType, sampleID, protType, window, interval):
    with open(hmmPathFile, 'r') as handle:
        results_dict = {}
        try:
            for record in SearchIO.parse(handle, hmm_string_fmt):
                hits = record.hits
                num_hits = len(hits)  #calculate how many hits per query
                if num_hits > 0:  #extract hits data
                    for i in range(0, num_hits):
                        hmm_name = hits[i].id
                        hmm_bitscore = hits[i].bitscore
                        if hmm_name not in results_dict: # add hits to results dictionary
                            hmmRec = hmmrecord.HMMRecord(hmm_name, sampleType, sampleID, protType, hmm_bitscore, window, interval)
                            results_dict[hmm_name] = hmmRec
            handle.close()
        except ValueError:
            print ("ERROR: duplicated queryIDs in hmmer result file:", hmmPathFile)
    return results_dict


def createPandaDF(hmm_dict, outfile):
    outputDF = pd.DataFrame()
    if len(hmm_dict) != 0:  # check if counter dict is not empty to add to df
        df = [(k, v.sampleType, v.sampleID, v.protType, v.bitscore, v.window, v.interval) for k, v in list(hmm_dict.items())]  # convert dictionary to list
        outputDF = outputDF.append(df)  # append list to our initalized dataframe
        outputDF.columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window","interval"]  # rename column names
        sorted_outputDF = outputDF.sort_values(by=['HMMScore'],ascending=[False])  # sort in decending order by Hit.counts column
        sorted_outputDF.to_csv(outfile, index=False, sep='\t', header=False)
    else:
        outputDF_empty_columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window", "interval"]  # rename column names
        outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
        outputDF_empty.to_csv(outfile, index=False, sep='\t', header=False)


"""
Function runs MSA against FASTA file. 
"""
def runMUSCLE(fastaFile,outputFile):
    logging.info('Running MUSCLE with {0}.'.format(fastaFile))
    cmd = "muscle -in " + fastaFile + " -out " + outputFile
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done Running MUSCLE with:" + fastaFile)

"""
Function to generate a 6 frame translated amino acid sequence from the nucleotide sequence. 
"""
def runTranSeq(fastaFile,frameCode,outputFile):
    logging.info('Running transeq with {0}.'.format(fastaFile))
    cmd = "transeq " + fastaFile + " " + outputFile + " -frame="+ frameCode +" -table=11 -sformat pearson"
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done Running transeq with:" + fastaFile)


"""
Function runs cd-hit on a FASTA file. 
"""
def runCDHit(fastaFile,outputFile,ncpu):
    logging.info('Running CD-Hit with {0}.'.format(fastaFile))
    cmd = "cd-hit-est -i " + fastaFile + " -o " + outputFile + " -c .95 -n 10 -d 0 -aS .95 -M 4098 -T "+ str(ncpu)
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    logging.info("Done Running Cd-Hit with:" + fastaFile)

"""
Interleave FR reads. 
"""
def interleave(iter1, iter2) :
    for (forward, reverse) in zip(iter1,iter2):
        yield forward
        yield reverse

def InterleaveReads(out_seq_directory, sampleName, r1FilePath,r2FilePath,seq_fmt):
    records_f = SeqIO.parse(open(r1FilePath, "rU"), seq_fmt)
    records_r = SeqIO.parse(open(r2FilePath, "rU"), seq_fmt)
    file_out = os.path.join(out_seq_directory, sampleName + ".fasta")
    handle = open(file_out, "w")
    count = SeqIO.write(interleave(records_f, records_r), handle, "fasta")
    handle.close()

def InterleaveReadsParallel(out_seq_directory,sampleNameList,r1FilePathList,r2FilePathList,seq_fmt):
    numOfprocess = len(sampleNameList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(InterleaveReads, zip(repeat(out_seq_directory),sampleNameList, r1FilePathList, r2FilePathList,
                                      repeat(seq_fmt.lower())))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

def TranseqReadsParallel(nSeqFileNameList,frameCode,protFileNameList):
    numOfprocess = len(nSeqFileNameList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(runTranSeq, zip(nSeqFileNameList, repeat(frameCode),protFileNameList))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Convert all files in a directory to corresponding 6 frame protein sequence
"""
def TranseqReadsDir(build_op_dir,nucl_seq_directory,ncpus):
    prot_seq_directory = os.path.join(build_op_dir, 'prot_seq_dir')
    os.makedirs(prot_seq_directory, 0o777, True)
    for subdir, dirs, files in os.walk(nucl_seq_directory):
        nSeqFileNameList = []
        protFileNameList = []
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*\.fasta$", file) and os.path.getsize(filePath) > 0:
                prot_file = prot_seq_directory + os.sep + os.path.basename(filePath)
                nSeqFileNameList.append(filePath)
                protFileNameList.append(prot_file)
            if len(nSeqFileNameList) >= ncpus:
                TranseqReadsParallel(nSeqFileNameList, "6", protFileNameList)
                nSeqFileNameList = []
                protFileNameList = []
        if len(nSeqFileNameList) > 0:
            TranseqReadsParallel(nSeqFileNameList, "6", protFileNameList)
    return prot_seq_directory

"""
Preprocess reads. Convert fastq files into interleaved paired-end FASTA files for BAST database creation...
"""

def PreProcessReadsPar(nucl_seq_directory, seq_fmt, pair_fmt, R1_file_suffix, R2_file_suffix, ouputDir,ncpu):
    if seq_fmt.lower() == "fasta" and (pair_fmt.lower() == "interleaved" or pair_fmt.lower() == "single"):
        return nucl_seq_directory

    out_seq_directory = os.path.join(ouputDir, 'nucl_seq_dir')
    if os.path.isdir(out_seq_directory):
        return out_seq_directory
    elif seq_fmt.lower() == "fasta" and pair_fmt.lower() == "split":
        os.makedirs(out_seq_directory, 0o777, True)
        for subdir, dirs, files in os.walk(nucl_seq_directory):
            sampleNameList = []
            r1FilePathList=[]
            r2FilePathList = []
            freeze_support()
            for file in files:
                filePath = os.path.join(subdir, file)
                regF = r".*" + R1_file_suffix + "$"
                regR = r".*" + R2_file_suffix + "$"
                if re.match(regF, file) and not re.match(regR, file) and os.path.getsize(filePath) > 0:
                    logging.info("Pre-processing:" + file)
                    sampleName = file.split(R1_file_suffix)[0]
                    r2FilePath = os.path.join(subdir, sampleName + R2_file_suffix)
                    sampleNameList.append(sampleName)
                    r1FilePathList.append(filePath)
                    r2FilePathList.append(r2FilePath)
                if len(sampleNameList)>=ncpu:
                    InterleaveReadsParallel(out_seq_directory,sampleNameList,r1FilePathList,r2FilePathList,seq_fmt.lower())
                    sampleNameList = []
                    r1FilePathList = []
                    r2FilePathList = []
            if len(sampleNameList) > 0:
                InterleaveReadsParallel(out_seq_directory, sampleNameList, r1FilePathList, r2FilePathList,
                                        seq_fmt.lower())
        return out_seq_directory
    elif seq_fmt.lower() == "fastq":
        os.makedirs(out_seq_directory, 0o777, True)
        if pair_fmt.lower() == "interleaved" or pair_fmt.lower() == "single":
            for subdir, dirs, files in os.walk(nucl_seq_directory):
                for file in files:
                    filePath = os.path.join(subdir, file)
                    if re.match(r".*.fastq$", file) and os.path.getsize(filePath) > 0:
                        logging.info("Pre-processing:" + file)
                        out_seq_path = os.path.join(out_seq_directory, file)
                        count = SeqIO.convert(filePath, "fastq", out_seq_path, "fasta")
        elif pair_fmt.lower() == "split":
            for subdir, dirs, files in os.walk(nucl_seq_directory):
                sampleNameList = []
                r1FilePathList = []
                r2FilePathList = []
                freeze_support()
                for file in files:
                    filePath = os.path.join(subdir, file)
                    regF = r".*" + R1_file_suffix + "$"
                    regR = r".*" + R2_file_suffix + "$"
                    if re.match(regF, file) and not re.match(regR, file) and os.path.getsize(filePath) > 0:
                        logging.info("Pre-processing:" + file)
                        sampleName = file.split(R1_file_suffix)[0]
                        r2FilePath = os.path.join(subdir, sampleName + R2_file_suffix)
                        sampleNameList.append(sampleName)
                        r1FilePathList.append(filePath)
                        r2FilePathList.append(r2FilePath)
                    if len(sampleNameList) >= ncpu:
                        InterleaveReadsParallel(out_seq_directory, sampleNameList, r1FilePathList, r2FilePathList,
                                                seq_fmt.lower())
                        sampleNameList = []
                        r1FilePathList = []
                        r2FilePathList = []
                if len(sampleNameList) > 0:
                    InterleaveReadsParallel(out_seq_directory, sampleNameList, r1FilePathList, r2FilePathList,
                                            seq_fmt.lower())
        else:
            logging.info("Invalid sequence pair format inputs.\n")
            exit(0)
        return out_seq_directory
    else:
        logging.info("Invalid file format inputs.\n")
        exit(0)


"""
Preprocess reads. Convert fastq files into interleaved paired-end FASTA files for BAST database creation...
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
                        logging.info("Pre-processing:" + file)
                        out_seq_path = os.path.join(out_seq_directory, file)
                        count = SeqIO.convert(filePath, "fastq", out_seq_path, "fasta")
        elif pair_fmt.lower() == "split":
            for subdir, dirs, files in os.walk(nucl_seq_directory):
                for file in files:
                    filePath = os.path.join(subdir, file)
                    regF = r".*"+R1_file_suffix+"$"
                    regR = r".*"+R2_file_suffix+"$"
                    if re.match(regF, file) and not re.match(regR, file) and os.path.getsize(filePath) > 0:
                        logging.info("Pre-processing:" + file)
                        sampleName = file.split(R1_file_suffix)[0]
                        r2FilePath = os.path.join(subdir,sampleName + R2_file_suffix)
                        file_out = os.path.join(out_seq_directory,sampleName + ".fasta")
                        records_f = SeqIO.parse(open(filePath, "rU"), seq_fmt.lower())
                        records_r = SeqIO.parse(open(r2FilePath, "rU"), seq_fmt.lower())
                        handle = open(file_out, "w")
                        count = SeqIO.write(interleave(records_f, records_r), handle, "fasta")
                        handle.close()
        else:
            logging.info("Invalid sequence pair format inputs.\n")
            exit(0)
        return out_seq_directory
    else:
        logging.info("Invalid file format inputs.\n")
        exit(0)

"""
Function searches all FASTA file in a directory against a HMM in parallel. 
"""
def RunHMMDirectoryParallelReduced(inputDir, hmmModel, ouputDir, ncpus=4):
    logging.info('Number of pool processes:{0}.'.format(ncpus))
    for subdir, dirs, files in os.walk(inputDir):
        fastaFileList=[]
        sampleStrList = []
        for file in files:
            filePath = os.path.join(subdir, file)
            if os.path.getsize(filePath) > 0:
                sampleStr = os.path.splitext(file)[0]
                fastaFileList.append(filePath)
                sampleStrList.append(sampleStr)
                if len(fastaFileList) >= ncpus:
                    HMMSearchParallelReduced(fastaFileList, hmmModel, ouputDir)
                    fastaFileList = []
                    sampleStrList = []
        if len(fastaFileList) > 0:
            HMMSearchParallelReduced(fastaFileList, hmmModel, ouputDir)

"""
Function to run make HMM search against FASTA files in parallel. 
"""
def HMMSearchParallelReduced(fastaFileList,hmmFile,ouputDir):
    numOfprocess = len(fastaFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(runHMMSearchReduced, zip(fastaFileList, repeat(hmmFile), repeat(ouputDir), repeat(1)))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Function searches FASTA file against HMM. 
"""
def runHMMSearchReduced(fastaFile, hmmFile, ouputDir, ncpus=4):
    hmmTblFileName = os.path.splitext(os.path.basename(fastaFile))[0] + ".tbl"
    hmmTblFilePath = os.path.join(ouputDir, hmmTblFileName)
    if not os.path.exists(hmmTblFilePath):
        logging.info('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
        cmd = "hmmsearch --cpu " + str(ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + hmmTblFilePath + " " + hmmFile + " "+ fastaFile + " > /dev/null"
        logging.info(cmd)
        subprocess.call(cmd, shell=True)
    else:
        logging.info("HMM search skipped... Using existing result for: " + fastaFile)
    result_dict = parseHMM(hmmTblFilePath, "hmmer3-tab", "", "", "", "", "")
    hmmSearchFileName = os.path.splitext(os.path.basename(fastaFile))[0] + ".txt"
    hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
    createPandaDF(result_dict, hmmSearchFilePath)
    logging.info("Done running HMM Build with:" + fastaFile)
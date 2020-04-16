from Bio import SearchIO
from Bio import SeqIO
import os
import subprocess
from metabgc.src.hmmrecord import HMMRecord
import pandas as pd
import re
from multiprocessing import Pool, freeze_support
from itertools import repeat
import shutil
import ntpath


"""
Constructs the expected BLAST DB path for DoniaLab; Specific to runs done in Donia environment.
"""
def ConstructDoniaDBPath(baseDir,sampleFileName):
    sampleType = sampleFileName.split('_')[2]
    sampleName = os.path.splitext(sampleFileName)[0]
    dbOut = os.path.join(baseDir,sampleType,sampleName,sampleName+"-raw-reads-fasta",sampleFileName)
    return dbOut

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
    print('Number of pool processes:{0}.'.format(ncpus))
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
def HMMSearchParallel(fastaFileList, hmmFile, ouputDir, sampleType,sampleStrList,protType,window,interval):
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
        print('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
        cmd = "hmmsearch --cpu " + str(ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + hmmTblFilePath + " " + hmmFile + " "+ fastaFile + " > /dev/null"
        print(cmd)
        subprocess.call(cmd, shell=True)
    else:
        print("HMM search skipped... Using existing result for: ", fastaFile)
    result_dict = parseHMM(hmmTblFilePath, "hmmer3-tab", sampleType, sampleStr, protType, window, interval)
    hmmSearchFileName = sampleStr + "__" + interval + ".txt"
    hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
    createPandaDF(result_dict, hmmSearchFilePath)
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
def runMakeBLASTDB(fastaFile, dbName, dbOpPath, type):
    dbOut = dbOpPath + os.sep + dbName
    cmd = "makeblastdb -in " + fastaFile + " -title " + dbName +" -dbtype nucl -out " + dbOut
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done building BLAST Build on:",fastaFile)
    return dbOut

"""
Function BLAST search a FASTA. 
"""
def runBLASTN(fastaFile, database, blastParamStr, outFile, ncpus=16):
    cmd = "blastn -num_threads " + str(ncpus) + " -query " + fastaFile + " -db " + database + " " + blastParamStr + " -outfmt \"6 sseqid slen sstart send qseqid qlen qstart qend pident evalue\" -out " + outFile
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done running BLAST Build on:",fastaFile)

"""
Function to run make BLAST db and search a FASTA file. 
"""
def MakeSearchBLASTN(dbFile, existingDbDir, dbOpPath, searchFile, blastParamStr, outFile):
    basename = os.path.basename(dbFile)
    if existingDbDir:
        print("Existing database path provided:" + existingDbDir)
        dbOut = ConstructDoniaDBPath(existingDbDir,basename)

    if not os.path.isfile(dbOut):
        print("Constructing BLAST DB for:" + dbFile)
        dbOpPath = dbOpPath + os.sep + basename
        os.makedirs(dbOpPath, 0o777, True)
        dbName = os.path.splitext(basename)[0]
        dbOut = runMakeBLASTDB(dbFile,dbName,dbOpPath,'nucl')
        runBLASTN(searchFile, dbOut, blastParamStr, outFile, 1)
        shutil.rmtree(dbOpPath)
    else:
        print("Found existing database path:" + dbOut)
        runBLASTN(searchFile, dbOut, blastParamStr, outFile, 1)

"""
Function to run make BLAST db and search a FASTA file. 
"""
def MakeSearchBLASTNParallel(dbFileList, existingDbDir,dbOpPath, searchFileList, blastParamStr, outFileList):
    numOfprocess = len(dbFileList)
    pool = Pool(processes=numOfprocess)
    pool.starmap(MakeSearchBLASTN, zip(dbFileList, repeat(existingDbDir),repeat(dbOpPath), searchFileList, repeat(blastParamStr), outFileList))
    pool.close()
    pool.join()
    pool.terminate()  # garbage collector

"""
Function to run BLAST against a directory. 
"""
def RunBLASTNDirectoryPar(dbDir, existingDbDir,queryFile, blastParamStr, ouputDir,ncpus=4):
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
                    MakeSearchBLASTNParallel(dbFileList,existingDbDir,ouputDir,searchFileList, blastParamStr,outFileList)
                    dbFileList = []
                    searchFileList = []
                    outFileList = []
        if len(dbFileList) > 0:
            MakeSearchBLASTNParallel(dbFileList,existingDbDir,ouputDir, searchFileList , blastParamStr, outFileList)

"""
Function to run BLAST against a directory. 
"""
def RunBLASTNDirectory(dbDir, queryFile, blastParamStr, ouputDir,ncpus=4):
    for subdir, dirs, files in os.walk(dbDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*\.fasta$", file) and os.path.getsize(filePath) > 0:
                runMakeBLASTDB(filePath, "TMPDB", ouputDir, 'nucl')
                sampleStr = os.path.splitext(file)[0]
                outputFileName = sampleStr + ".txt"
                outputFilePath = os.path.join(ouputDir, outputFileName)
                runBLASTN(queryFile, ouputDir+os.sep+"TMPDB", blastParamStr,outputFilePath, ncpus)

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
                            hmmRec = HMMRecord(hmm_name, sampleType, sampleID, protType, hmm_bitscore, window, interval)
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
    cmd = "transeq " + fastaFile + " " + outputFile + " -frame="+ frameCode +" -table=11 -sformat pearson"
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Done Running transeq with:",fastaFile)


"""
Function runs cd-hit on a FASTA file. 
"""
def runCDHit(fastaFile,outputFile,ncpu):
    print('Running CD-Hit with {0}.'.format(fastaFile))
    cmd = "cd-hit-est -i " + fastaFile + " -o " + outputFile + " -c .95 -n 10 -d 0 -aS .95 -T "+ str(ncpu)
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
                prot_file = prot_seq_directory + os.sep + ntpath.basename(filePath)
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
                    print("Pre-processing:", file)
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
                        print("Pre-processing:", file)
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
                        print("Pre-processing:", file)
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
            print("Invalid sequence pair format inputs.\n")
            exit(0)
        return out_seq_directory
    else:
        print("Invalid file format inputs.\n")
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

"""
Function searches all FASTA file in a directory against a HMM in parallel. 
"""
def RunHMMDirectoryParallelReduced(inputDir, hmmModel, ouputDir, ncpus=4):
    print('Number of pool processes:{0}.'.format(ncpus))
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
        print('Running HMM Search with {0} against {1}.'.format(fastaFile, hmmFile))
        cmd = "hmmsearch --cpu " + str(ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + hmmTblFilePath + " " + hmmFile + " "+ fastaFile + " > /dev/null"
        print(cmd)
        subprocess.call(cmd, shell=True)
    else:
        print("HMM search skipped... Using existing result for: ", fastaFile)
    result_dict = parseHMM(hmmTblFilePath, "hmmer3-tab", "", "", "", "", "")
    hmmSearchFileName = os.path.splitext(os.path.basename(fastaFile))[0] + ".txt"
    hmmSearchFilePath = os.path.join(ouputDir, hmmSearchFileName)
    createPandaDF(result_dict, hmmSearchFilePath)
    print("Done Running HMM Build with:",fastaFile)
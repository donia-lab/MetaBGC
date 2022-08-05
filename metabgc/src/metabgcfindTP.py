from metabgc.src.utils import *
import pandas as pd
from pandas import DataFrame
import os
from metabgc.src.extractfastaseq import RunExtractDirectoryPar
from metabgc.src.extractfastaseq import RunExtractDescription

def mbgcfindtp(alnFile, prot_seq_directory, nucl_seq_directory, out_dir, do_alignment):
    # setup paths
    identifyReadIds = os.path.join(out_dir,'CombinedReadIds.txt')
    allHMMResult = os.path.join(out_dir,'CombinedHMMResults.txt')
    fasta_seq_dir = os.path.join(out_dir,'extract_seq')
    aa_multiFastaFile = os.path.join(out_dir,'CombinedHMMHits.faa')
    nc_multiFastaFile = os.path.join(out_dir,'CombinedHMMHits.fna')

    alnOutput=alnFile.split('.fasta')[0] +".faa"
    if do_alignment:
        runMUSCLE(alnFile, alnOutput)
    else:
        alnOutput = alnFile

    hmmFile=alnFile.split('.fasta')[0] +".hmm"
    modelName=alnFile.split('.fasta')[0]
    runHMMBuild(alnOutput, hmmFile, modelName)
    hmm_search_directory = os.path.join(out_dir, "hmm_search")
    os.makedirs(hmm_search_directory, 0o777, True)
    RunHMMDirectoryParallelReduced(prot_seq_directory, hmmFile, hmm_search_directory, 4)

    with open(allHMMResult, 'w') as outfile:
        for subdir, dirs, files in os.walk(hmm_search_directory):
            for file in files:
                filePath = os.path.join(subdir, file)
                sample = os.path.splitext(file)[0]
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(sample + '\t' +line)

    df_HMMMatch = pd.read_csv(allHMMResult, delimiter="\t",
                              names=["Sample", "readID", "sampleType", "Sample2", "protType", "HMMScore", "window",
                                     "interval"])
    df_HMMMatch = df_HMMMatch.drop(columns=["sampleType", "Sample2", "protType", "window", "interval"])
    df_HMMMatch.to_csv(identifyReadIds, index=False, sep='\t')

    os.makedirs(fasta_seq_dir, 0o777, True)
    RunExtractDirectoryPar(prot_seq_directory, identifyReadIds, fasta_seq_dir,
                           aa_multiFastaFile, "faa", True, True, 1)
    record_desc_dict = RunExtractDescription(aa_multiFastaFile,"fasta")

    RunExtractDirectoryPar(nucl_seq_directory, identifyReadIds, fasta_seq_dir,
                           nc_multiFastaFile, "fna", False, False, 1)

    df_record_desc = DataFrame(list(record_desc_dict.items()), columns=['readID', 'Description'])
    df_HMMMatch = pd.merge(df_HMMMatch,df_record_desc,on=['readID'],how='inner')
    df_HMMMatch.to_csv(allHMMResult, index=False, sep='\t')







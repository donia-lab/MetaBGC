from metabgc.src.utils import *
import pandas as pd
import os
import sys
from metabgc.src.extractfastaseq import RunExtractDirectoryPar

if __name__ == '__main__':
    # setup paths
    alnFile = sys.argv[1]
    hmm_search_directory = sys.argv[2]
    prot_seq_directory = sys.argv[3]
    out_dir = sys.argv[4]
    identifyReadIds = os.path.join(out_dir,'CombinedReadIds.txt')
    allHMMResult = os.path.join(out_dir,'CombinedHMMResults.txt')
    fasta_seq_dir = os.path.join(out_dir,'extract_seq')
    multiFastaFile = os.path.join(out_dir,'CombinedHMMHits.faa')

    alnOutput=alnFile.split('.fasta')[0] +".faa"
    runMUSCLE(alnFile, alnOutput)
    hmmFile=alnFile.split('.fasta')[0] +".hmm"
    modelName=alnFile.split('.fasta')[0]
    runHMMBuild(alnOutput, hmmFile, modelName)
    os.makedirs(hmm_search_directory, 0o777, True)
    RunHMMDirectoryParallel(prot_seq_directory, hmmFile, hmm_search_directory, 4)

    with open(allHMMResult, 'w') as outfile:
        for subdir, dirs, files in os.walk(hmm_search_directory):
            for file in files:
                filePath = os.path.join(subdir, file)
                sample = os.path.splitext(file)[0]
                if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(sample + '\t' +line)

    df_HMMMatch = pd.read_csv(allHMMResult, delimiter="\t", names=["Sample", "readID", "sampleType", "Sample2", "protType", "HMMScore", "window","interval"])
    df_HMMMatch = df_HMMMatch.drop(columns=["sampleType", "Sample2", "protType", "window","interval"])
    df_HMMMatch.to_csv(identifyReadIds, index=False, sep='\t')
    os.makedirs(fasta_seq_dir, 0o777, True)
    RunExtractDirectoryPar(prot_seq_directory, identifyReadIds, fasta_seq_dir, multiFastaFile, "faa", 1)



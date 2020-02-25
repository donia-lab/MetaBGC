from metabgc.src.utils import *

if __name__ == '__main__':
    alnFile='/tigress/DONIA/abiswas/MetaBGCRuns/TP_Finder/Rum-like-rSAM/Rum-like-rSAMs-homologs.fasta'
    hmmFile=alnFile.split('.fasta')[0] +".hmm"
    modelName=alnFile.split('.fasta')[0]
    hmm_search_directory='/tigress/DONIA/abiswas/MetaBGCRuns/TP_Finder/Rum-like-rSAM/hmm_search'
    prot_seq_directory='/tigress/DONIA/abiswas/MetaBGCRuns/TP_Finder/prot_seq'
    runHMMBuild(alnFile, hmmFile, modelName)
    os.makedirs(hmm_search_directory, 0o777, True)
    RunHMMDirectoryParallel(prot_seq_directory, hmmFile, hmm_search_directory, 20)
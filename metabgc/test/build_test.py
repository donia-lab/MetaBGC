import metabgc.src.metabgcbuild as build
import metabgc.src.evaluate_sphmms as evaluate
from Bio import SeqIO


def test_ungappedseqsearch():
    assert build.ungappedseqsearch("--MSE-HDTDV---LVGGSM", "TDV-LV") == [8, 16]
    assert build.ungappedseqsearch("HLAGEPLETCLRYGAIAGAYACTIPATRAGAIDRAALLRPAA--------",
                                   "ACTIPATRAGAIDRAALLRPAA--------") == [20, 42]

def test_GenerateSpHMM():
    prot_family_name = "AbcK"
    prot_aln_file = "AbcK/data/TP_Homolog_Alignment.afa"
    tp_prot_file = "AbcK/data/TPGenes.faa"
    gene_pos_file = "AbcK/output/build/Gene_Interval_Pos.txt"
    gene_pos_file_aa = "AbcK/output/build/Gene_Interval_Pos_AA.txt"
    hmm_directory = "AbcK/output/build/spHMMs"
    build.gensphmmfiles(prot_family_name, prot_aln_file, tp_prot_file, hmm_directory, gene_pos_file, gene_pos_file_aa)


def test_evaluate_sphmm_1():
    HMMRunFile = "AbcK/data/CombinedHmmSearch.txt"
    BLAST_TP_NoCov_File = "AbcK/data/CombinedBLASTSearch.txt"
    GeneIntervalPosFile = "AbcK/output/build/Gene_Interval_Pos.txt"
    HMM_Model_Name = 'AbcK'
    F1_Threshold = 0.50
    HMMOutDir = "AbcK/output/build/spHMMs"
    HMMHighPerfOutDir = "AbcK/output/build/HiPerf"
    evaluate.evaluate_sphmms(HMMRunFile,
                             BLAST_TP_NoCov_File,
                             GeneIntervalPosFile,
                             HMM_Model_Name,
                             F1_Threshold,
                             HMMOutDir,
                             HMMHighPerfOutDir)

def test_evaluate_sphmm_2():
    HMMRunFile = "RdhA/build/CombinedHmmSearch.txt"
    BLAST_TP_NoCov_File = "RdhA/build/CombinedBLASTSearch.txt"
    GeneIntervalPosFile = "RdhA/build/Gene_Interval_Pos.txt"
    HMM_Model_Name = 'RdhA'
    F1_Threshold = 0.50
    HMMOutDir = "RdhA/build/spHMMs"
    HMMHighPerfOutDir = "RdhA/build/HiPer_spHMMs"
    evaluate.evaluate_sphmms(HMMRunFile,
                             BLAST_TP_NoCov_File,
                             GeneIntervalPosFile,
                             HMM_Model_Name,
                             F1_Threshold,
                             HMMOutDir,
                             HMMHighPerfOutDir)


import pytest
import metabgc.src.metabgcbuild as build
import metabgc.src.evaluate_sphmms as evaluate
from Bio import SeqIO
import os

def test_ungappedseqsearch():
    assert build.ungappedseqsearch("--MSE-HDTDV---LVGGSM","TDV-LV") == [8,16]
    assert build.ungappedseqsearch("HLAGEPLETCLRYGAIAGAYACTIPATRAGAIDRAALLRPAA--------","ACTIPATRAGAIDRAALLRPAA--------") == [20,42]

def test_gengeneposlist():
    prot_aln_file ="AbcK/data/AcbK-homologs-for-HMM.fasta"
    hmm_directory ="AbcK/output/build/spHMMs"
    prot_family_name ="AbcK"
    alnOutput ="AbcK/data/AcbK-homologs-plus-three-synthetic-TPs.fasta"
    gene_pos_file ="AbcK/output/build/Gene_Interval_Pos.txt"
    module_dir = os.path.dirname(pytest.__file__)
    hmmDict = build.gensphmmfiles("AbcK", prot_aln_file, hmm_directory)
    protAlnSeqs = list(SeqIO.parse(prot_aln_file, "fasta"))
    build.gengeneposlist(prot_family_name, protAlnSeqs, hmmDict, alnOutput, gene_pos_file)

def test_GenerateSpHMM():
    prot_family_name="AbcK"
    prot_aln_file="AbcK/data/TP_Homolog_Alignment.afa"
    hmm_directory="AbcK/output/build/spHMMs"
    tp_prot_file="AbcK/data/TPGenes.faa"
    gene_pos_file="AbcK/output/build/Gene_Interval_Pos.txt"
    gene_pos_file_aa = "AbcK/output/build/Gene_Interval_Pos_AA.txt"
    build.gensphmmfiles(prot_family_name, prot_aln_file, tp_prot_file, hmm_directory, gene_pos_file, gene_pos_file_aa)

def test_evaluate_sphmm():
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
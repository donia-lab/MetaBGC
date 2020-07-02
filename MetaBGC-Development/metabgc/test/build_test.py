import pytest
import metabgc.src.metabgcbuild as build
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
    hmm_directory="AbcK/output/build/HiPerf"
    tp_prot_file="AbcK/data/TPGenes.faa"
    gene_pos_file="AbcK/output/build/Gene_Interval_Pos.txt"
    gene_pos_file_aa = "AbcK/output/build/Gene_Interval_Pos_AA.txt"
    build.gensphmmfiles(prot_family_name, prot_aln_file, tp_prot_file, hmm_directory, gene_pos_file, gene_pos_file_aa)
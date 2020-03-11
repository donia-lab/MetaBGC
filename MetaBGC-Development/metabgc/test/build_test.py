import pytest
from metabgc.src.metabgcbuild import ungappedseqsearch
from metabgc.src.metabgcbuild import gensphmmfiles
from metabgc.src.metabgcbuild import gengeneposlist
from Bio import AlignIO
from Bio import SeqIO
import os

def test_ungappedseqsearch():
    assert ungappedseqsearch("--MSE-HDTDV---LVGGSM","TDV-LV") == [8,16]
    assert ungappedseqsearch("HLAGEPLETCLRYGAIAGAYACTIPATRAGAIDRAALLRPAA--------","ACTIPATRAGAIDRAALLRPAA--------") == [20,42]

def test_gengeneposlist():
    prot_aln_file ="AbcK/data/AcbK-homologs-for-HMM.fasta"
    hmm_directory ="AbcK/output/build/spHMMs"
    prot_family_name ="AbcK"
    alnOutput ="AbcK/data/AcbK-homologs-plus-three-synthetic-TPs.fasta"
    gene_pos_file ="AbcK/output/build/Gene_Interval_Pos.txt"
    module_dir = os.path.dirname(pytest.__file__)
    hmmDict = gensphmmfiles("AbcK", prot_aln_file, hmm_directory)
    protAlnSeqs = list(SeqIO.parse(prot_aln_file, "fasta"))
    gengeneposlist(prot_family_name, protAlnSeqs, hmmDict, alnOutput, gene_pos_file)
# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

MetaBGC is a read-based pipeline for the detection of biosynthetic gene clusters (BGCs) directly in human microbiome-derived metagenomic sequencing data, without the need for bacterial cultivation or sequencing. 

## Getting Started

These instructions will get you run MetaBGC on your local machine.

### Prerequisites

To run MetaBGC, you will need the following dependencies:

* [Python](https://www.python.org/downloads/) (version >= 3.6)
* [Biopython](https://biopython.org/wiki/Download) v1.72
* [scikit-learn](https://scikit-learn.org/stable/install.html) v0.20.1
* [Pandas](https://pandas.pydata.org/pandas-docs/stable/install.html) v0.19.2 
* [rpy2](https://pypi.org/project/rpy2/) v2.9.1
* [EMBOSS Transeq](http://emboss.sourceforge.net/download/) version 6.6.0.0
* [HMMER](http://hmmer.org/download.html) version 3.1b2
* [CD-HIT-EST](https://github.com/weizhongli/cdhit/releases) version 4.7
* [ncbi-blast-2.7.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
* [R 3.6.1](http://lib.stat.cmu.edu/R/CRAN/)
* [R Studio] (https://www.rstudio.com/products/rstudio/download/)

### Program Structure

MetaBGC conists of three main modules: 

**MetaBGC-Build** - This module builds, evaluates, and selects high performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. If the BGC of interest commonly contains proteins of the one of the pre-built high performance spHMMs, this step can be skipped.    

**MetaBGC-Identify** - This modules includes running segmented profile Hidden Markov Models (spHMMs) on translated metagenomic reads and parsing the results to identify reads that are biosynthetic.

**MetaBGC-Quantify** - This module includes de-replicating the results from MetaBGC-Identify module and quantifying.

**MetaBGC-Cluster** - This module generates bins of reads from MetaBGC-Quantify module.

### Running MetaBGC-Build to build spHMMs

1. Construct a YAML file with required input files, output files, and input parameters. Some YAML files are included for Lantibiotics, Siderophore and 4 cyclase that were used for benchmarking. These can be run with the data files in [benchmark data folder](https://drive.google.com/drive/folders/1MV88Qv4qMmYAtU159uFFwIczZYzPIRCK?usp=sharing). The config parameters are:
```
default:
    DataRoot: "/home/metabgc/benchmark_data/siderophores"	# Root data directory. All the input paths and output paths are relative to this directory.(REQUIRED)
    InputFiles.ProtAlnFasta: "siderophore.fasta"	# Alignment of the protein homologs in FASTA format. (REQUIRED)
    InputFiles.HMMRun: "Combined_spHMM_Search.txt"	# 
    InputFiles.BLAST_FP: "blastx-fp-subfive.txt"	# 
    InputFiles.BLAST_TP_NoCov: "TP-reads/All_TP_Blast.txt"
    InputFiles.GeneIntervalPos: "TP-reads/Gene_Interval_Pos.txt"
    InputFiles.GeneIntervalBlast: "TP-reads/Gene_Interval_Blast.txt"
    InputFiles.HMMUniqueBlast: "hmm-unique-analysis/HMM_Unique_Blast.txt"
    InputParam.HMM_Model_Name: "siderophore"	# HMM model name.  
    InputParam.Seq_start: 1	# Start position to build models.
    InputParam.Seq_end: 489	# End position to build models.
    InputParam.F1_Threshold: 0.50	# F1 cutoff threshold.
    OutputFiles.HMMOutDir: "output/spHMM" 
    OutputFiles.HMMHighPerfOutDir: "output/spHMM/HighPerf"
    OutputFiles.SampleReadIdDir: "output/Sample_ReadIds"
    OutputFiles.F1_Cutoff: "output/siderophores_F1_Cutoff.txt"
    OutputFiles.Plot: "output/siderophores_synthetic_F1-final.eps"
    OutputFiles.ModelData: "output/siderophores-final_models-data.txt"
    OutputFiles.RawModelData: "output/siderophores-final_data-raw.txt"
``` 
2. Replace the config file in the R notebook MetaBGC-Build/SPHMM_Model.Rmd
3. Run the R notebook. The output high performance spHMMs will be in OutputFiles.HMMHighPerfOutDir  

### Running MetaBGC-Identify to detect biosynthetic-like reads

1. Translate metagenomic reads from nucleic acid to amino acid using EMBOSS:6.6.0.0 transeq tool in all six open reading frames: `-frame=6 - table=0 -sformat pearson` 
2. Run `hmmsearch` command from HMMER on the translated metagenomic data from step 1 using the heuristic filters: `--F1 0.02 --F2 0.02 --F3 0.02` 
3. Parse results from step 2 using the parser script. 
4. Filter the results with the cutoffs for each spHMM. 
5. Extract read ids per sample from filtered results and extract sequences for read ids to create a fasta file for each sample. 

### Running MetaBGC-Quantify to profile biosynthetic-like reads

1. Combine fasta files from MetaBGC-Identify module and de-replicate at 95% identity and 95% alignment coverage using CD-HIT- EST with the following parameters: `-c .95 -n 10 -d 0 -aS .95`
2. Quantify de-replicated fasta file using BLASTn against all sample metagenomes using the following parameters: `-task blastn -dust no -max_target_seqs 1000000 -perc_identity 95.0 -qcov_hsp_perc 50 window_size 11`
3. Produce an abundance profile using the results from step2 BLASTn. 

### Running MetaBGC-Cluster to generate BGCs
1. Abundance profiles from the MetaBGC-Quantify module are clustered using DBSCAN clustering method with Pearson correlation distance metric and the following parameters: `eps .2 min_samples 1`
2. Each bin from DBSCAN must contain at least 50 reads or more to count as a cluster to examine.

## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](LICENSE.md) file for details.


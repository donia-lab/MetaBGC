# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

MetaBGC is a read-based pipeline for the detection of biosynthetic gene clusters (BGCs) directly in human microbiome-derived metagenomic sequencing data, without the need for bacterial cultivation or sequencing. 

## Getting Started

These instructions will get you run MetaBGC on your local machine.

### Prerequisites

To run MetaBGC, you will need the following dependencies and data preprocessing:

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
* [R Studio](https://www.rstudio.com/products/rstudio/download/)
* Before using MetaBGC to either identify Type II Polyketides or other protein families of interest in metagenomic dataset(s), a user must prepare their dataset(s) by translating the metagenomic reads from nucleic acid to amino acid using **EMBOSS:6.6.0.0 transeq** tool in all six open reading frames with the following parameters: `-frame=6 - table=0 -sformat pearson`. 



### Program Structure

MetaBGC consists of three main modules:

**MetaBGC-Build** - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. If the BGC of interest commonly contains proteins of one of the pre-built high-performance spHMMs, this step can be skipped.

**MetaBGC-Identify** - This module reformats HMMER results, filters the data using spHMMs cutoffs and parses out the read IDs for each sample that passes during the filtering stage.

**MetaBGC-Quantify** - This module takes reads identified and that scored above a defined cutoff by the MetaBGC-Identify module as "biosynthetic." A user must then combine these biosynthetic reads into a multi-FASTA file to run de-replication and BLAST which quantifies these de-replicated reads in all samples of the entire metagenomic dataset(s). An abundance matrix is generated for all unique reads against all samples. 

**MetaBGC-Cluster** - In this module, the abundance matrix file produced from MetaBGC-Quantify is the input for MetaBGC-Cluster. This module uses Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to cluster reads with similar coverage profiles across different metagenomic samples into distinct bins. 

### Running MetaBGC-Build to build spHMMs

1. Construct a YAML file with required input files, output files, and input parameters. Some YAML files are included for Lantibiotics, Siderophore and 4 cyclases that were used for benchmarking. These can be run with the data files in the [benchmark data folder](https://github.com/donia-lab/MetaBGC-TIIPKS/tree/master/benchmark_data). A sample template EmptyTemplate.yaml is provided in the MetaBGC-Build directory for adding new dataset(s). The config parameters are:
	1.  DataRoot - REQUIRED: Root data directory. All the input paths and output paths are relative to this directory.
	2.  InputFiles.ProtAlnFasta - REQUIRED: Alignment of protein homologs.  
	3.  InputFiles.HMMRun - REQUIRED: HMMER search of the synthetic reads against all the spHMM models. 
	4.  InputFiles.BLAST_TP_NoCov - REQUIRED: TP reads matched to the proteins with BLAST. 
	5.  InputFiles.BLAST_IGNORE - OPTIONAL: Reads ignored in BLAST that was missed by antiSMASH.
	6.  InputFiles.GeneIntervalPos - REQUIRED: antiSMASH output with the interval positions of each gene. 
	7.  InputFiles.GeneIntervalBlast - OPTIONAL: BLAST result of aligning the reads to the antiSMASH intervals in each protein. 
	8.  InputFiles.HMMUniqueBlast - OPTIONAL: Filtered blast search output that uniquely matches the domain.
	9.  InputFiles.HMM_Cutoff_Scores - REQUIRED: File with HMM cutoffs to compare to BLAST interval reads.
	10. InputFiles.ScaffoldCheck - OPTIONAL: List of scaffold ids for the median, +5, and -5 that are known to match the scaffolds that have the protein of interest. Reduces running time. 
	11. InputParam.HMM_Model_Name - REQUIRED: Name of the PFAM. 
	12. InputParam.Seq_start - REQUIRED: Start position in the PFAM alignment to start building spHMMs. 
	13. InputParam.Seq_end - REQUIRED: End position in the PFAM alignment to start building spHMMs. 
	14. InputParam.F1_Threshold - REQUIRED: F1 score threshold for eliminating low-performance models and tuning spHMMs. 
	15. OutputFiles.HMMOutDir - REQUIRED: Output directory where the spHMMs are saved. 
	16. OutputFiles.HMMHighPerfOutDir - REQUIRED: Output directory where the high performing spHMMs are saved. 
	17. OutputFiles.SampleReadIdDir - REQUIRED: Output directory where read ids matched to spHMMs are saved. 
	18. OutputFiles.F1_Cutoff - REQUIRED: Output file with a list of spHMMs that pass the F1 threshold.  
	19. OutputFiles.Plot - REQUIRED: Output plot of the F1 scores for each spHMM with threshold boundary lines.
	20. OutputFiles.ModelData - REQUIRED: Output FP/TP for HMM reads detected in final models with final cutoffs that keep reads duplicated with the highest HMM Score. In case of a tie, then choose the first one. 
	21. OutputFiles.RawModelData - REQUIRED: Output FP/TP for HMM reads detected in final models with final cutoffs duplicate reads in intervals. 
 
2. Replace the config file in the R notebook MetaBGC-Build/SPHMM_Model.Rmd
3. Run the R notebook. The output high-performance spHMMs will be in OutputFiles.HMMHighPerfOutDir  

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


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
* Before using MetaBGC to either identify Type II Polyketides (Type II PKS) or other protein families of interest in metagenomic dataset(s), a user must prepare their dataset(s) by translating the metagenomic reads from nucleic acid to amino acid using **EMBOSS:6.6.0.0 transeq** tool in all six open reading frames with the following parameters: `-frame=6 - table=0 -sformat pearson`. 



### Program Structure

MetaBGC consists of three main modules:

**MetaBGC-Build** - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. If the BGC of interest commonly contains proteins of one of the pre-built high-performance spHMMs, this step can be skipped.

**MetaBGC-Identify** - This module reformats HMMER results, filters the data using spHMMs HMMER Score cutoffs and parses out the read IDs for each sample that passes during the filtering stage.

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

1. For identifying Type II PKS BGCs you can find the pre-built high-performing [spHMM models](https://github.com/donia-lab/MetaBGC-TIIPKS/tree/master/models) and run HMMER on your translated metagenomic dataset(s) using the `hmmsearch` command with the heuristic filters: `--F1 0.02 --F2 0.02 --F3 0.02`. Please create an HMMER result parent directory with each spHMM results as a sub-directory for the next step.  

2. To parse and reformat HMMER results use the [spHMM parser](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Identify/spHMM_parser.py) with the following inputs:

	```
	A. --hmmscan_file_dir, required=True: The HMMER results parent directory 
	B. --outdir, required=True: The desired parsed output directory name
	C. --cyclase_type,required=True: Protein family name (ie, TcmN or siderophore)
	D. --window, required=False, default='30_10': Window used for segmentation
	E. --interval, required=True: The name of the spHMM interval parsing (ie, 0_30)
	F. --sampleID, required=True: The name of the Sample
	G. --sampleType, required=True: The type of sample (ie,bodysite, isolation_source, cohort)
	```
	
3. To filter the HMMER results from *step 2* you need to combine the parsed results and create a new combined HMMER result file and then use [MetaBGC-Identify](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Identify/MetaBGC-Identify.py). For Type II Polyketides you can use the pre-defined HMMER Score cutoffs [here](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/models/spHMM-cutoffs.txt).

	```
	A. --hmm_file, required=True: Combined spHMM results file  
	B. --outdir, required=True: PATH to save filtered results
	C. --cutoff_file, required=True: Text file with HMM Score cutoffs per interval
	D. --fasta_dir, required=True: PATH to save FASTA files for each sample and their reads that passed the defined HMMER score cutoffs.
	```
	
 **If your protein of interest is not Type II PKS, then use the cutoffs defined from MetaBGC-build. Please be aware because these cutoffs were designated using synthetic genomes, you may need to manually re-tune the spHMM cutoffs after running MetaBGC on a metagenomic dataset(s).**

4. MetaBGC-Identify will produce a text file of each sample and their read IDs that passed in *step 3*. Use the [extractFASTAfromLIST-2MultiFASTA.py](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Identify/extractFASTAfromLIST-2MultiFASTA.py) script to parse out the nucleotide sequences for those sample's reads that passed. 

### Running MetaBGC-Quantify to profile biosynthetic-like reads

1. A user must combine each sample FASTA file from the MetaBGC-Identify module into one multi-FASTA which then can be used to de-replicate at 95% identity and 95% alignment coverage using **CD-HIT-EST** with the following parameters: `-c .95 -n 10 -d 0 -aS .95`.

	**If a user is interested in adding metadata to the headers of the multi-FASTA file such as sampletype, cyclasetype you can use this [script](https://github.com/donia-lab/MetaBGC-TIIPKS/tree/master/MetaBGC-Quantify/add_metadata_to_FASTA.py).**

2. To quantify these de-replicated reads, a users must use the multi-FASTA de-replicated file as query using BLASTn against all sample metagenomes using the following parameters: `-task blastn -dust no -max_target_seqs 1000000 -perc_identity 95.0 -qcov_hsp_perc 50 window_size 11`
**Please name BLAST tabular files with the following synthax 
{SAMPLE_\__ID}_\__{BLAST_\__extension}. For example: V1.UC48.0_\__reads_\__against_\__combined-reads**
3. To produce an abundance profile using the results from *step 2* a user must then use [MetaBGC-Quantify.py](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Quantify/MetaBGC-Quantify.py) with the following parameters: 

	```
		A. --blast_extension, required=True: extension to combine BLAST results to search for all the BLAST tabular files in a directory  
		B. --sample_extension, required=True: extension to parse sample name out of BLAST tabular file
		C. --blast_dirpath, required=True: PATH to BLAST results
		D. --tabular_colnames, nargs='+', required=False, default = "sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue"
		E. --outfile, required=False, default='combined-blast_quantifier-results.txt'
		F.--outdir, required=True: PATH to deposit output results
		G. --cohort_name', required=True: Dataset name 
	```

### Running MetaBGC-Cluster to generate BGCs bins 
1. To generate BGC bins, users will use the abundance profile file **(abundance_table.txt)** produced from the MetaBGC-Quantify module as input for [MetaBGC-Cluster.py](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Cluster/MetaBGC-Cluster.py) with the following suggested parameters: `eps .2 min_samples 1`

2. We suggest for synthetic dataset to examine bins that contain at least 50 reads or more and at least 10 reads or more for metagenomic data. 

**Please be aware these are suggested parameters and may have to be tuned depending on the metagenomic dataset(s).**


## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
If there are any questions please contact:

**Mohamed Abou Donia, Ph.D.**
 <br /> Assistant Professor
 <br /> Department of Molecular Biology
 <br /> Princeton University
 <br /> Email: <donia@princeton.edu>




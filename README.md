# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

MetaBGC is a read-based pipeline for the detection of biosynthetic gene clusters (BGCs) directly in human microbiome-derived metagenomic sequencing data, without the need for bacterial cultivation or sequencing. 

## Getting Started

These instructions will get you setup to MetaBGC on your local Linux or Apple computer. 

### Prerequisites

To run MetaBGC, you will need the following dependencies and data preprocessing:

* [Python](https://www.python.org/downloads/) (version >= 3.6)
* [Biopython](https://biopython.org/wiki/Download) v1.72
* [scikit-learn](https://scikit-learn.org/stable/install.html) v0.20.1
* [Pandas](https://pandas.pydata.org/pandas-docs/stable/install.html) v0.19.2 
* [rpy2](https://pypi.org/project/rpy2/) v2.9.1
* [fuzzywuzzy](https://github.com/seatgeek/fuzzywuzzy) 
* [python-Levenshtein](https://pypi.org/project/python-Levenshtein)
* [MUSCLE 3.8.31](https://www.drive5.com/muscle/downloads.htm)
* [EMBOSS Transeq](http://emboss.sourceforge.net/download/) version 6.6.0.0
* [HMMER](http://hmmer.org/download.html) version 3.1b2
* [CD-HIT-EST](https://github.com/weizhongli/cdhit/releases) version 4.7
* [ncbi-blast-2.7.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
* [R 3.6.1](http://lib.stat.cmu.edu/R/CRAN/)
* [EMBOSS-6.5.7](ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz)


### Program Structure

MetaBGC consists of four main modules:

**MetaBGC-Build** - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. If the BGC of interest commonly contains proteins of one of the pre-built high-performance spHMMs, this step can be skipped.

**MetaBGC-Identify** - This module reformats HMMER results, filters the data using spHMMs HMMER Score cutoffs and parses out the read IDs for each sample that passes during the filtering stage.

**MetaBGC-Quantify** - This module takes reads identified and that scored above a defined cutoff by the MetaBGC-Identify module as "biosynthetic." A user must then combine these biosynthetic reads into a multi-FASTA file to run de-replication and BLAST which quantifies these de-replicated reads in all samples of the entire metagenomic dataset(s). An abundance matrix is generated for all unique reads against all samples. 

**MetaBGC-Cluster** - In this module, the abundance matrix file produced from MetaBGC-Quantify is the input for MetaBGC-Cluster. This module uses Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to cluster reads with similar coverage profiles across different metagenomic samples into distinct bins. 

### Running MetaBGC-Build to build spHMMs

1. To build and evaluate spHMMs for the protein family of choice being used to target the BGC of interest, the MetaBGC-Build.py script has to be executed with required input files. To select high performance spHMMs, a synthetic dataset of reads must be generated with reads from true positive genes spiked in to test the performance of each spHMM.

	```
    A. --prot_alignment, required=True: Alignment of the protein homologs in FASTA format.
	B. --prot_family_name, required=True: Name of the protein family. 
         This is used as prefix for spHMM files.
	C. --cohort_name, required=True: Name of the cohort of read samples used for evaluation.
	D. --nucl_seq_directory, required=True:  Directory of reads of metagenome samples. The filenames are used as sample names.  
    E. --seq_fmt, required=True: {fasta, fastq} Sequence file format and extension.
    F. --pair_fmt, required=True: {single, split, interleaved} Sequence pair format.
    G. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    H. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	I. --tp_genes_nucl, required=True: Nucleotide sequence in multi-FASTA format with true positive genes in the synthetic dataset used to generate the synthetic reads.
    J. --F1_Thresh, required=False: Threshold of the F1 score for selection of high performing spHMMs (Def.=0.5).
	K. --output_directory, required=True: Directory to save results. 
	L. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The high-performing spHMMs will be saved in the HiPer_spHMMs folder in the output directory specified. The HiPer_spHMMs folder should have the following files:
	```
    A. plot.eps : F1 score plot of all the spHMMs and the F1 cutoff threshold. 
    B. *.hmm : A set of spHMMs that perform above the F1 cutoff threshold.
    C. F1_Cutoff.txt: HMM search cutoff scores to be used for each high-performing spHMM interval. 
    	**If your protein of interest is not Type II PKS, then use the cutoffs defined from MetaBGC-Build. Please be aware because these cutoffs were designated using synthetic genomes, you may need to manually re-tune the spHMM cutoffs after running MetaBGC on a metagenomic dataset(s).**
    ```
    
### Running MetaBGC-Identify to detect biosynthetic-like reads

1. For identifying BGCs using the spHMMs constructed in the build step the MetaBGC-Identify.py script should be executed with the required input files. To identify Type II PKS BGCs you can find the pre-built high-performing [spHMM models](https://github.com/donia-lab/MetaBGC-TIIPKS/tree/master/models). 
	```
	A. --sphmm_directory, required=True: The high performing HMMER model directory generated from MetaBGC-Build.
	B. --nucl_seq_directory, required=True: Directory of reads of metagenome samples. The filenames are used as sample names.  
    C. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    D. --pair_fmt, required=True: {single, split, interleaved} Sequence pair format.
    G. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    H. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	G. --prot_family_name, required=True: Name of the protein family.
	H. --cohort_name, required=True: Name of the cohort of read samples used for model building.
	I. --output_directory, required=True: Directory to save results.
	J. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```

2. MetaBGC-Identify will produce a FASTA file, **CombinedIDFASTASeqs.fasta**, of the nucleotide sequences for those sample's reads that passed the cutoff in the output directory. 

### Running MetaBGC-Quantify to profile biosynthetic-like reads

1. To quantify these reads identified in the previous step, the MetaBGC-Quantify.py script should be executed with the following parameters: 

	```
    A. --identify_fasta, required=True: Path to the CombinedIDFASTASeqs.fasta file produced by MetaBGC-Identify.
    B. --nucl_seq_directory, required=True: Directory of reads of metagenome samples. The filenames are used as sample names.  
    C. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    D. --pair_fmt, required=True: {single, split, interleaved} Sequence pair format.
    G. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    H. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	G. --cohort_name, required=True: Name of the cohort of read samples used for model building.
	H. --output_directory, required=True: Directory to save results.
	I. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The output of the quantify script is an abundance profile file **abundance_table.txt**. 

### Running MetaBGC-Cluster to generate BGCs bins 
1. To generate BGC bins, users will use the abundance profile file, **abundance_table.txt**, produced from the MetaBGC-Quantify module as input for [MetaBGC-Cluster.py](https://github.com/donia-lab/MetaBGC-TIIPKS/blob/master/MetaBGC-Cluster/MetaBGC-Cluster.py) with the following suggested parameters: `eps .2 min_samples 1`

2. We suggest for synthetic dataset to examine bins that contain at least 50 reads or more and at least 10 reads or more for metagenomic data. 

	**Please be aware these are suggested parameters and may have to be tuned depending on the metagenomic dataset(s).**


## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](LICENSE.md) file for details.

Please feel free to report bugs at our [github issue page](https://github.com/donia-lab/MetaBGC-TIIPKS/wiki).

## Contact
If there are any questions please contact:

**Mohamed Abou Donia, Ph.D.**
 <br /> Assistant Professor
 <br /> Department of Molecular Biology
 <br /> Princeton University
 <br /> Email: <donia@princeton.edu>




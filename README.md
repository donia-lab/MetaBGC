# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

MetaBGC is a read-based algorithm for the detection of biosynthetic gene clusters (BGCs) directly in metagenomic sequencing data.

## Getting Started

These instructions will get you setup to MetaBGC on your local Linux or Apple environment. 

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
* [EMBOSS-6.5.7](http://emboss.sourceforge.net/)

### Installation 

Add the MetaBGC-V1_Scripts folder to PYTHONPATH.

```
export PYTHONPATH="<PATH>/MetaBGC/MetaBGC-V1_Scripts"
```

The python scripts for each module in available in the folder with the same name and can be run directly from commandline. 


### Program Structure

MetaBGC consists of four main modules:

**MetaBGC-Build** - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. Pre-built high-performance spHMMs exist for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs). If any of these protein families is to be used, this step can be skipped.

**MetaBGC-Identify** - This module runs on translated metagenomic reads from a cohort of samples using a selected set of high-performance spHMMs and their pre-set score cutoffs, as determined in MetaBGC-Build. The results are parsed into a list of identified biosynthetic reads in fasta format.

**MetaBGC-Quantify** - This module de-replicates all biosynthetic reads discovered by MetaBGC-Identify from all metagenomic samples in the cohort into a unified set of “unique biosynthetic reads”. An abundance profile martrix is then generated for all “unique biosynthetic reads” by quantifying them yingin all samples of the metagenomic cohort.

**MetaBGC-Cluster** - This module uses Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to cluster “unique biosynthetic reads” with similar abundance profiles across different metagenomic samples into distinct bins, based on the abundance profile martrix obtained in MetaBGC-Quantify. 

### Running MetaBGC-Build to build spHMMs

1. To build and evaluate spHMMs for the protein family of interest, the MetaBGC-Build.py script has to be executed with required input files. To select high performance spHMMs, a synthetic metagenomic dataset must be generated with reads from true positive genes spiked in to test the performance of each spHMM.

	```
    A. --prot_alignment, required=True: Alignment of homologs from the protein family of interest in FASTA format.
	B. --prot_family_name, required=True: Name of the protein family. 
         This is used as prefix for spHMM files.
	C. --cohort_name, required=True: Name of the cohort of synthetic metagenomic samples used for evaluation.
	D. --nucl_seq_directory, required=True:  Directory of reads for the synthetic metagenomic samples. The filenames are used as sample names.  
    E. --seq_fmt, required=True: {fasta, fastq} Sequence file format and extension.
    F. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    G. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    H. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	I. --tp_genes_nucl, required=True: Nucleotide sequence of the full-length true-positive genes in the synthetic dataset in multi-FASTA format.
    J. --F1_Thresh, required=False: Threshold of the F1 score for selection of high performance spHMMs (Def.=0.5).
	K. --output_directory, required=True: Directory to save results. 
	L. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The high-performance spHMMs will be saved in the HiPer_spHMMs folder in the output directory specified. The HiPer_spHMMs folder should have the following files:
	```
    A. plot.eps : F1 score plot of all the spHMMs and the F1 cutoff threshold. 
    B. *.hmm : A set of spHMMs that perform above the F1 cutoff threshold.
    C. F1_Cutoff.txt: HMM search cutoff scores to be used for each high-performance spHMM interval. 
    	**Because synthetic datasets do not fully represent real data, please be aware that some of the spHMM cutoffs may need to be further tuned after running MetaBGC on a real metagenomic dataset, as was done with the Type II polyketide cyclase cutoffs in the original MetaBGC publication.
		**
    ```
    
### Running MetaBGC-Identify to detect biosynthetic reads

1. For identifying biosynthetic reads from the protein family of interest using the spHMMs constructed in MetaBGC-Build, the build step the MetaBGC-Identify.py script should be executed with the required input files. To use our pre-built high-performance spHMMs for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs). Please find the high performance input folders [here](https://github.com/donia-lab/MetaBGC/tree/master/MetaBGC-V1/MetaBGC-Build_Outputs).

	```
	A. --sphmm_directory, required=True: The high performance spHMM directory generated from MetaBGC-Build.
	B. --nucl_seq_directory, required=True: Directory of reads for the metagenomic samples to be analyzed. The filenames are used as sample names.  
    C. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    D. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    E. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    F. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	G. --prot_family_name, required=True: Name of the protein family.
	H. --cohort_name, required=True: Name of the cohort of metagenomic samples used for analysis.
	I. --output_directory, required=True: Directory to save results in.
	J. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```

2. MetaBGC-Identify will produce a FASTA file, **identified-biosynthetic-reads.fasta**, comprised of all biosynthetic reads identified in the metagenomic samples of the analyzed cohort, based on the specified cutoffs.

### Running MetaBGC-Quantify to profile biosynthetic-like reads

1. To quantify the biosynthetic reads identified by MetaBGC-Identify, the MetaBGC-Quantify.py script should be executed with the following parameters: 

	```
    A. --identify_fasta, required=True: Path to the CombinedIDFASTASeqs.fasta file produced by MetaBGC-Identify.
    B. --nucl_seq_directory, required=True: Directory of reads for the metagenomic samples to be analyzed. The filenames are used as sample names.  
    C. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    D. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    E. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    F. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	G. --cohort_name, required=True: Name of the cohort of metagenomic samples used for analysis.
	H. --output_directory, required=True: Directory to save results.
	I. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The output of the quantify script is an abundance profile file **unique-biosynthetic-reads-abundance-table.txt**. 

### Running MetaBGC-Cluster to generate BGCs bins 
1. To generate BGC bins of “unique biosynthetic reads”, users should use the abundance profile file, **unique-biosynthetic-reads-abundance-table.txt**, produced from the MetaBGC-Quantify module as input for [MetaBGC-Cluster.py](https://github.com/donia-lab/MetaBGC/tree/master/MetaBGC-V1/MetaBGC-V1_Scripts/MetaBGC-Cluster/MetaBGC-Cluster.py) with the following suggested parameters: `eps .2 min_samples 1` (these are suggested parameters and may have to be tuned depending on the specific dataset and protein family analyzed)

2. For synthetic datasets, we suggest examining bins that contain at least 50 reads, and for real datasets, we suggest examining bins that contain at least 10 reads (these are suggested parameters and may have to be tuned depending on the specific dataset and protein family analyzed). The resulting bins can be utilized in downstream analyses, such as targeted or untargeted assemblies to obtain the complete BGC, bin abundance calculations to determine the distribution of a given BGC in the entire cohort, etc. Please see the original MetaBGC publication for example analyses. 


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




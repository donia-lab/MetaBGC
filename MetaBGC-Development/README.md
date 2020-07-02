# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

[![DOI](https://zenodo.org/badge/187888803.svg)](https://zenodo.org/badge/latestdoi/187888803)

MetaBGC is a read-based algorithm for the detection of biosynthetic gene clusters (BGCs) directly in metagenomic sequencing data.

## Publication 

See [A metagenomic strategy for harnessing the chemical repertoire of the human microbiome](https://doi.org/10.1126/science.aax9176) for detailed description and for information on citing MetaBGC.

## Getting Started

These instructions will get you setup to run MetaBGC on your local Linux or Apple environment. 

## Important Notes

* Users can test the pipeline with the toy dataset provided in [https://github.com/donia-lab/MetaBGC#quick-start](https://github.com/donia-lab/MetaBGC#quick-start).
* Running your own read libraries against the test cyclase model provided may not work considering that cyclases are rare in metagenomic datasets.
* For non-cyclase spHMMs you can de novo create your own using the build module. In addition, we are currently developing high-performance spHMMs for several other biosynthetic classes and will be releasing them in a few months (as pre-built models) in a follow-up release/publication.
  
## Bioconda Distribution 

Coming soon...

## Manual Installation

### Prerequisites

To run MetaBGC, please make sure you have the following dependencies installed and in PATH.

* [Python](https://www.python.org/downloads/) (version >= 3.6)
* [MUSCLE 3.8.31](https://www.drive5.com/muscle/downloads.htm)
* [HMMER](http://hmmer.org/download.html) version 3.1b2
* [CD-HIT-EST](https://github.com/weizhongli/cdhit/releases) version 4.7
* [ncbi-blast-2.7.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
* [EMBOSS-6.5.7](http://emboss.sourceforge.net/)
* [R 3.6.1](http://lib.stat.cmu.edu/R/CRAN/): MetaBGC will try to install the following R packages at runtime. They can be installed manually if the user does not have access to the default R lib path.  
  * tidyverse
  * ggsci
  * ggpubr
  * dplyr
  * ggplot2

### Installing the Package

All general purpose users should obtain the package from PyPI:
```
pip install metabgc
```
All the internal python dependencies are specified in the setup will be installed. Help on the commandline parameters is provided. 
```
metabgc --help
```
The latest development build is available for other developers and may often contain bugs and is not stable. To install the latest development build from GitHub:
```
git clone https://github.com/donia-lab/MetaBGC.git
cd MetaBGC/MetaBGC-Development
pip install .
metabgc --help
```

### Quick Start

To run a toy search example, please download the already constructed spHMMs and some search samples from [here](https://drive.google.com/file/d/1-6nhr7WpWFbAxj89F-Q4fmB83nBhBDI-/view?usp=sharing) . To build your own spHMM database you will need to construct simulated read libraries as described in the publication [here](https://doi.org/10.1126/science.aax9176).  

```
OP_PATH=<set a path>
cd ${OP_PATH}
tar -zxvf OxyN_Build.tar.gz
metabgc search --sphmm_directory ${OP_PATH}/build/HiPer_spHMMs --prot_family_name Cyclase_OxyN --cohort_name OxyN --nucl_seq_directory ${OP_PATH}/build/nucl_seq_dir --seq_fmt FASTA --pair_fmt interleaved --output_directory ${OP_PATH}/output --cpu 4
```

The run will take about 60GB of memory and 220 minutes to run using 4 threads on a 2.5 GHz Ivybridge Intel processor. Please reduce the number of threads if you are running out of memory, but execution will take longer.  

### Program Structure

MetaBGC consists of four main modules:

**Build** - ```metabgc build --help``` - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. Pre-built high-performance spHMMs exist for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs).  If any of these protein families is to be used, this step can be skipped. The models are here: https://github.com/donia-lab/MetaBGC/tree/master/Models

**Identify** - ```metabgc identify --help``` -  This module runs on translated metagenomic reads from a cohort of samples using a selected set of high-performance spHMMs and their pre-set score cutoffs, as determined in MetaBGC-Build. The results are parsed into a list of identified biosynthetic reads in fasta format.

**Quantify** - ```metabgc quantify --help``` - This module de-replicates all biosynthetic reads discovered by MetaBGC-Identify from all metagenomic samples in the cohort into a unified set of **unique biosynthetic reads**. An abundance profile martrix is then generated for all unique biosynthetic reads by quantifying them in all samples of the metagenomic cohort.

**Cluster** - ```metabgc cluster --help``` - This module uses Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to cluster unique biosynthetic reads with similar abundance profiles across different metagenomic samples into distinct bins, based on the abundance profile martrix obtained in MetaBGC-Quantify. 

**Search** - ```metabgc search --help``` -**NEW** option to run Identify, Quantify, and Cluster together as a single command.  

### Running Build to construct the spHMMs

1. To build and evaluate spHMMs for the protein family of interest, the ```metabgc build``` command has to be executed with required input files. To select high performance spHMMs, a synthetic metagenomic dataset must be generated with reads from true positive genes spiked in to test the performance of each spHMM. To generate synthetic metagenomes for the build process, one can choose any publicly available tools, including ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm), wgsim (https://github.com/lh3/wgsim), etc. If using ART, as what we did in this paper, you can generate reads for individual genomes and concatenate the resulting reads to simulate a metagenome. To control the abundance of each genome, you need to control number of reads to generate from each genome. Users may also refer to pipelines such as CAMISIM (https://github.com/CAMI-challenge/CAMISIM) to facilitate this process.

	```
    1. --prot_alignment, required=True: Alignment of homologs from the protein family of interest in FASTA format.
	2. --prot_family_name, required=True: Name of the protein family. 
         This is used as prefix for spHMM files.
	3. --cohort_name, required=True: Name of the cohort of synthetic metagenomic samples used for evaluation.
	4. --nucl_seq_directory, required=True: Directory of reads for the synthetic metagenomic samples. The filenames are used as sample names.
	5. --prot_seq_directory, required=False: Directory with translated synthetic read files of the cohort. Computed if not provided.
    6. --seq_fmt, required=True: {fasta, fastq} Sequence file format and extension.
    7. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    8. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    9. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	10. --tp_genes_nucl, required=True: Nucleotide sequence of the full-length true-positive genes in the synthetic dataset in multi-FASTA format.
    11. --F1_Thresh, required=False: Threshold of the F1 score for selection of high performance spHMMs (Def.=0.5).
    12. --blastn_search_directory, required=False: Directory with BLAST search of the synthetic read files against the TP genes. Computed if not provided. To compute seperately, please see job_scripts in development.
    13. --hmm_search_directory, required=False: Directory with HMM searches of the synthetic read files against all the spHMMs. Computed if not provided. To compute seperately, please see job_scripts in development.
	14. --output_directory, required=True: Directory to save results. 
	15. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The high-performance spHMMs will be saved in the ```HiPer_spHMMs``` folder in the output directory specified. The ```HiPer_spHMMs``` folder should have the following files:

	```
    1. F1_Plot.eps : F1 score plot of all the spHMMs and the F1 cutoff threshold. 
    2. *.hmm : A set of spHMMs that perform above the F1 cutoff threshold.
    3. <prot_family_name>_F1_Cutoff.tsv: HMM search cutoff scores to be used for each high-performance spHMM interval.
    4. <prot_family_name>_Scores.tsv: FP, TP and FN scores of the the HMM search for all the spHMMs.
    5. <prot_family_name>_FP_Reads.tsv: False positive reads in the HMM search but not in the BLAST search. 
    ```
  	>**Because synthetic datasets do not fully represent real data, please be aware that some of the spHMM cutoffs may need to be further tuned after running MetaBGC on a real metagenomic dataset, as was done with the Type II polyketide cyclase cutoffs in the original MetaBGC publication.**
    	
    
### Running Identify to detect biosynthetic reads

1. For identifying biosynthetic reads from the protein family of interest using the spHMMs constructed in Build step, the ```metabgc identify``` command should be executed with the required input files. To use our pre-built high-performance spHMMs for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs), please find the high performance input folders [here](https://github.com/donia-lab/MetaBGC/tree/master/MetaBGC-V1/MetaBGC-Build_Outputs).

	```
	1. --sphmm_directory, required=True: The high performance spHMM directory generated from MetaBGC-Build.
	2. --nucl_seq_directory, required=True: Directory of reads for the metagenomic samples to be analyzed. The filenames are used as sample names.
	3. --prot_seq_directory, required=False: Directory with translated synthetic read files of the cohort. Computed if not provided.  
    4. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    5. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    6. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    7. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	8. --prot_family_name, required=True: Name of the protein family.
	9. --cohort_name, required=True: Name of the cohort of metagenomic samples used for analysis.
	10. --hmm_search_directory, required=False: Directory with HMM searches of the synthetic read files against all the spHMMs. Computed if not provided. To compute seperately, please see job_scripts in development.
	11. --output_directory, required=True: Directory to save results in.
	12. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```

2. Identify will produce a FASTA file, **identified-biosynthetic-reads.fasta**, comprised of all biosynthetic reads identified in the metagenomic samples of the analyzed cohort, based on the specified cutoffs.

### Running Quantify to de-replicate and quantify biosynthetic reads

1. To de-replicate and quantify the biosynthetic reads found by Identify, the ```metabgc quantify``` command should be executed with the following parameters: 

	```
    1. --identify_fasta, required=True: Path to the identified-biosynthetic-reads.fasta file produced by MetaBGC-Identify.
    2. --nucl_seq_directory, required=True: Directory of reads for the metagenomic samples to be analyzed. The filenames are used as sample names.  
    3. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    4. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    5. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    6. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	7. --cohort_name, required=True: Name of the cohort of metagenomic samples used for analysis.
	8. --blastn_search_directory, required=False: Directory with BLAST search of the synthetic read files against the TP genes. Computed if not provided. To compute seperately, please see job_scripts in development.
	9. --output_directory, required=True: Directory to save results.
	10. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The output of the quantify command is are abundance profile files:
 ```
    1. unique-biosynthetic-reads-abundance-table.txt : Contains the read level abundance matrix.
    2. unique-biosynthetic-reads-abundance-table-wide.txt : Contains the sample and read level abundance values.  
 ``` 

### Running Cluster to generate biosynthetic read bins
1. To generate BGC bins of **unique biosynthetic reads**, users should use the abundance profile file, **unique-biosynthetic-reads-abundance-table.txt**, produced by Quantify as input for ```metabgc cluster```. The input parameters to the script are:

	```
	1. --table, required=True: Path of tab-delimited abundance table,unique-biosynthetic-reads-abundance-table.txt, produced by MetaBGC-Quantify.
	2. --table_wide, required=True: Path of tab-delimited sample and read level abundance table,unique-biosynthetic-reads-abundance-table-wide.txt, produced by MetaBGC-Quantify.
	3. --identify_fasta, required=True: Path to the identified-biosynthetic-reads.fasta file produced by MetaBGC-Identify.
	4. --max_dist, required=False: Maximum Pearson distance between two reads to be in the same cluster. (Def.=0.1)
	5. --min_samples, required=False: Minimum number of samples required for a cluster. If min_samples > 1, noise are labelled as -1. (Def.=1).
	6. --min_reads_bin, required=False: Minimum number of reads required in a bin to be considered in analytics output files. (Def.=10)
	7. --min_abund_bin, required=False: Minimum total read abundance required in a bin to be considered in analytics output files. (Def.=10)
	8. --cpu, required=False: Number of CPU threads. (Def.=1)
	```

2. The cluster command produces the following files:
```
    1. unique-biosynthetic-reads-abundance-table_DBSCAN.json: Bin labels assigned DBSCAN, comprised of all the biosynthetic reads clustered in json format.
    2. BinSummary.txt : Summary of the bins and other statistics. 
    3. ReadLevelAbundance.tsv: Abundance of each biosynthetic read in each assigned bin with >= min_reads_bin reads. 
    4. SampleAbundanceMatrix.tsv : Abundance matrix of each sample against the bin ids for each bin with >= min_reads_bin reads.
    5. bin_fasta: Directory with FASTA files containing reads belonging each bin.  
```   

3. For synthetic datasets, we suggest examining bins that contain at least 50 reads, and for real datasets, we suggest examining bins that contain at least 10 reads (these are suggested parameters and may have to be tuned depending on the specific dataset and protein family analyzed). The resulting bins can be utilized in downstream analyses, such as targeted or untargeted assemblies to obtain the complete BGC, bin abundance calculations to determine the distribution of a given BGC in the entire cohort, etc. Please see the original MetaBGC publication for example analyses. 

### Running Search to detect biosynthetic read bins

1. For detecting biosynthetic reads from the protein family of interest using the spHMMs constructed in Build step, the ```metabgc search``` command should be executed with the required input files to run Identify, Quantify and Cluster in ine command. To use our pre-built high-performance spHMMs for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs), please find the high performance input folders [here](https://github.com/donia-lab/MetaBGC/tree/master/MetaBGC-V1/MetaBGC-Build_Outputs).

	```
	1. --sphmm_directory, required=True: The high performance spHMM directory generated from MetaBGC-Build.
	2. --nucl_seq_directory, required=True: Directory of reads for the metagenomic samples to be analyzed. The filenames are used as sample names.
	3. --prot_seq_directory, required=False: Directory with translated synthetic read files of the cohort. Computed if not provided.  
    4. --seq_fmt, required=True: {fasta,fastq} Sequence file format and extension.
    5. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    6. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    7. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	8. --prot_family_name, required=True: Name of the protein family.
	9. --cohort_name, required=True: Name of the cohort of metagenomic samples used for analysis.
	10. --blastn_search_directory, required=False: Directory with BLAST search of the synthetic read files against the TP genes. Computed if not provided. To compute seperately, please see job_scripts in development.
	11. --hmm_search_directory, required=False: Directory with HMM searches of the synthetic read files against all the spHMMs. Computed if not provided. To compute seperately, please see job_scripts in development.
    12. --max_dist, required=False: Maximum Pearson distance between two reads to be in the same cluster. (Def.=0.1)
	13. --min_samples, required=False: Minimum number of samples required for a cluster. If min_samples > 1, noise are labelled as -1. (Def.=1).
	14. --min_reads_bin, required=False: Minimum number of reads required in a bin to be considered in analytics output files. (Def.=10)
	15. --min_abund_bin, required=False: Minimum total read abundance required in a bin to be considered in analytics output files. (Def.=10)
	16. --output_directory, required=True: Directory to save results in.
	17. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```

2. Search will produce the same output as the Cluster command described above and all the intermediate files from each step. 

## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](https://github.com/donia-lab/MetaBGC/blob/master/MetaBGC-V1/LICENSE) file for details.

Please feel free to report bugs at our [github issue page](https://github.com/donia-lab/MetaBGC-TIIPKS/wiki).

## Contact
If there are any questions please contact:

**Mohamed Abou Donia, Ph.D.**
 <br /> Assistant Professor
 <br /> Department of Molecular Biology
 <br /> Princeton University
 <br /> Email: <donia@princeton.edu>

# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

[![PyPI](https://img.shields.io/pypi/v/metabgc.svg)](https://pypi.python.org/pypi/metabgc)
[![DOI](https://zenodo.org/badge/187888803.svg)](https://zenodo.org/badge/latestdoi/187888803)
[![PyPI download](https://img.shields.io/pypi/dm/metabgc?label=pypi%20downloads)](https://pypistats.org/packages/metabgc)
![PyPI Format](https://img.shields.io/pypi/format/metabgc)
![PyPI Python Version](https://img.shields.io/pypi/pyversions/metabgc)
![Programming languages](https://img.shields.io/github/languages/top/donia-lab/metabgc)
![License](https://img.shields.io/github/license/donia-lab/metabgc)


MetaBGC is a read-based algorithm for the detection of biosynthetic gene clusters (BGCs) directly in metagenomic sequencing data.

## Publication 

See [A metagenomic strategy for harnessing the chemical repertoire of the human microbiome](https://doi.org/10.1126/science.aax9176) for detailed description and for information on citing MetaBGC.

## Getting Started

These instructions will get you setup to run MetaBGC on your local Linux or Apple environment. 

## Important Notes

* **NEW** : All R dependency has been removed.
* **NEW** : Users can test the build pipeline with the toy dataset provided in [https://github.com/donia-lab/MetaBGC#quick-start](https://github.com/donia-lab/MetaBGC#quick-start).
* **NEW** : Docker container and Bioconda distribution is here!
* Users can test the search pipeline with the toy dataset provided in [https://github.com/donia-lab/MetaBGC#quick-start](https://github.com/donia-lab/MetaBGC#quick-start).
* Running your own read libraries against the test cyclase model provided may not work considering that cyclases are rare in metagenomic datasets.
* For non-cyclase spHMMs you can de novo create your own using the build module. In addition, we are currently developing high-performance spHMMs for several other biosynthetic classes and will be releasing them in a few months (as pre-built models) in a follow-up release/publication.

## Docker Containers
To install metabgc using Docker, please install Docker for your operating system. Once Docker is up follow the steps below:

### Download Source Code
Docker container files are provided for releases > 2.0.0. Go to the latest release on GitHub and download the source code tar ball and uncompress it.
```
wget https://github.com/donia-lab/MetaBGC/archive/2.0.0.tar.gz
tar -zxvf 2.0.0.tar.gz
```

### Build the Docker Image
Change to the source code directory and bild the container. Then you can run the container from commandline to view metabgc help. 
```
cd 2.0.0/MetaBGC-Development
docker build --tag metabgc .
```
 
## Bioconda Distribution 

MetaBGC can be installed using Bioconda too. The following commands install the dependencies and then install PyPI metabgc from PyPI:
```
conda create --name metabgc_env python=3.7
conda install muscle=3.8.31
conda install -c bioconda hmmer=3.1b2 cd-hit blast emboss=6.5.7
pip install metabgc
```

## Manual Installation 

### Prerequisites

To run MetaBGC, please make sure you have the following dependencies installed and in PATH.

* [Python](https://www.python.org/downloads/) (version >= 3.6)
* [MUSCLE 3.8.31](https://www.drive5.com/muscle/downloads.htm)
* [HMMER](http://hmmer.org/download.html) version 3.1b2
* [CD-HIT-EST](https://github.com/weizhongli/cdhit/releases) version 4.7
* [ncbi-blast-2.7.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
* [EMBOSS-6.5.7](http://emboss.sourceforge.net/)
* [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

### Installing the Package using pip3

All general purpose users should obtain the package from PyPI:
```
pip install metabgc
```
All the internal python dependencies are specified in the setup will be installed. Help on the commandline parameters is provided. 
```
metabgc --help
```

## Quick Start
The quick start guide provides 2 sample datasets for building a new database and searching an existing database. The example commands are provided to run the processes and can be used as a template.  
### Build a Database
To run a toy build example, please download a file from [here](https://drive.google.com/file/d/198bEperQcWhbl7ubn5OHq5wwppzXozkb). 

```
OP_PATH=<set a path>
cd ${OP_PATH}
tar -zxvf OxyN_Build_v3.tar.gz
mkdir output
metabgc build --prot_alignment Cyclase_OxyN.fasta --prot_family_name OxyN --cohort_name OxyN --nucl_seq_directory toy_reads --seq_fmt FASTA --pair_fmt interleaved --tp_genes_nucl TP_Genes.fasta --f1_thresh 0.5 --output_directory output --cpu 4
```
This should take about 1 hour using 8 threads on a 2.5 GHz Ivybridge Intel processor. The directory also has a SLURM job script ``` runBuildTest.sh```, that can be submitted to a cluster after changing the paths.  

### Perform a Search 
To run a toy search example, please download the already constructed spHMMs and some search samples from [here](https://drive.google.com/file/d/1-UG7aW72pDXUlvGeZWHEE49dtx2HSVzN) . To build your own spHMMs you will need to construct simulated read libraries using the ```metabgc synthesize``` command.  

```
OP_PATH=<set a path>
cd ${OP_PATH}
tar -zxvf OxyN_Search_v2.tar.gz
cd OxyN_Search_v3
metabgc search --sphmm_directory ${OP_PATH}/build/HiPer_spHMMs --prot_family_name Cyclase_OxyN --cohort_name MetaHit --nucl_seq_directory ${OP_PATH}/nucl_seq_dir --seq_fmt FASTA --pair_fmt interleaved --output_directory ${OP_PATH}/output --cpu 4
```
The run will take about 60GB of memory and 220 minutes to run using 4 threads on a 2.5 GHz Ivybridge Intel processor. Please reduce the number of threads if you are running out of memory, but execution will take longer.  

## Modules and Commandline Options 

MetaBGC consists of various modules to run the search pipeline. 

**Synthesize** - ```metabgc synthesize --help``` - This module creates a synthetic read dataset for building high-performance segmented profile Hidden Markov Models (spHMMs) using the build module.  

**FindTP** - ```metabgc findtp --help``` - This module finds the true-positive genes in the synthetic dataset genomes for the protein family of interest in multi-FASTA format to be used by the ```metabgc build module```.

**Build** - ```metabgc build --help``` - This module builds, evaluates, and selects high-performance segmented profile Hidden Markov Models (spHMMs) for a new protein family that is commonly found in the BGC of interest. Pre-built high-performance spHMMs exist for cyclases/aromatases commonly found in TII-PKS BGCs (OxyN, TcmN, TcmJ, and TcmI types), LanC_like proteins (found in lantibiotic BGCs), and IucA/IucC proteins (found in siderophore BGCs).  If any of these protein families is to be used, this step can be skipped. The models are here: https://github.com/donia-lab/MetaBGC/tree/master/Models. In addition, we are currently developing high-performance spHMMs for several other biosynthetic classes and will be releasing them in a few months (as pre-built models) in a follow-up release/publication. 

**Identify** - ```metabgc identify --help``` -  This module runs on translated metagenomic reads from a cohort of samples using a selected set of high-performance spHMMs and their pre-set score cutoffs, as determined in MetaBGC-Build. The results are parsed into a list of identified biosynthetic reads in fasta format.

**Quantify** - ```metabgc quantify --help``` - This module de-replicates all biosynthetic reads discovered by MetaBGC-Identify from all metagenomic samples in the cohort into a unified set of **unique biosynthetic reads**. An abundance profile martrix is then generated for all unique biosynthetic reads by quantifying them in all samples of the metagenomic cohort.

**Cluster** - ```metabgc cluster --help``` - This module uses Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to cluster unique biosynthetic reads with similar abundance profiles across different metagenomic samples into distinct bins, based on the abundance profile martrix obtained in MetaBGC-Quantify. 

**Search** - ```metabgc search --help``` - This is a combined option to run Identify, Quantify, and Cluster together as a single command.  

**Analytics** - ```metabgc analytics --help``` - This is a combined option to run final analytics as produced by the Donia Lab using metadata information of the cohorts. 


### Running Synthesize to Create a Synthetic Dataset 

1. Creating a synthetic dataset for running the build model requires a large set of background genomes, and a few genomes which are positive for the protein family of interest. To generate synthetic metagenomes for the build process, we use ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). The ```metabgc synthesize``` module uses ART to generate reads for individual genomes and concatenate the resulting reads to simulate a metagenome. To control the abundance of each genome, you need to control number of reads to generate from each genome. Users may also refer to pipelines such as CAMISIM (https://github.com/CAMI-challenge/CAMISIM) to facilitate this process.

2. The ```metabgc synthesize``` command has to be executed with required input files and the following parameters:

	```
	1. --indir1, required=True: Input directory of background genomes in FASTA for simulation.
	2. --indir2, required=True: Input directory protein-family-positive genomes in FASTA for simulation.
	3. --system, required=True: Illumina sequencing system (HS10, HS25, MSv3, etc.). Same options as -ss in art_illumina (Def.=HS20). 
	4. --length, required=True: Read length in bp (Def.=100). 
	5. --mflen, required=True: Mean fragment size in bp (Def.=400).
	6. --mflensd, required=True: Standard dev of fragment size in bp (Def.=20).
	7. --num_reads, required=True: Total number of read pairs (Def.=51100000).
	8. --samples, required=True: Number of samples to generate (Def.=70).
	9. --prop, required=True: Proportion of background genomes to draw from the input directory of background genomes prodided by --indir1. Should be between 0 and 1 (Def.=0.9).
	10. --output_directory, required=True: Directory to save results. 
	11. --base_name, required=False: Prefix of sample name (Def.='S'). 
	12. --cpu, required=False: Number of processes to run. Should be smaller than --samples (Def.=1).
	13. --seed, required=False: Random seed (Def.=915).
	```
 3. Generates sample metagenome read files in the output directory. 


### Running FindTP to find the True Positive Genes

1. The ```metabgc findtp``` finds TP genes in the protein family positive fasta files used for simulation. The command has to be executed with required input files and the following parameters:

	```
	1. --alnFile, required=True: Input protein family genomes in a FASTA format or a MUSCLE alignment. If it is a FASTA file then set --do_alignment True. 
	2. --prot_seq_directory, required=True: 
	3. --output_directory, required=True: Directory to save results.
	4. --do_alignment, required=True: Set it up to do an alignment if the --alnFile is a FASTA.
	```

2. The output produced are:
	```
	1. CombinedHMMHits.faa: It is a FASTA file with the TP sequences. 
	2. CombinedHMMResults.txt: The HMM score summary of the TP genes and genome it matches with the match score. 
	```


### Running Build to Construct the spHMMs

1. To build and evaluate spHMMs for the protein family of interest, the ```metabgc build``` command has to be executed with required input files. To select high performance spHMMs, a synthetic metagenomic dataset must be generated with reads from true positive genes spiked in to test the performance of each spHMM. To generate synthetic metagenomes for the build process use the ```metagbc synthesize module```.

	```
    1. --prot_alignment, required=True: Alignment of homologs from the protein family of interest in FASTA format.
	2. --prot_family_name, required=True: Name of the protein family. This is used as prefix for spHMM files.
	3. --cohort_name, required=True: Name of the cohort of synthetic metagenomic samples used for evaluation.
	4. --nucl_seq_directory, required=True: Directory of reads for the synthetic metagenomic samples. The filenames are used as sample names.
	5. --prot_seq_directory, required=False: Directory with translated synthetic read files of the cohort. Computed if not provided.
    6. --seq_fmt, required=True: {fasta, fastq} Sequence file format and extension.
    7. --pair_fmt, required=True: {single, split, interleaved} Paired-end information.
    8. --R1_file_suffix, required=False: Suffix including extension of the file name specifying the forward reads. Not specified for single or interleaved reads. Example: .R1.fastq
    9. --R2_file_suffix, required=False: Suffix including extension of the file name specifying the reverse reads. Not specified for single or interleaved reads. Example: .R2.fastq 
	10. --tp_genes_nucl, required=True: Nucleotide sequence of the full-length true-positive genes in the synthetic dataset in multi-FASTA format. This can be generated using the "metabgc findtp" module.
    11. --F1_Thresh, required=False: Threshold of the F1 score for selection of high performance spHMMs (Def.=0.5).
    12. --blastn_search_directory, required=False: Directory with BLAST search of the synthetic read files against the TP genes. Computed if not provided. To compute seperately, please see job_scripts in development.
    13. --hmm_search_directory, required=False: Directory with HMM searches of the synthetic read files against all the spHMMs. Computed if not provided. To compute seperately, please see job_scripts in development.
	14. --output_directory, required=True: Directory to save results. 
	15. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
2. The high-performance spHMMs will be saved in the ```HiPer_spHMMs``` folder in the output directory specified. The ```HiPer_spHMMs``` folder should have the following files:

	```
    1. F1_Plot.png : F1 score plot of all the spHMMs and the F1 cutoff threshold. 
    2. *.hmm : A set of spHMMs that perform above the F1 cutoff threshold.
    3. <prot_family_name>_F1_Cutoff.tsv: HMM search cutoff scores to be used for each high-performance spHMM interval.
    4. <prot_family_name>_Scores.tsv: FP, TP and FN read counts of the the HMM search for all the spHMMs.
    5. <prot_family_name>_FP_Reads.tsv: The false positive reads. These are reads identified in the spHMM search but are derived from the TP genes in the BLAST search. 
    ```
  	>**Because synthetic datasets do not fully represent real data, please be aware that some of the spHMM cutoffs may need to be further tuned after running MetaBGC on a real metagenomic dataset, as was done with the Type II polyketide cyclase cutoffs in the original MetaBGC publication.**
    	
    
### Running Identify to Detect Biosynthetic Reads

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


### Running Quantify to De-replicate and Quantify Biosynthetic Reads

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


### Running Cluster to Generate Biosynthetic Read Bins
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
	5. bin_fasta: Directory with FASTA files containing reads belonging to each bin.  
	```

3. For synthetic datasets, we suggest examining bins that contain at least 50 reads, and for real datasets, we suggest examining bins that contain at least 10 reads (these are suggested parameters and may have to be tuned depending on the specific dataset and protein family analyzed). The resulting bins can be utilized in downstream analyses, such as targeted or untargeted assemblies to obtain the complete BGC, bin prevalence calculations to determine the distribution of a given BGC in the entire cohort, etc. Please see the original MetaBGC publication for example analyzes. 


### Running Search to Perform Identify, Quantify and Cluster 

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


### Running Analytics to Produce Reports

1. The analytics module produces useful reports from the clustered bins by merging relevant cohort metadata with the reads in the clusters. The metadata gathered for an example cohort is available in the ```metabgc/metadata``` folder. The paths in these files are only relevant to Princeton Research Computing clusters. The users will need to download relevant files of these publicly available metagenome datasets and prepare similar metadata files. The publications from Donia Lab refer to the datasets.



2. The ```metabgc analytics``` command should be executed with the required parameters and input files:

	```
	1. --metabgc_output_dir, required=True: Output directory of a successful metabgc search run.
	2. --cohort_metadata_file, required=True: Cohort metadata file such as metabgc/metadata/combined_cohort_v2.csv. 
	3. --assembly_metadata_file, required=True: CSV file with sample assembly paths such as metabgc/metadata/DoniaSampleAssemblyFileList.csv.
	4. --output_directory, required=True: Directory to save results.
	5. --cpu, required=False: Number of CPU threads to use (Def.=4). 
	```
3. The reports produced by analytics modules are:

	```
	1. ReadLevelBinStats.tsv : The metadata information of each identified and clustered read. Information includes bin_id, sample_id, sunject_id, and BLAST based annotation etc.
	2. MaxBinSample.tsv: The maximum contributing sample of each bin. 
	3. BinSubjectCounts_Bodysite.tsv: The subject counts in each bin grouped by the sample collection sites.
	4. BinSubjectCounts_BodysiteAgg.tsv: The subject counts in each bin grouped by the major sample collection sites.
	5. BinSampleCounts_Bodysite.tsv: The sample counts in each bin grouped by the sample collection sites.
	6. BinSampleCounts_BodysiteAgg.tsv: The sample counts in each bin grouped by the major sample collection sites. 
	7. BodysiteBinVisitsAbundances.tsv: The abundance of each subject in each bin across all visits. 
	8. SampleAbundanceBodysite_Stacked.tsv: The read abundance of each subject in each bin. 
	```

## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](https://github.com/donia-lab/MetaBGC/blob/master/MetaBGC-V1/LICENSE) file for details.

Please feel free to report bugs at our [github issue page](https://github.com/donia-lab/MetaBGC-TIIPKS/wiki).

## Contact
If there are any questions please contact:

**Mohamed Abou Donia, Ph.D.**
 <br /> Associate Professor
 <br /> Department of Molecular Biology
 <br /> Princeton University
 <br /> Email: <donia@princeton.edu>

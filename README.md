# MetaBGC: Metagenomic identifier of Biosynthetic Gene Clusters 

MetaBGC is a read-based pipeline for the detection of biosynthetic gene clusters (BGCs) directly in human microbiome-derived metagenomic sequencing data, without the need for bacterial cultivation or sequencing. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

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

### Program Structure

MetaBGC conists of three main modules: 

**MetaBGC-Identify** - This modules includes running segmented profile Hidden Markov Models (spHMMs) on translated metagenomic reads and parsing the results to identify reads that are biosynthetic.

**MetaBGC-Quantify** - This module includes de-replicating the results from MetaBGC-Identify module and quantifying 

**MetaBGC-Cluster** - This module generates bins of reads 

### Running MetaBGC-Identify to detect biosynthetic-like reads

1. Translate metagenomic reads from nucleic acid to amino acid using EMBOSS:6.6.0.0 transeq tool in all six open reading frames: `-frame=6 - table=0 -sformat pearson`. 
2. Run `hmmsearch` command from HMMER on the translated metagenomic data from step 1 using the heuristic filters: `--F1 0.02 --F2 0.02 --F3 0.02`. 
3. Parse results from step 2 using the parser script. 
4. Filter the results with the cutoffs for each spHMM. 
5. Extract read ids per sample from filtered results and extract sequences for read ids to create a fasta file for each sample. 

### Running MetaBGC-Quantify to profile biosynthetic-like reads

1. Combine fasta files from MetaBGC-Identify module and de-replicate at 95% identity and 95% alignment coverage using CD-HIT- EST with the following parameters: `-c .95 -n 10 -d 0 -aS .95`. 
2. Quantify de-replicated fasta file using BLASTn against all sample metagenomes using the following parameters: `-task blastn -dust no -max_target_seqs 1000000 -perc_identity 95.0 -qcov_hsp_perc 50 window_size 11`. 
3. Produce an abundance profile using the results from step2 BLASTn. 

### Running MetaBGC-Cluster to generate BGCs
1. Abundance profiles from the MetaBGC-Quantify module are clustered using DBSCAN clustering method with Pearson correlation distance metric and the following parameters: `.eps 2 min_samples 1.`
2. Each bin from DBSCAN must contain at least 50 reads or more to count as a cluster to examine.


## Assembly of bins from human metagenomes

Samples with the highest abundance (as calculated by MetaBGC-Quantify) of a given bin were assembled using SPAdes 3.8 (12), and scaffolds of interest were identified using BLASTn with the bin reads as a query and the metagenomic scaffolds as a database. Contigs matching the bins were annotated by identifying the open reading frames in Geneious (28) and comparing their protein sequence to the non-redundant protein database (nr) on NCBI using BLASTp. When several or short scaffolds were produced by SPAdes, these scaffolds were extracted and used as seeds for targeted assemblies using PRICE (29).


## License

This project is licensed under the GNU General Public License V3 - see the [LICENSE.md](LICENSE.md) file for details.


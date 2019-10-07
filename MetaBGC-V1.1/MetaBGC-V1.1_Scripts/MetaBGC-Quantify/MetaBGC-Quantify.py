#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
from rpy2.robjects.vectors import StrVector
import os
import pandas as pd
import glob
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from Utils.Utils import PreProcessReadsPar
from Utils.Utils import RunBLASTNDirectoryPar
from Utils.Utils import runCDHit

CPU_THREADS=4

def combine_blast_results(blast_dir_path, outdir, cohort_name):
	filenames = [f for f in os.listdir(blast_dir_path) if os.isfile(os.join(blast_dir_path, f))]
	tabular_colnames = "sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue"
	df_colnames = tabular_colnames.split()
	list_of_dfs = [pd.read_csv(filename, names=df_colnames, header=None, delim_whitespace=True) for filename in filenames]
	for dataframe, filename in zip(list_of_dfs, filenames):
		dataframe['Sample'] = os.path.basename(filename).split(".txt")[0]
		dataframe['cohort'] = cohort_name
	combined_df = pd.concat(list_of_dfs, ignore_index=True)
	combinedBLAST = os.path.join(outdir,"CobminedQuantifyBLAST.txt")
	combined_df.to_csv(combinedBLAST, index=False, sep='\t', header=False)
	return combinedBLAST

def create_clustering_file(outdir,outfile):
	base = rpackages.importr('base')
	packageNames = ('tidyverse')
	utils = rpackages.importr('utils')
	utils.chooseCRANmirror(ind=1)
	packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
	if len(packnames_to_install) > 0:
		utils.install_packages(StrVector(packnames_to_install))
	tidyverse = rpackages.importr('tidyverse')
	robjects.r['options'](warn=-1)
	create_file = robjects.r('''
		function(results_file,outdir) {
			all_domains_blast_df <- read_tsv(results_file, col_names = F)
			names(all_domains_blast_df) <- c("sseqid", "slen","sstart", "send", "qseqid", "qlen", "qstart", "qend", "qcovs", "pident"," evalue", "Sample", "cohort")
			all_domains_blast_df_count <- all_domains_blast_df %>% group_by(Sample, qseqid) %>% count() %>% ungroup()
			all_domains_blast_df_count_table <- all_domains_blast_df_count %>% spread(., Sample, n, fill =0 )
			abundFile = file.path(outdir, "unique-biosynthetic-reads-abundance-table.txt")
			abundWideFile = file.path(outdir, "unique-biosynthetic-reads-abundance-table-wide.txt")
			write_tsv(all_domains_blast_df_count_table, abundFile, col_names = T)
			write_tsv(all_domains_blast_df_count, abundWideFile, col_names = T)
		}
		''')

	create_file(outfile,outdir)

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--identify_fasta', required=True, help= "Path to the file produced by MetaBGC-Identify.")
	parser.add_argument('--nucl_seq_directory', required=True, help= "Directory with nuleotide fasta files.")
	parser.add_argument('--seq_fmt', required=True, help="Sequence file format and extension.: {fasta,fastq}.")
	parser.add_argument('--pair_fmt', required=True, help="Sequence pair format: {single, split, interleaved}.")
	parser.add_argument('--R1_file_suffix', required=False,
						help="Suffix including extension of the file name specifying R1 reads. Not specified for single or interleaved reads.")
	parser.add_argument('--R2_file_suffix', required=False,
						help="Suffix including extension of the file name specifying R2 reads. Not specified for single or interleaved reads.")
	parser.add_argument('--cohort_name', required=True, help="Name of the sample/cohort.")
	parser.add_argument('--output_directory', required=True, help="Directory to save results.")
	parser.add_argument('--cpu', required=False, help="Number of threads. Def.: 4")
	args = parser.parse_args()

	if args.cpu is not None:
		CPU_THREADS = int(args.cpu)

	quant_op_dir = args.output_directory

	nucl_seq_directory = PreProcessReadsPar(args.nucl_seq_directory, args.seq_fmt, args.pair_fmt,
										 args.R1_file_suffix.strip(), args.R2_file_suffix.strip(),
										 quant_op_dir, CPU_THREADS)

	cdHitFile = os.path.join(quant_op_dir,"CombinedIDFASTASeqs_Drep.fasta")
	runCDHit(args.identify_fasta,cdHitFile,CPU_THREADS)

	blastn_search_directory = os.path.join(quant_op_dir, 'quantify_blastn_result')
	os.makedirs(blastn_search_directory, 0o777, True)
	RunBLASTNDirectoryPar(nucl_seq_directory, cdHitFile, blastn_search_directory,CPU_THREADS)
	combinedBLASTPath = combine_blast_results(blastn_search_directory, quant_op_dir, args.cohort_name)
	create_clustering_file(quant_op_dir, combinedBLASTPath)

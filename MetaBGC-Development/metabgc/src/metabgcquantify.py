#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
#from rpy2.robjects.vectors import StrVector
#import rpy2.robjects.packages as rpackages
#import rpy2.robjects as robjects
from metabgc.src.utils import *

CPU_THREADS=4

def combine_blast_results(blast_dir_path, outdir, cohort_name):
	filenames = [os.path.join(blast_dir_path, f) for f in os.listdir(blast_dir_path) if os.path.isfile(os.path.join(blast_dir_path, f)) and os.path.getsize(os.path.join(blast_dir_path, f)) > 0]
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

def create_clustering_file(outdir,blast_result):
	rpackages.importr('base')
	utils = rpackages.importr('utils')
	packageNames = ('tidyverse')
	packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
	if len(packnames_to_install) > 0:
		utils.install_packages(StrVector(packnames_to_install))
	rpackages.importr('tidyverse')
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

	create_file(blast_result,outdir)
	return os.path.join(outdir,"unique-biosynthetic-reads-abundance-table.txt")


def mbgcquantify(identify_fasta, prot_family_name, cohort_name, nucl_seq_directory,
             seq_fmt, pair_fmt, r1_file_suffix, r2_file_suffix,
             output_directory, cpu):

	if cpu is not None:
		CPU_THREADS = int(cpu)
	if r1_file_suffix is None:
		r1_file_suffix = ""
	if r2_file_suffix is None:
		r2_file_suffix = ""

	quant_op_dir = output_directory
	nucl_seq_directory = PreProcessReadsPar(nucl_seq_directory, seq_fmt, pair_fmt,
										 r1_file_suffix.strip(), r2_file_suffix.strip(),
										 quant_op_dir, CPU_THREADS)
	cdHitFile = os.path.join(quant_op_dir,"CombinedIDFASTASeqs_Drep.fasta")
	runCDHit(identify_fasta,cdHitFile,CPU_THREADS)
	blastn_search_directory = os.path.join(quant_op_dir, 'quantify_blastn_result')
	os.makedirs(blastn_search_directory, 0o777, True)
	RunBLASTNDirectoryPar(nucl_seq_directory, cdHitFile, blastn_search_directory,CPU_THREADS)
	combinedBLASTPath = combine_blast_results(blastn_search_directory, quant_op_dir, cohort_name)
	abund_file = create_clustering_file(quant_op_dir, combinedBLASTPath)
	return abund_file
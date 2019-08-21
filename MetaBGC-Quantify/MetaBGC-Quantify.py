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


def combine_blast_results(blast_extension, sample_extension, blast_dirpath, tabular_colnames, outfile, outdir, cohort_name):
	#os.chdir(blast_dirpath)
	filenames = glob.glob(os.path.join(blast_dirpath, '*{}').format(blast_extension))
	#filenames = glob.glob("*blast_extension")

	list_of_dfs = [pd.read_csv(filename, names=tabular_colnames, header=None, delim_whitespace=True) for filename in filenames]
	for dataframe, filename in zip(list_of_dfs, filenames):
		dataframe['Sample'] = os.path.basename(filename).split(sample_extension)[0]
		dataframe['cohort'] = cohort_name

	combined_df = pd.concat(list_of_dfs, ignore_index=True)
	os.chdir(outdir)
	combined_df.to_csv(outfile, index=False, sep='\t', header=False)

def create_clustering_file(outdir, outfile ):
	os.chdir(outdir)
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
		function(results_file) {
			all_domains_blast_df <- read_tsv(results_file, col_names = F)
			names(all_domains_blast_df) <- c("sseqid", "slen","sstart", "send", "qseqid", "qlen", "qstart", "qend", "qcovs", "pident"," evalue", "Sample", "cohort")
			all_domains_blast_df_count <- all_domains_blast_df %>% group_by(Sample, qseqid) %>% count() %>% ungroup()
			all_domains_blast_df_count_table <- all_domains_blast_df_count %>% spread(., Sample, n, fill =0 )
			write_tsv(all_domains_blast_df_count_table, "abundance_table.txt", col_names = T)
			write_tsv(all_domains_blast_df_count, "abundance_table-wide.txt", col_names = T)

		}
		''')

	create_file(outfile)

def main(blast_extension, sample_extension, blast_dirpath, tabular_colnames, outfile, outdir, cohort_name):
	df_colnames = tabular_colnames.split()
	combine_blast_results(blast_extension, sample_extension, blast_dirpath, df_colnames, outfile, outdir, cohort_name)
	create_clustering_file(outdir, outfile)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--blast_extension', required=True, help= "extension to combine BLAST results")  
	parser.add_argument('--sample_extension', required=True, help= "extension to parse sample name")
	parser.add_argument('--blast_dirpath', required=True, help= "output directory name")
	parser.add_argument('--tabular_colnames', nargs='+', required=False, default = "sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue") 
	parser.add_argument('--outfile', required=False, default='combined-blast_quantifier-results.txt')
	parser.add_argument('--outdir', required=True, help="directory to output results")
	parser.add_argument('--cohort_name', required=True, help="dataset_name")


	args = parser.parse_args()

	main(args.blast_extension, args.sample_extension, args.blast_dirpath, args.tabular_colnames, args.outfile , args.outdir, args.cohort_name)

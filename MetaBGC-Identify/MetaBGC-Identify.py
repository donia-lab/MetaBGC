#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at fcamacho@princeton.edu).
#####################################################################################
from rpy2.robjects.vectors import StrVector
import os
import pandas as pd
import glob
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

#Filter HMM results using predetermined spHMM models score cutoffs 
def filter_spHMM_data(spHMM_df, cutoff_df):

	filter_spHMM_df = pd.DataFrame()
	for index, row in cutoff_df.iterrows():
		cyclase_model = cutoff_df.at[index,'cyclase_type']
		model_interval = cutoff_df.at[index,'interval']
		cutoff = cutoff_df.at[index,'cutoff']
		#filter spHMM dataframe with designated cutoffs 
		filter_df = spHMM_df.query("HMMScore >= @cutoff & cyclaseType == @cyclase_model & interval == @model_interval").reset_index(drop=True)
		filter_spHMM_df = filter_spHMM_df.append(filter_df)
	return(filter_spHMM_df)

#reformat IDs and filter duplicate readIDs and keep the highest HMM Score 
def create_reformat_data(input_df, outdir):
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
	reformat_df = robjects.r('''
		function(hmmdf) {
			hmmdfRecoded <- separate(hmmdf, readID, into = c("read","F_R_read_frame"), sep = "_", extra = "merge") %>%
			select(-c(F_R_read_frame))
			hmmdfRecodedDFUnique<-aggregate(HMMScore ~ read + Sample + sampleType + cyclaseType , hmmdfRecoded, max)
			colnames(hmmdfRecodedDFUnique)<-c("readID","Sample", "sampleType", "cyclaseType","HMMScore")
			write_tsv(hmmdfRecodedDFUnique, "spHMM-filtered-results.txt", col_names = T)
		}
		''')
	input_r_df = pandas2ri.py2ri(input_df) # convert pandas df to R datafame

	data_filter = reformat_df(input_r_df)
	return(data_filter)

def parse_sample_reads(reformat_df, reads_outdir):
	os.chdir(reads_outdir)
	parseReads = robjects.r('''
		function(HMMdf) {
			samples<-unique(HMMdf$Sample)
			for (s in 1:length(samples)){
				currentSample<-samples[s]
				currentSampleResults<- HMMdf %>% filter(Sample == currentSample)
				currentSampleReads<- unique(currentSampleResults$readID)
				fileName<-paste0(currentSample,paste("-","detected-reads",sep =""), ".txt", sep ="")
				write.table(currentSampleReads,fileName, quote = F, row.names = F, col.names = F )
		}}
		''')
	parseReads(reformat_df)

def main(hmm_file, outdir, cutoff_file, fasta_dir):

	spHMM_df = pd.read_csv(hmm_file, sep ="\t", names = ["readID", "sampleType", "Sample", "cyclaseType", "HMMScore", "window","interval"])
	cutoff_df = pd.read_csv(cutoff_file, sep ="\t", header=0)
	spHMM_df_filtered = filter_spHMM_data(spHMM_df, cutoff_df)
	spHMM_df_filtered_reformat = create_reformat_data(spHMM_df_filtered, outdir)
	parse_sample_reads(spHMM_df_filtered_reformat, fasta_dir)

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--hmm_file', required=True, help= "hmm results file")  
	parser.add_argument('--outdir', required=True, help= "directory to save results")
	parser.add_argument('--cutoff_file', required=True, help= "metadata file for score cutoffs")
	parser.add_argument('--fasta_dir', required=True, help= "directory to save sample reads fasta files")

	args = parser.parse_args()

	main(args.hmm_file, args.outdir, args.cutoff_file, args.fasta_dir)
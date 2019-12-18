#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
from metabgc.src.extractfastaseq import ExtractFASTASeq
from metabgc.src.utils import *
from rpy2.robjects.vectors import StrVector
import os
import pandas as pd
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from pathlib import Path

#Filter HMM results using predetermined spHMM models score cutoffs 
def filter_spHMM_data(spHMM_df, cutoff_df):
	filter_spHMM_df = pd.DataFrame()
	for index, row in cutoff_df.iterrows():
		model_interval = cutoff_df.at[index,'interval']
		cutoff = cutoff_df.at[index,'cutoff']
		#filter spHMM dataframe with designated cutoffs 
		filter_df = spHMM_df.query("HMMScore >= @cutoff & interval == @model_interval").reset_index(drop=True)
		filter_spHMM_df = filter_spHMM_df.append(filter_df)
	return(filter_spHMM_df)

#reformat IDs and filter duplicate readIDs and keep the highest HMM Score 
def create_reformat_data(input_df, outdir):
	rpackages.importr('base')
	utils = rpackages.importr('utils')
	packageNames = ('tidyverse')
	packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
	if len(packnames_to_install) > 0:
		utils.install_packages(StrVector(packnames_to_install))
	rpackages.importr('tidyverse')

	robjects.r['options'](warn=-1)
	reformat_df = robjects.r('''
		function(hmmdf,outDir) {
			hmmdfRecoded <- separate(hmmdf, readID, into = c("read","F_R_read_frame"), sep = "_", extra = "merge") %>%
			select(-c(F_R_read_frame))
			hmmdfRecodedDFUnique<-aggregate(HMMScore ~ read + Sample + sampleType + protType , hmmdfRecoded, max)
			colnames(hmmdfRecodedDFUnique)<-c("readID","Sample", "sampleType", "protType","HMMScore")
			write_tsv(hmmdfRecodedDFUnique, file.path(outDir, "spHMM-filtered-results.txt"), col_names = T)
			return(hmmdfRecodedDFUnique)
		}
		''')
	# convert pandas df to R datafame
	with localconverter(ro.default_converter + pandas2ri.converter):
		input_r_df = ro.conversion.py2rpy(input_df)

	data_filter = reformat_df(input_r_df,outdir)
	return(data_filter)

def parse_sample_reads(reformat_df, reads_outdir):
	parseReads = robjects.r('''
		function(HMMdf,reads_outdir) {
			samples<-unique(HMMdf$Sample)
			for (s in 1:length(samples)){
				currentSample<-samples[s]
				currentSampleResults<- HMMdf %>% filter(Sample == currentSample)
				currentSampleReads<- unique(currentSampleResults$readID)
				fileName<-paste0(currentSample,paste("-","detected-reads",sep =""), ".txt", sep ="")
				parsedFileName <- file.path(reads_outdir, fileName)
				write.table(currentSampleReads,parsedFileName, quote = F, row.names = F, col.names = F )
		}}
		''')
	parseReads(reformat_df,reads_outdir)

def runidentify(hmm_file, outdir, cutoff_file, fasta_dir):
	spHMM_df = pd.read_csv(hmm_file, sep ="\t", names = ["readID", "sampleType", "Sample", "protType", "HMMScore", "window","interval"])
	cutoff_df = pd.read_csv(cutoff_file, sep ="\t", header=0)
	spHMM_df_filtered = filter_spHMM_data(spHMM_df, cutoff_df)
	spHMM_df_filtered_reformat = create_reformat_data(spHMM_df_filtered, outdir)
	parse_sample_reads(spHMM_df_filtered_reformat, fasta_dir)


def mbgcidentify(sphmm_directory, cohort_name, nucl_seq_directory,prot_seq_directory,
					seq_fmt, pair_fmt, r1_file_suffix, r2_file_suffix,
					prot_family_name, hmm_search_directory, output_directory, cpu):

	if cpu is not None:
		CPU_THREADS = int(cpu)

	identify_op_dir = output_directory

	if r1_file_suffix is None:
		r1_file_suffix = ""
	if r2_file_suffix is None:
		r2_file_suffix = ""

	if hmm_search_directory is None:
		hmm_search_directory = os.path.join(identify_op_dir, 'hmm_identify_search')
	allHMMResult = hmm_search_directory + os.sep + "CombinedHmmSearch.txt"
	identify_directory = os.path.join(identify_op_dir, 'identify_result')
	fasta_dir = os.path.join(identify_op_dir, 'fasta_result')
	identifyReadIds = fasta_dir + os.sep + "CombinedReadIds.txt"
	multiFastaFile = identify_op_dir + os.sep + "identified-biosynthetic-reads.fasta"

	nucl_seq_directory = PreProcessReadsPar(nucl_seq_directory,
											seq_fmt, pair_fmt,
											r1_file_suffix.strip(),
											r2_file_suffix.strip(),
											identify_op_dir,
											CPU_THREADS)

	# Translate nucleotide seq
	if not os.path.isdir(prot_seq_directory):
		prot_seq_directory = TranseqReadsDir(identify_op_dir, nucl_seq_directory, CPU_THREADS)

	# HMMER search
	os.makedirs(hmm_search_directory, 0o777, True)
	for filename in os.listdir(sphmm_directory):
		fileBase = Path(filename).resolve().stem
		if filename.endswith(".hmm"):
			hmmInterval = fileBase.split("__")[2]
			hmmfilename = os.path.join(sphmm_directory,filename)
			RunHMMDirectory(prot_seq_directory, hmmfilename, cohort_name, prot_family_name, "30_10", hmmInterval,
						hmm_search_directory, CPU_THREADS)


	with open(allHMMResult, 'w') as outfile:
		for subdir, dirs, files in os.walk(hmm_search_directory):
			for file in files:
				filePath = os.path.join(subdir, file)
				if re.match(r".*txt$", file) and os.path.getsize(filePath) > 0:
					with open(filePath) as infile:
						for line in infile:
							outfile.write(line)


	os.makedirs(identify_directory, 0o777, True)
	os.makedirs(fasta_dir, 0o777, True)
	cutoff_file = os.path.join(sphmm_directory, prot_family_name + "_F1_Cutoff.tsv")
	runidentify(allHMMResult, identify_directory, cutoff_file, fasta_dir)


	with open(identifyReadIds, 'w') as outfile:
		for filename in os.listdir(fasta_dir):
			if filename.endswith(".txt"):
				filePath = os.path.join(fasta_dir, filename)
				with open(filePath) as infile:
					for line in infile:
						outfile.write(line)

	ExtractFASTASeq(nucl_seq_directory,multiFastaFile,identifyReadIds)
	return multiFastaFile


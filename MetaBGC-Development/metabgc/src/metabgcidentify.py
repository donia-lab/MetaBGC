#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
from metabgc.src.extractfastaseq import RunExtractDirectoryPar
from metabgc.src.utils import *
import os
import pandas as pd
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

def runidentify(hmm_file, cutoff_file,filteredTableFile,identifyReadIdFile):
	spHMM_df = pd.read_csv(hmm_file, sep ="\t", names = ["readID", "sampleType", "Sample", "protType", "HMMScore", "window","interval"])

	# Reformat IDs to remove frame identifier from transeq
	for i, row in spHMM_df.iterrows():
		translated_read_id = spHMM_df.at[i, 'readID']
		if translated_read_id.rfind('_') != -1:
			nucl_read_id = translated_read_id[:translated_read_id.rfind('_')]
			spHMM_df.at[i, 'readID'] = nucl_read_id

	cutoff_df = pd.read_csv(cutoff_file, sep ="\t", header=0)
	spHMM_df_filtered = filter_spHMM_data(spHMM_df, cutoff_df)

	# Filter duplicate readIDs and keep the highest HMM Score
	spHMM_df_filtered_uniq = spHMM_df_filtered.groupby(['readID','Sample','sampleType','protType'])['HMMScore'].max().reset_index()
	spHMM_df_filtered_uniq_interval = pd.merge(spHMM_df_filtered, spHMM_df_filtered_uniq, how='inner')

	spHMM_df_filtered_uniq_interval.to_csv(filteredTableFile, sep="\t",index=False)
	identifyReadIdList = list(set(spHMM_df_filtered_uniq_interval.readID.values.tolist()))

	with open(identifyReadIdFile, 'w') as outfile:
		for readID in identifyReadIdList:
			outfile.write(readID + '\n')

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
	fasta_seq_dir = os.path.join(identify_op_dir, 'fasta_seq_result')
	identifyReadIds = identify_directory + os.sep + "CombinedReadIds.txt"
	filteredHMMResult = identify_directory + os.sep + "spHMM-filtered-results.txt"
	multiFastaFile = identify_directory + os.sep + "identified-biosynthetic-reads.fasta"

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
	cutoff_file = os.path.join(sphmm_directory, prot_family_name + "_F1_Cutoff.tsv")
	runidentify(allHMMResult, cutoff_file, filteredHMMResult, identifyReadIds)

	os.makedirs(fasta_seq_dir, 0o777, True)
	RunExtractDirectoryPar(nucl_seq_directory, filteredHMMResult, fasta_seq_dir, multiFastaFile, CPU_THREADS)
	return multiFastaFile


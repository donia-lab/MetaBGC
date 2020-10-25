#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

from Bio import SearchIO
import os
import pandas as pd


class readID():
	def __init__(self, acc_numb, sampleType, sampleID, protType, bitscore, window, interval):
		self.acc_numb = acc_numb
		self.sampleType = sampleType
		self.sampleID = sampleID
		self.protType = protType
		self.bitscore = bitscore
		self.window = window
		self.interval = interval

def parseHMM(hmmFile, hmmDir, sampleType, sampleID, protType, window, interval):
	with open(os.path.join(hmmDir,hmmFile), 'rU') as handle:
		results_dict = {}
		try:
			for record in SearchIO.parse(handle, 'hmmer3-tab'):

				hits = record.hits
				num_hits = len(hits)  #calculate how many hits per query
				if num_hits > 0:  #extract hits data 
					for i in range(0, num_hits):
						hmm_name = hits[i].id
						hmm_bitscore = hits[i].bitscore
						if hmm_name not in results_dict: # add hits to results dictionary 
							readID_Class = readID(hmm_name, sampleType, sampleID, protType, hmm_bitscore, window, interval)
							results_dict[hmm_name] = readID_Class
						# the bio module wont allow duplicated query IDs
			handle.close()
		except ValueError:
			print ("ERROR: duplicated queryIDs in hmmer result file:", hmmFile)
	return results_dict

def createPandaDF(cyclase_dict, outdir, outfile):
	os.chdir(outdir)
	outputDF = pd.DataFrame()
	if len(cyclase_dict) != 0:  # check if counter dict is not empty to add to df
		df = [(k, v.sampleType, v.sampleID, v.protType, v.bitscore, v.window, v.interval) for k, v in
			  list(cyclase_dict.items())]  # convert dictionary to list

		outputDF = outputDF.append(df)  # append list to our initalized dataframe
		outputDF.columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window","interval"]  # rename column names
		sorted_outputDF = outputDF.sort_values(by=['HMMScore'],
											   ascending=[False])  # sort in decending order by Hit.counts column
		sorted_outputDF.to_csv(outfile, index=False, sep='\t', header=False)
	else:
		outputDF_empty_columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window", "interval"]  # rename column names
		outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
		outputDF_empty.to_csv(outfile, index=False, sep='\t', header=False)

def main(hmmscan_file_dir, outdir, cyclase_type, window, sampleType):
	os.makedirs(outdir,0o777,True)
	for hmmscan_file in os.listdir(hmmscan_file_dir):
		baseOutFile = hmmscan_file.split('.')[0]
		outfile = baseOutFile + "-PARSED.txt"
		sample_tok = baseOutFile.split('_')
		interval = sample_tok[-2] + '_' + sample_tok[-1]
		sample_tok.pop()
		sample_tok.pop()
		sampleID = '_'.join(sample_tok)
		protType = cyclase_type
		result_cyclase_Dict = parseHMM(hmmscan_file, hmmscan_file_dir, sampleType, sampleID, protType, window, interval)
		createPandaDF(result_cyclase_Dict, outdir, outfile)
if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--hmmscan_file_dir', required=True, help= "hmm results directory PATH")  
	parser.add_argument('--outdir', required=True, help= "output directory name")
	parser.add_argument('--cyclase_type',required=True, help="cyclase model name")
	parser.add_argument('--window', required=False, default='30_10')
	# parser.add_argument('--interval', required=True, help="spHMM interval")
	# parser.add_argument('--sampleID', required=True, help="sample name")
	parser.add_argument('--sampleType', required=True, help="type of sample, bodysite, isolation_source")


	args = parser.parse_args()

	main(args.hmmscan_file_dir, args.outdir, args.cyclase_type, args.window, args.sampleType)

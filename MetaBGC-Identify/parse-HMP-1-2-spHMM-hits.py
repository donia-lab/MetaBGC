#!/usr/bin/env python

from Bio import SearchIO
import os
import pandas as pd


class readID():
	def __init__(self, acc_numb, sampleType, sampleID, cyclaseType, bitscore, window, interval):
		self.acc_numb = acc_numb
		self.sampleType = sampleType
		self.sampleID = sampleID
		self.cyclaseType = cyclaseType
		self.bitscore = bitscore
		self.window = window
		self.interval = interval

def parseHMM(hmmFile, hmmDir, sampleType, sampleID, cyclaseType, window, interval):
	with open(os.path.join(hmmDir,hmmFile), 'rU') as handle:

		cyclaseHit_dict = {}
		try:
			for record in SearchIO.parse(handle, 'hmmer3-text'):

				hits = record.hits
				num_hits = len(hits)  # how many fattyacid hits per queries?
				if num_hits > 0:  # if there are more than 0 hits per query then we need to extract the info
					for i in range(0, num_hits):
						hmm_name = hits[i].id
						hmm_bitscore = hits[i].bitscore
						if hmm_name not in cyclaseHit_dict:  # if read id (PF01135.16) not in the dictionary add it and initalize the
							readIDClass = readID(hmm_name, sampleType, sampleID, cyclaseType, hmm_bitscore, window, interval)

							cyclaseHit_dict[hmm_name] = readIDClass
						# no else here because the bio module wont allow duplicated query IDs
			handle.close()

		except ValueError:
			print ("ERROR: duplicated queryIDs in hmmer result file:", hmmFile)

	return cyclaseHit_dict

def createPandaDF(cyclase_dict, outdir, outfile):
	os.chdir(outdir)
	outputDF = pd.DataFrame()
	if len(cyclase_dict) != 0:  # check if counter dict is not empty to add to df
		df = [(k, v.sampleType, v.sampleID, v.cyclaseType, v.bitscore, v.window, v.interval) for k, v in
			  list(cyclase_dict.items())]  # convert dictionary to list

		outputDF = outputDF.append(df)  # append list to our initalized dataframe
		outputDF.columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMM.Score", "window","interval"]  # rename column names
		sorted_outputDF = outputDF.sort_values(by=['HMM.Score'],
											   ascending=[False])  # sort in decending order by Hit.counts column
		sorted_outputDF.to_csv(outfile, index=False, sep='\t', header=False)  # write dataframe to csv format (text file)
	else:
		outputDF_empty_columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMM.Score", "window", "interval"]  # rename column names
		outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
		outputDF_empty.to_csv(outfile, index=False, sep='\t', header=False)  # write dataframe to csv format (text file)

def main(hmmscan_file_dir, outdir, cyclase_type, window, interval):
	os.makedirs(outdir)
	for hmmscan_file in os.listdir(hmmscan_file_dir):
		baseOutFile = hmmscan_file.split('.txt')[0]
		outfile = baseOutFile + "-PARSED.txt"
		metadata = str.split(baseOutFile, "-against") # parse file name to retrieve metadata on file, sample, sampletype
		sampleID = metadata[0]
		sampleType = metadata[0]
		cyclaseType = cyclase_type
		result_cyclase_Dict = parseHMM(hmmscan_file, hmmscan_file_dir, sampleType, sampleID, cyclaseType, window, interval)

		createPandaDF(result_cyclase_Dict, outdir, outfile)
if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--hmmscan_file_dir', required=True) # --hmm file dir 
	parser.add_argument('--outdir', required=True) #output directory 
	parser.add_argument('--cyclase_type',required=True)
	parser.add_argument('--window', required=False, default='30_10')
	parser.add_argument('--interval', required=True)


	args = parser.parse_args()

	main( args.hmmscan_file_dir, args.outdir, args.cyclase_type, args.window, args.interval)

#!/usr/bin/env python

from Bio import SearchIO
import logging
import pandas as pd


class readID():
	def __init__(self, acc_numb, sampleType, sampleID, cyclaseType, bitscore):
		self.acc_numb = acc_numb
		self.sampleType = sampleType
		self.sampleID = sampleID
		self.cyclaseType = cyclaseType
		self.bitscore = bitscore


def parseHMM(hmmfile, sampleType, sampleID, cyclaseType):
	with open(hmmfile, 'rU') as handle:
		polyketideType_dict = {}
		try:
			for record in SearchIO.parse(handle, 'hmmer3-text'):

				hits = record.hits
				num_hits = len(hits)  # how many fattyacid hits per queries?
				if num_hits > 0:  # if there are more than 0 hits per query then we need to extract the info
					for i in range(0, num_hits):
						hmm_name = hits[i].id
						hmm_bitscore = hits[i].bitscore
						if hmm_name not in polyketideType_dict:  # if read id (PF01135.16) not in the dictionary add it and initalize the
							readIDClass = readID(hmm_name, sampleType, sampleID, cyclaseType, hmm_bitscore)
							polyketideType_dict[hmm_name] = readIDClass
						# nothing else here because the bio module wont allow duplicated query IDs
			handle.close()
		except ValueError:
			print ("ERROR: duplicated queryIDs in hmmer result file:", hmmfile)
	return polyketideType_dict

def createPandaDF(polyketideTypeDict, hmmscore_cutoff, outfile ):
	outputDF = pd.DataFrame()
	if len(polyketideTypeDict) != 0:  # check if counter dict is not empty to add to df
		df = [(k, v.sampleType, v.sampleID, v.cyclaseType, v.bitscore) for k, v in
			  list(polyketideTypeDict.items())]  # convert dictionary to list
		outputDF = outputDF.append(df)  # append list to our initalized dataframe
		outputDF.columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMM.Score"]  # rename column names

		if hmmscore_cutoff > 0: # appply hmmscore cutoff if numerical value is greater than 0 
			filteredHMMDF =  outputDF[outputDF['HMM.Score'] >= hmmscore_cutoff]
		else:
			filteredHMMDF = outputDF # no filter applied 


		sorted_outputDF = filteredHMMDF.sort_values(by=['HMM.Score'],
											   ascending=[False])  # sort in decending order by Hit.counts column
		sorted_outputDF.to_csv(outfile, index=False, sep='\t')  # write dataframe to csv format (text file)
	else:
		outputDF_empty_columns = ["readID", "sampleType", "sampleID", "cyclaseType", "HMM.Score"]  # rename column names
		outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
		outputDF_empty.to_csv(outfile, index=False, sep='\t')  # write dataframe to csv format (text file)


def main(sampleID, sampleType, hmmfile, outfile, cyclaseType, hmmscore_cutoff):

	result_polyketide_Dict = parseHMM(hmmfile, sampleType, sampleID, cyclaseType)
	createPandaDF(result_polyketide_Dict, hmmscore_cutoff, outfile)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--sampleID', required=False, default="sample1")
	parser.add_argument('--sampleType', required=False, default = "metagenomic") 
	parser.add_argument('--hmmfile', required=True)
	parser.add_argument('--outfile', required=True)
	parser.add_argument('--cyclaseType', required=False, default="TypeIIPKS")
	parser.add_argument('--hmmscore_cutoff', required=False, default=0, type = int)
	parser.add_argument('-v', '--verbose', action='store_true')



	args = parser.parse_args()

	logger = logging.getLogger('hmm-result_parser')
	if args.verbose:
		logLevel = getattr(logging, 'DEBUG')
	else:
		logLevel = getattr(logging, 'INFO')

	logger.setLevel(logLevel)
	ch = logging.StreamHandler()
	ch.setLevel(logLevel)

	formatter = logging.Formatter('%(asctime)s - %(funcName)s - %(lineno)d - %(levelname)s - %(message)s', datefmt='%Y-%m-%d-%H%M%S')
	ch.setFormatter(formatter)
	logger.addHandler(ch)


	logger.info('sampleID: %s', args.sampleID)
	logger.info('sampleType: %s', args.sampleType)
	logger.info('hmmfile: %s', args.hmmfile)
	logger.info('outfile: %s', args.outfile)
	logger.info('cyclaseType: %s' % args.cyclaseType)
	logger.info('hmmscore_cutoff: %d' % args.hmmscore_cutoff)

	main(args.sampleID, args.sampleType, args.hmmfile, args.outfile, args.cyclaseType, args.hmmscore_cutoff)

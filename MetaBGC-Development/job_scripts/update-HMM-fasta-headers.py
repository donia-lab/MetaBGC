#!/usr/bin/python

"""
Created on Fri July 27 4:38:22 2018
@author: francinecamacho
"""
import os
from Bio import SeqIO
import pandas as pd 

def updateFastaHeader(fasta_dir, df, outdir ):
	
	for sampleFasta in os.listdir(fasta_dir):
		if sampleFasta.endswith(".fasta"):
			fasta_file = os.path.join(fasta_dir, sampleFasta)
			sampleTypeSeqIO = SeqIO.parse(open(fasta_file),'fasta')
			newFastaFile = sampleFasta.split(".fasta")[0] + "-updatedIDs.fasta" 

			with open(os.path.join(outdir, newFastaFile), "w") as output_handle:
				for record in sampleTypeSeqIO:
					recordID = record.id
					recordID_metaData = df[df.readID == recordID]
					if recordID_metaData.shape[0] > 1:
						row_max = recordID_metaData['HMM.Score'].idxmax()
						recordID_metaData = df.iloc[row_max].to_frame().T
						
					recordSample = recordID_metaData.Sample.to_string(index = False)
					recordHMMScore = recordID_metaData['HMM.Score'].to_string(index = False) 
					recordDomainType = recordID_metaData.cyclaseType.to_string(index = False)
					recordSampleType = recordID_metaData.sampleType.to_string(index = False)
					newRecordID = recordID + "__"+ recordSample +"__"+ recordSampleType + "__"+ recordHMMScore + "__"+ recordDomainType
					record.id = newRecordID
					record.description = ""
					SeqIO.write(record, output_handle, "fasta")


def main(fasta_dir, metadata_file, outdir):

	df = pd.read_csv(metadata_file,sep="\t")
	updateFastaHeader( fasta_dir, df, outdir)


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('--fasta_dir', type=str, required=True, help='directory where fasta files are located')
	parser.add_argument('--metadata_file', type=str, required=True, help='HMM results metadata file')
	parser.add_argument('--outdir', type=str,
						help='directory to output fasta file with metadata in header')


	args = parser.parse_args()

	main(args.fasta_dir, args.metadata_file, args.outdir) 

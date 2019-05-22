#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at fcamacho@princeton.edu).
#####################################################################################

import os
from Bio import SeqIO
import pandas as pd 

#Function to create a new fasta file with each record's (detected read) relevant metadata information.
def updateFastaHeader(fasta_file, metadata_df, out_file ):
	fasta_SeqIO = SeqIO.parse(open(fasta_file),'fasta')
	with open(out_file, "w") as output_handle:
		for record in fasta_SeqIO:
			record_id = record.id
			record_id_metadata = metadata_df[metadata_df.readID == record_id] 
			record_sample = record_id_metadata.Sample.to_string(index = False)
			record_sampleType = record_id_metadata.sampleType.to_string(index = False)  
			record_cyclase = record_id_metadata.cyclaseType.to_string(index = False)
			record_cohort = record_id_metadata.cohort.to_string(index = False)
			updated_record_id = record_id + "__"+ record_sample +"__"+ record_sampleType +"__"+ record_cyclase + "__"+record_cohort
			record.id = updated_record_id
			record.description = ""
			SeqIO.write(record, output_handle, "fasta")

def main(fasta_file, metadata_file, outfile):

	read_metadata = pd.read_csv(metadata_file,sep="\t")
	updateFastaHeader(fasta_file, read_metadata, outfile)
	print("Process Completed!")

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('-f','--fasta_file', type=str, required=True, help='fasta file of detected reads')
	parser.add_argument('-m', '--metadata_file', type=str, required=True, help='metadata file for reads')
	parser.add_argument('-o','--outfile', type=str, required=True,
						help='updated file with metadata in fasta file')


	args = parser.parse_args()

	main(args.fasta_file, args.metadata_file, args.outfile) 

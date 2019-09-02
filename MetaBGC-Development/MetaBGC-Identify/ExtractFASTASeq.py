#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

from Bio import SeqIO
import os
import re

"Function to parse fasta file based on text file with fasta header ids"
def parseFastaFile(fasta_file, output_file, id_list):
	fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
	with open(output_file, "a+") as f:
		for seq in fasta_sequences:
			if seq.id in id_list:
			   SeqIO.write(seq, f, "fasta")
	f.close()

def ExtractFASTASeq(fasta_dir,output_file,text_file):
	with open(text_file, "r") as ids: 
		id_list = ids.read() 
	ids.close()

	for subdir, dirs, files in os.walk(fasta_dir):
		for file in files:
			filePath = os.path.join(subdir, file)
			if re.match(r".*\.fasta$", file) and not re.match(r".*\.translated.fasta$", file) and os.path.getsize(
				filePath) > 0:
				parseFastaFile(filePath, output_file, id_list)

# if __name__ == '__main__':
# 	import argparse
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--fasta_file', required=True)
# 	parser.add_argument('--output_file', required=True)
# 	parser.add_argument('--text_file', required=True)
#
#
# 	args = parser.parse_args()
#
# 	main(args.fasta_file, args.output_file, args.text_file)

#!/usr/bin/env python

#####################################################################################
#@author: francinecamacho
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################
from metabgc.src.utils import *
from metabgc.src.blastrunlib import *



def combine_blast_results(blast_dir_path, combinedBLASTFile, cohort_name):
	filenames = [os.path.join(blast_dir_path, f) for f in os.listdir(blast_dir_path) if os.path.isfile(os.path.join(blast_dir_path, f)) and os.path.getsize(os.path.join(blast_dir_path, f)) > 0]
	tabular_colnames = "sseqid slen sstart send qseqid qlen qstart qend qcovs pident evalue"
	df_colnames = tabular_colnames.split()
	list_of_dfs = [pd.read_csv(filename, names=df_colnames, header=None, delim_whitespace=True) for filename in filenames]
	for dataframe, filename in zip(list_of_dfs, filenames):
		dataframe['Sample'] = os.path.basename(filename).split(".txt")[0]
		dataframe['cohort'] = cohort_name
	combined_df = pd.concat(list_of_dfs, ignore_index=True)
	combined_df.to_csv(combinedBLASTFile, index=False, sep='\t', header=False)
	numOfRows = combined_df.shape[0]
	return 	numOfRows

def create_clustering_file(blast_result,abundFile,abundWideFile):
	all_domains_blast_df = pd.read_csv(blast_result, sep="\t",
						   names=["sseqid", "slen","sstart", "send", "qseqid", "qlen", "qstart", "qend", "qcovs", "pident"," evalue", "Sample", "cohort"])
	all_domains_blast_df_count = all_domains_blast_df.groupby(['Sample','qseqid']).qseqid.agg('count').to_frame('count').reset_index()
	all_domains_blast_df_count_table = all_domains_blast_df_count.pivot_table(index='qseqid', columns='Sample', values='count',fill_value=0)
	all_domains_blast_df_count_table.to_csv(abundFile,sep='\t')
	all_domains_blast_df_count.to_csv(abundWideFile, index=False,sep='\t')

def mbgcquantify(identify_fasta, prot_family_name, cohort_name, nucl_seq_directory,
             seq_fmt, pair_fmt, r1_file_suffix, r2_file_suffix,blast_db_directory_map_file,
			 blastn_search_directory,output_directory, cpu):
	try:
		CPU_THREADS = 4
		if cpu is not None:
			CPU_THREADS = int(cpu)
		if r1_file_suffix is None:
			r1_file_suffix = ""
		if r2_file_suffix is None:
			r2_file_suffix = ""
		# Check if BLAST DB directory mapping file is provided or not
		if blast_db_directory_map_file is None:
			blast_db_directory_map_file = ""

		if blastn_search_directory is None:
			blastn_search_directory = os.path.join(output_directory, 'quantify_blastn_result')
		combinedBLASTFile = os.path.join(output_directory, "CombinedQuantifyBLAST.txt")
		abundFile = os.path.join(output_directory, "unique-biosynthetic-reads-abundance-table.txt")
		abundWideFile = os.path.join(output_directory, "unique-biosynthetic-reads-abundance-table-wide.txt")


		nucl_seq_directory = PreProcessReadsPar(nucl_seq_directory, seq_fmt, pair_fmt,
											 r1_file_suffix.strip(), r2_file_suffix.strip(),
											 output_directory, CPU_THREADS)

		cdHitFile = os.path.join(output_directory,"CombinedIDFASTASeqs_Drep.fasta")
		try:
			runCDHit(identify_fasta,cdHitFile,CPU_THREADS)
		except:
			print("Metabgc-quantify has failed during clustering of identified reads.")
			raise

		if not os.path.exists(blastn_search_directory):
			os.makedirs(blastn_search_directory, 0o777, True)
			RunPCMakeDBandBlastN(nucl_seq_directory, blast_db_directory_map_file,
								 cdHitFile, "blastn", "-dust no -max_target_seqs 1000000 -perc_identity 95.0 -qcov_hsp_perc 50 -window_size 11 -outfmt \"6 sseqid slen sstart send qseqid qlen qstart qend pident evalue\" ",
								 blastn_search_directory, CPU_THREADS)
		else:
			print("Metabgc-quantify is using the existing BLASTN hits found.")

		blastCount = combine_blast_results(blastn_search_directory, combinedBLASTFile, cohort_name)
		if blastCount == 0:
			print("Metabgc-quantify could not find any reads during quanify BLAST search.")
		else:
			print("Metabgc-quantify found " + str(blastCount) + " BLAST hits.")
		create_clustering_file(combinedBLASTFile, abundFile, abundWideFile)
		return abundFile, abundWideFile
	except:
		print("Metabgc-quantify has failed because no reads could be quantified. Please check your inputs and contact support on : https://github.com/donia-lab/MetaBGC")
		exit()
import os
from Bio import SeqIO
import sys
import shutil
import pandas as pd
from typing import TextIO, Tuple, Dict, Any

from Bio.SeqRecord import SeqRecord
import regex

def fuzzy_match(seq_a, seq_b, max_mismatches=2):
    shorter, longer = (seq_a, seq_b) if len(seq_a) <= len(seq_b) else (seq_b, seq_a)
    return bool(regex.search(f'({regex.escape(shorter)}){{s<={max_mismatches}}}', longer))

def exact_match(seq_a, seq_b):
    return seq_a == seq_b or seq_a in seq_b or seq_b in seq_a

def extract_genes_with_hit(genbank_file: str, list_of_hit_genes: list) -> [dict, dict]:
    # Initialize an empty list to store the genes within the range
    genes_within_range = {}
    dict_of_genes = {}
    # Parse the GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        record_name = record.name.replace('.', '_').split(sep='_')[0]
        dict_of_genes[record_name] = {}
        genes_within_range[record_name] = {}
        for feature in record.features:
            if feature.type == "CDS":
                locus_tag_id = feature.qualifiers["locus_tag"][0]
                dict_of_genes[record_name][locus_tag_id] = feature.qualifiers["product"][0]
                gene_start = int(feature.location.start)
                gene_end = int(feature.location.end)
                if "translation" in feature.qualifiers:
                    gene_sequence = feature.qualifiers["translation"][0]
                    # Check if the gene is in the hit list
                    for seq_record in list_of_hit_genes:
                        if fuzzy_match(str(seq_record.seq), gene_sequence):
                            genes_within_range[record_name][locus_tag_id] = [gene_start, gene_end]
                else:
                    print("No translation found for gene with locus tag: " + locus_tag_id)
                    #for seq_record in list_of_hit_genes:
                    #    if seq_record.seq == feature.location.extract(record).seq:
                    #        genes_within_range[sample_name][locus_tag_id] = [gene_start, gene_end]
                    
    return dict_of_genes, genes_within_range


def load_pgap_cds(file_name: str) -> tuple[dict[Any, Any], dict[Any, Any]]:
    cds_dict = {}
    id_locus_mapping_dict = {}
    for record in SeqIO.parse(file_name, "fasta"):
        desc_tokens = record.description.split(sep=' ')
        for desc_token in desc_tokens:
            if 'locus_tag' in desc_token:
                locus_token = desc_token.replace('[', '')
                locus_token = locus_token.replace(']', '')
                locus_tag = locus_token.split('=')[1]
                cds_dict[locus_tag] = record
                id_locus_mapping_dict[locus_tag] = record.id
                break
    return cds_dict, id_locus_mapping_dict


def load_pgap_faa(file_name: str) -> dict:
    faa_dict = {}
    for record in SeqIO.parse(file_name, "fasta"):
        locus_tag = record.id.split(sep='|')[2]
        faa_dict[locus_tag] = record
    return faa_dict

def split_fasta_file(fastaFile: str, seq_count: int, output_dir) -> None:
    write_rec = []
    output_file_list = []
    records = list(SeqIO.parse(fastaFile, "fasta"))
    out_prefix = os.path.splitext(os.path.basename(fastaFile))[0]
    file_count = 0
    for seq_rec in records:
        if len(write_rec) < seq_count:
            write_rec.append(seq_rec)
        else:
            output_file = os.path.join(output_dir, out_prefix + '-' + str(file_count) + ".faa")
            count = SeqIO.write(write_rec, output_file, "fasta")
            output_file_list.append(output_file)
            file_count = file_count + 1
            write_rec.clear()
    if len(write_rec) > 0:
        output_file = os.path.join(output_dir, out_prefix + '-' + str(file_count) + ".faa")
        count = SeqIO.write(write_rec, output_file, "fasta")
        output_file_list.append(output_file)
        file_count = file_count + 1
        write_rec.clear()

def print_hit_neighborhood(hit_gene_list: list, result_folder: str, extract_boundary: int,
                            sample_id: str, last_record_id: int, output_dir: str, f_out_product: TextIO,
                            f_out_id: TextIO, hit_seq_handle: TextIO) -> int:
    """
    Process blast results, extract neighborhood information, and write output records.

    Parameters:
    - hit_gene_list (list): List of hit genes.
    - result_folder (str): Path to the folder containing annotation files.
    - extract_boundary (int): Number of genes to extract in each direction from the hit gene.
    - bin_id (str): Identifier for the bin.
    - sample_id (str): Identifier for the sample.
    - last_record_id (int): Last record ID used.
    - output_dir (str): Directory to store output files.
    - f_out_product (TextIO): Output file handle for product records.
    - f_out_id (TextIO): Output file handle for ID records.
    - hit_seq_handle (TextIO): Output file handle for hit sequences.

    Returns:
    int: Updated last record ID.
    """

    print('Processing:{0}'.format(sample_id))

    gbk_file = os.path.join(result_folder, 'annot.gbk')
    cds_file = os.path.join(result_folder, 'annot_cds_from_genomic.fna')
    cds_dict, id_locus_mapping_dict = load_pgap_cds(cds_file)
    faa_file = os.path.join(result_folder, 'annot.faa')
    faa_dict = load_pgap_faa(faa_file)

    dict_of_genes, hit_genes = extract_genes_with_hit(gbk_file, hit_gene_list)
    # list_of_genes = list(dict_of_genes.keys())

    record_id = last_record_id
    central_hit_dict = {}
    hit_gene_found = False
    for scaff_name, genes in hit_genes.items():
        for gene in genes.keys():
            hit_start = genes[gene][0]
            hit_end = genes[gene][1]
            hit_gene_found = True
            if gene not in faa_dict:
                print("No faa record found for hit gene with locus tag: " + gene)
                continue
            list_of_genes = list(dict_of_genes[scaff_name].keys())
            gene_index = list_of_genes.index(gene)
            start_index = (gene_index - extract_boundary) if (gene_index - extract_boundary) > 0 else 0

            up_gene_neigh_list = list_of_genes[start_index:gene_index]
            gene_found_ctr = len(up_gene_neigh_list)
            up_gene_neigh_list = [''] * (extract_boundary-len(up_gene_neigh_list)) + up_gene_neigh_list

            end_index = (gene_index + extract_boundary) if (gene_index + extract_boundary) < len(list_of_genes) \
                else len(list_of_genes) - 1
            down_gene_neigh_list = list_of_genes[gene_index:end_index+1]
            gene_found_ctr = gene_found_ctr + len(down_gene_neigh_list)
            down_gene_neigh_list = down_gene_neigh_list + [''] * (extract_boundary + 1 - len(down_gene_neigh_list))

            gene_product_list = []
            gene_id_list = []
            faa_seq_list = []
            cds_seq_list = []
            for elem in up_gene_neigh_list:
                if elem in dict_of_genes[scaff_name]:
                    gene_product_list.append(dict_of_genes[scaff_name][elem])
                    gene_id_list.append(id_locus_mapping_dict[elem])
                    if elem in faa_dict:
                        faa_seq_list.append(faa_dict[elem])
                    cds_seq_list.append(cds_dict[elem])
                else:
                    gene_product_list.append('')
                    gene_id_list.append('')

            for elem in down_gene_neigh_list:
                if elem in dict_of_genes[scaff_name]:
                    gene_product_list.append(dict_of_genes[scaff_name][elem])
                    gene_id_list.append(id_locus_mapping_dict[elem])
                    if elem in faa_dict:
                        faa_seq_list.append(faa_dict[elem])
                    cds_seq_list.append(cds_dict[elem])
                else:
                    gene_product_list.append('')
                    gene_id_list.append('')

            gene_prod_str = '\t'.join(elem for elem in gene_product_list)
            gene_id_str = '\t'.join(elem for elem in gene_id_list)
            record_id = record_id + 1

            # Write product record
            f_out_product.write(sample_id + '\t' + str(record_id) +
                        '\t' + str(gene_found_ctr) + '\t' + str(len(genes)) + '\t' +
                                str(len(hit_genes)) + '\t' + scaff_name + '\t' + str(hit_start) + '\t' + str(hit_end) +
                        '\t' + gene_prod_str + '\n')

            # Write ID record
            f_out_id.write(sample_id + '\t' + str(record_id) + 
                        '\t' + str(gene_found_ctr) + '\t' + str(len(genes)) + '\t' +
                                str(len(hit_genes)) + '\t' + scaff_name + '\t' + str(hit_start) + '\t' + str(hit_end) +
                        '\t' + gene_id_str + '\n')

            # Save the central hit for annotating
            if gene in faa_dict:
                central_hit_dict[sample_id + '_' + str(record_id)] = faa_dict[gene]
            # Write the cds and fna files of the record
            with open(os.path.join(output_dir, sample_id + '_' + str(record_id) + '.faa'), "w") \
                    as output_handle:
                SeqIO.write(faa_seq_list, output_handle, "fasta")

            with open(os.path.join(output_dir, sample_id + '_' + str(record_id) + '.fna'), "w") \
                    as output_handle:
                SeqIO.write(cds_seq_list, output_handle, "fasta")

    # Write out the central hit sequences
    for k, v in central_hit_dict.items():
        v.description = v.id
        v.id = k
        SeqIO.write(v, hit_seq_handle, "fasta")
    return record_id, hit_gene_found

def write_mt_depth_summary(df_extracted_complete_ids, df_all_depth, output_dir):
    # Write depth summary file
    df_melted_all = pd.melt(df_extracted_complete_ids, id_vars=['BinID', 'SampleName', 'RecordID',
                                                            'CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample',
                                                            'MultiHit', 'ScaffoldName', 'HitStart', 'HitEnd'],
                            value_vars=columns_to_convert,
                            var_name='gene_position',
                            value_name='genes')
    df_melted_all = df_melted_all[df_melted_all['gene_position'] == '0']
    df_melted_all = df_melted_all.merge(df_all_depth, on=['genes'])
    df_melted_all.drop(columns=['CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample', 'MultiHit', 'HitStart',
                                'HitEnd', 'genes'], inplace=True)
    pivot_columns = df_melted_all.columns[5:]
    df_melted_all = pd.melt(df_melted_all, id_vars=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'gene_position'],
                       value_vars=pivot_columns, var_name='mt_samples', value_name='rpkm')
    df_melted_all['gene_position'] = df_melted_all['gene_position'].astype(int)
    df_melted_all = df_melted_all[df_melted_all['rpkm'] > 0]

    df_depth_summary = (df_melted_all.groupby(['BinID', 'SampleName', 'RecordID',
                                              'ScaffoldName', 'gene_position']).size().sort_values(ascending=False).
                        reset_index(name='sample_count'))
    df_depth_summary.to_csv(os.path.join(output_dir, 'mt_depth_summary.tsv'), sep='\t', index=False)

def write_mt_breath_summary(df_extracted_complete_ids, df_all_breath, output_dir):
    # Write breath summary file
    df_melted_all = pd.melt(df_extracted_complete_ids, id_vars=['BinID', 'SampleName', 'RecordID',
                                                            'CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample',
                                                            'MultiHit', 'ScaffoldName', 'HitStart', 'HitEnd'],
                            value_vars=columns_to_convert,
                            var_name='gene_position',
                            value_name='genes')
    df_melted_all = df_melted_all[df_melted_all['gene_position'] == '0']
    df_melted_all = df_melted_all.merge(df_all_breath, on=['genes'])
    df_melted_all.drop(columns=['CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample', 'MultiHit', 'HitStart',
                                'HitEnd', 'genes'], inplace=True)
    pivot_columns = df_melted_all.columns[5:]
    df_melted_all = pd.melt(df_melted_all, id_vars=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'gene_position'],
                       value_vars=pivot_columns, var_name='mt_samples', value_name='rpkm')
    df_melted_all['gene_position'] = df_melted_all['gene_position'].astype(int)
    df_melted_all = df_melted_all[df_melted_all['rpkm'] > 0]

    df_breath_summary = (df_melted_all.groupby(['BinID', 'SampleName', 'RecordID',
                                              'ScaffoldName', 'gene_position']).size().sort_values(ascending=False).
                        reset_index(name='sample_count'))
    df_breath_summary.to_csv(os.path.join(output_dir, 'mt_breath_summary.tsv'), sep='\t', index=False)


# if __name__ == '__main__':
#     record_id = 0
#     with open('test.tsv', 'w') as f_out, open('hit_seq.faa', 'w') as hit_seq_handle:
#         record_id = print_hit_neighborhood(r'C:\Users\ab50\Documents\data\DoniaLab\scaffold_gene_extract\87_612472597-scaffolds.txt',
#                                r'C:\Users\ab50\Documents\data\DoniaLab\scaffold_gene_extract\test\819_SRR413715-scaffolds_results',
#                                 20, '819', 'SRR413715-scaffolds', record_id, ".", f_out, hit_seq_handle)

if __name__ == '__main__':
    extract_boundary = 20
    pfam_name = 'paola_rumc'
    annot_dir  = '/scratch/gpfs/DONIA/abiswas/annotations/' + pfam_name  + '/annotations'
    blast_hit_genes  = '/scratch/gpfs/DONIA/abiswas/annotations/paola_rumc/sequences/uhgp_hits_rumc.fasta'
    output_dir = '/scratch/gpfs/DONIA/abiswas/annotations/' + pfam_name  + '/extract'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    mt_quantification = False
    mt_breath_file = 'Rum_quantified_breath.csv'
    mt_depth_file  = 'Rum_quantified_depth.csv'
    output_file = os.path.join(output_dir,'gene_neighbor_matrix_' + str(extract_boundary) + '.tsv')
    output_file_id = os.path.join(output_dir,'gene_neighbor_matrix_identifiers_' + str(extract_boundary) + '.tsv')
    output_file_blastp = os.path.join(output_dir,'blastp_gene_neighbor_matrix_' + str(extract_boundary) + '.tsv')
    output_seq_file = os.path.join(output_dir,'hits_genes.faa')
    blastp_file = os.path.join(output_dir,'hits_genes_blast.out')
    gene_files_dir = os.path.join(output_dir, "gene_files")
    if not os.path.exists(gene_files_dir):
        os.makedirs(gene_files_dir)
    record_id = 0
    hit_records_dict = {}
    for seq_id, seq_record in SeqIO.to_dict(SeqIO.parse(blast_hit_genes, "fasta")).items():
        seq_key = seq_id.replace('.', '_').split(sep='_')[0]
        if seq_key in hit_records_dict:
            hit_records_dict[seq_key].append(seq_record)
        else:
            hit_records_dict[seq_key] = [seq_record]

    annot_found_counter = 0
    hit_found_counter = 0
    hit_missing_list = []
    missing_annot_list = []
    with (open(output_file, 'w') as f_out_prod, open(output_file_id, 'w') as f_out_id,
          open(output_seq_file, 'w') as hit_seq_handle):
        f_out_prod.write('SampleName\tRecordID\tNumberOfGenes\tHitsInSample\tMultiHit\tScaffoldName\tHitStart\tHitEnd\t')
        f_out_id.write('SampleName\tRecordID\tNumberOfGenes\tHitsInSample\tMultiHit\tScaffoldName\tHitStart\tHitEnd\t')

        for i in range(-20, 21):
            f_out_prod.write(str(i) + '\t')
            f_out_id.write(str(i) + '\t')
        f_out_prod.write('\n')
        f_out_id.write('\n')

        for sample_annot in sorted(os.listdir(annot_dir)):
            sample_annot_tok = sample_annot.replace('.', '_').split(sep='_')
            sample_name = sample_annot_tok[0]
            gbk_file = os.path.join(annot_dir, sample_annot, 'annot.gbk')
            cds_file = os.path.join(annot_dir, sample_annot, 'annot_cds_from_genomic.fna')
            faa_file = os.path.join(annot_dir, sample_annot, 'annot.faa')
            if os.path.exists(gbk_file) and os.path.exists(cds_file) and os.path.exists(faa_file):
                annot_found_counter = annot_found_counter + 1
                print('Processing:' + sample_annot)
                record_id, hit_gene_found = print_hit_neighborhood(hit_records_dict[sample_name], os.path.join(annot_dir, sample_annot),
                                    20, sample_name, record_id, gene_files_dir,
                                                   f_out_prod, f_out_id, hit_seq_handle)
                if hit_gene_found:
                    hit_found_counter = hit_found_counter + 1
                else:
                    hit_missing_list.append(sample_annot)
            else:
                print('Annotation files missing for sample: ' + sample_annot)
                missing_annot_list.append(sample_annot)   
    # Generate splits to do blastp later
    split_out_dir = os.path.join(output_dir, "blast_splits")
    if not os.path.exists(split_out_dir):
        os.makedirs(split_out_dir)
    split_fasta_file(output_seq_file, 300, split_out_dir)

    print('Total annotations found: ' + str(annot_found_counter))
    print('Total hits found: ' + str(hit_found_counter))
    for missing in hit_missing_list:
        print('No hit gene found for sample: ' + missing)
    for missing in missing_annot_list:
        print('Annotation files missing for sample: ' + missing)

    if mt_quantification:
        # Load the extracted genes
        df_extracted_all = pd.read_csv(output_file, sep='\t')
        df_extracted_all = df_extracted_all.iloc[:, :-1]
        cols = ['BinID', 'SampleName', 'RecordID']
        df_extracted_all['combined_key'] = df_extracted_all[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        # Load the blastp output
        df_extracted_blastp = pd.read_csv(blastp_file, sep='\t', names=['qseqid', 'sseqid', 'pident', 'evalue',
                                                                        'staxids', 'sscinames', 'scomnames',
                                                                        'sskingdoms', 'stitle'])
        df_extracted_blastp.drop(['sseqid', 'pident', 'evalue',
                                  'staxids', 'scomnames', 'sskingdoms',
                                  'stitle'], axis=1, inplace=True)

        df_extracted_all = df_extracted_all.merge(df_extracted_blastp,
                                                  left_on='combined_key', right_on='qseqid', how='left')
        df_extracted_all.drop(['combined_key', 'qseqid'], axis=1, inplace=True)
        df_extracted_all.to_csv(output_file_blastp, sep='\t', index=False)

        # Pick the best central hit for bin,sample using the 'CentralGeneReadSupport'
        df_sorted = df_extracted_all.sort_values(['BinID', 'SampleName', 'MultiHit', 'CentralGeneReadSupport', 'NumberOfGenes'],
                                                 ascending=[True, True, True, False, False])
        # Use groupby to group by the 'BinID','SampleName' column and find the index of the first occurrence after sorting
        indices = df_sorted.groupby(['BinID', 'SampleName']).head(1).index
        df_extracted_complete = df_extracted_all.loc[indices]

        # Pick the best central hit for bin using the 'NumberOfGenes'
        df_sorted = df_extracted_complete.sort_values(['BinID', 'NumberOfGenes'], ascending=[True, False])
        # Use groupby to group by the 'BinID' column and find the index of the first occurrence after sorting
        indices = df_sorted.groupby('BinID').head(1).index
        df_extracted_complete = df_extracted_complete.loc[indices]

        complete_hit_records = [f"{binid}_{samplename}_{recordid}" for binid, samplename, recordid in
                          zip(df_extracted_complete['BinID'], df_extracted_complete['SampleName'],
                              df_extracted_complete['RecordID'])]

        # Generate list of complete hits
        complete_seq_file = os.path.join(output_dir,'max_sample_hits_genes.faa')
        with open(complete_seq_file, 'w') as hit_seq_handle:
            for record in SeqIO.parse(output_seq_file, "fasta"):
                if record.id in complete_hit_records:
                    SeqIO.write(record, hit_seq_handle, "fasta")
        complete_hits = 'max_sample_gene_neighbor_matrix_' + str(extract_boundary) + '.tsv'
        df_extracted_complete.to_csv(os.path.join(output_dir, complete_hits), sep='\t', index=False)

        # Copy the scaffold file
        scaff_files_dir = os.path.join(output_dir, "scaffold_files")
        if not os.path.exists(scaff_files_dir):
            os.makedirs(scaff_files_dir)
        for index, row in df_extracted_complete.iterrows():
            file_name = str(row['BinID']) + '_' + row['SampleName']
            src_scf_file = os.path.join(annot_dir, file_name + '_results', file_name + '.fasta')
            dest_scf_file = os.path.join(scaff_files_dir, file_name + '.fasta')
            shutil.copy(src_scf_file, dest_scf_file)

        # Copy the GBK file
        gbk_files_dir = os.path.join(output_dir, "gbk_files")
        if not os.path.exists(gbk_files_dir):
            os.makedirs(gbk_files_dir)
        for index, row in df_extracted_complete.iterrows():
            file_name = str(row['BinID']) + '_' + row['SampleName']
            src_gbk_file = os.path.join(annot_dir, file_name + '_results', 'annot.gbk')
            dest_gbk_file = os.path.join(gbk_files_dir, file_name + '.gbk')
            shutil.copy(src_gbk_file, dest_gbk_file)

        # Generate MT rpkm files
        mt_files_dir = os.path.join(output_dir, "mt_files")
        if not os.path.exists(mt_files_dir):
            os.makedirs(mt_files_dir)
        df_extracted_all_ids = pd.read_csv(output_file_id, sep='\t')
        df_extracted_all_ids = df_extracted_all_ids.iloc[:, :-1]
        df_extracted_complete_ids = df_extracted_all_ids[df_extracted_all_ids['RecordID'].isin(df_extracted_complete['RecordID'].tolist())]
        df_extracted_complete = df_extracted_complete.fillna('')
        if os.path.exists(mt_depth_file) and os.path.exists(mt_breath_file):
            df_all_depth = pd.read_csv(mt_depth_file, sep=',')
            df_all_breath = pd.read_csv(mt_breath_file, sep=',')
            columns_to_convert = [str(number) for number in list(range(-20, 21))]

            write_mt_depth_summary(df_extracted_complete_ids, df_all_depth, output_dir)

            write_mt_breath_summary(df_extracted_complete_ids, df_all_breath, output_dir)

            for index, selected_row in df_extracted_complete.iterrows():

                row_df = df_extracted_all_ids[(df_extracted_all_ids['BinID'] == selected_row['BinID']) &
                                                    (df_extracted_all_ids['SampleName'] == selected_row['SampleName']) &
                                                    (df_extracted_all_ids['RecordID'] == selected_row['RecordID'])]
                df_melted = pd.melt(row_df, id_vars=['BinID', 'SampleName', 'RecordID',
                                                     'CentralGeneReadSupport', 'NumberOfGenes',
                                                     'HitsInSample','MultiHit', 'ScaffoldName',
                                                     'HitStart', 'HitEnd'],
                                    value_vars=columns_to_convert,
                                    var_name='gene_position',
                                    value_name='genes')

                # Write the depth files
                df_depth = df_melted.merge(df_all_depth, on=['genes'], how='left')
                df_depth = df_depth.fillna(0)
                df_depth.drop(columns=['CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample', 'MultiHit', 'HitStart',
                                       'HitEnd', 'genes'], inplace=True)

                pivot_columns = df_depth.columns[5:]
                df_depth = pd.melt(df_depth, id_vars=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'gene_position'],
                                   value_vars=pivot_columns, var_name='mt_samples', value_name='rpkm')
                df_depth['gene_position'] = df_depth['gene_position'].astype(int)
                # Sort the DataFrame based on the column
                df_depth = df_depth.pivot(index=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'mt_samples'],
                                          columns='gene_position', values='rpkm').reset_index()
                depth_file_name = (str(selected_row['BinID']) + '_' + selected_row['SampleName'] + '_'
                                   + str(selected_row['RecordID']) + '_depth.tsv')
                with open(os.path.join(mt_files_dir, depth_file_name), 'w') as f_out:
                    # Write the column names
                    header_row = '\t'.join(map(str, df_depth.columns))
                    f_out.write(header_row)
                    f_out.write('\n')
                    # Write the gene names
                    f_out.write('\t\t\t\t\t')
                    f_out.write('\t'.join(map(str, selected_row[10:51])))
                    f_out.write('\n')
                    # Write the summary row
                    f_out.write('\t\t\t\t\t')
                    for col in df_depth.columns[5:46]:
                        zero_count = (df_depth[col] != 0).sum()
                        f_out.write(str(zero_count) + '\t')
                    f_out.write('\n')
                    for index, row in df_depth.iterrows():
                        f_out.write('\t'.join(map(str, row)))
                        f_out.write('\n')



                # Write the breath files
                df_breath = df_melted.merge(df_all_breath, on=['genes'], how='left')
                df_breath = df_breath.fillna(0)
                df_breath.drop(columns=['CentralGeneReadSupport', 'NumberOfGenes', 'HitsInSample', 'MultiHit', 'HitStart',
                                       'HitEnd', 'genes'], inplace=True)
                pivot_columns = df_breath.columns[5:]
                df_breath = pd.melt(df_breath, id_vars=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'gene_position'],
                                   value_vars=pivot_columns, var_name='mt_samples', value_name='rpkm')
                df_breath['gene_position'] = df_breath['gene_position'].astype(int)
                # Sort the DataFrame based on the column
                df_breath = df_breath.pivot(index=['BinID', 'SampleName', 'RecordID', 'ScaffoldName', 'mt_samples'],
                                          columns='gene_position', values='rpkm').reset_index()
                breath_file_name = (str(selected_row['BinID']) + '_' + selected_row['SampleName'] + '_'
                                    + str(selected_row['RecordID']) + '_breath.tsv')
                with open(os.path.join(mt_files_dir, breath_file_name), 'w') as f_out:
                    # Write the column names
                    header_row = '\t'.join(map(str, df_breath.columns))
                    f_out.write(header_row)
                    f_out.write('\n')
                    # Write the gene names
                    f_out.write('\t\t\t\t\t')
                    f_out.write('\t'.join(map(str, selected_row[10:51])))
                    f_out.write('\n')
                    # Write the summary row
                    f_out.write('\t\t\t\t\t')
                    for col in df_breath.columns[5:46]:
                        zero_count = (df_breath[col] != 0).sum()
                        f_out.write(str(zero_count) + '\t')
                    f_out.write('\n')
                    for index, row in df_breath.iterrows():
                        f_out.write('\t'.join(map(str, row)))
                        f_out.write('\n')



import os
import sys

if len(sys.argv) > 3:
    hmm_dir = sys.argv[1]
    output_dir = sys.argv[2]
    process_count = int(sys.argv[3])
    os.makedirs(output_dir, exist_ok=True)
    print(f"Generating all the HMMER commands...")
    fastq_dir = "/scratch/gpfs/DONIA/abiswas/MetaBGCRuns/DATAFILES/search/ALL/prot"
    # Loop over all HMM models
    with open('hmmsearch_cmds.txt', 'w') as f_out:
        for hmm_filename in os.listdir(hmm_dir):
            if hmm_filename.endswith(".hmm"):
                hmm_filepath = os.path.join(hmm_dir, hmm_filename)
                # Loop over all fasta files
                for fa_filename in os.listdir(fastq_dir):
                    fa_filepath = os.path.join(fastq_dir, fa_filename)
                    interval_str = os.path.splitext(hmm_filename)[0]
                    interval_str = interval_str.split('__')[-1]
                    output_filename = os.path.splitext(fa_filename)[0] + '__' + interval_str + '.tbl'
                    output_filepath = os.path.join(output_dir, output_filename)
                    if not os.path.isfile(output_filepath):
                        f_out.write('hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ' + output_filepath + ' ' + hmm_filepath + ' ' + fa_filepath + ' > /dev/null\n')
    # Create the splits
    os.makedirs("hmmersearch_splits", exist_ok=True)
    with open('hmmsearch_cmds.txt', 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    chunk_size = len(lines) // process_count + (len(lines) % process_count > 0)
    for i in range(process_count):
        chunk_lines = lines[i * chunk_size:(i + 1) * chunk_size]
        if not chunk_lines:
            break
        with open(f"hmmersearch_splits/split_{i + 1:04d}.txt", "w") as out:
            out.write("\n".join(chunk_lines) + "\n")
else:
   print("Not all parameters provided. Three parameters are required. The directory with " + \
          "spHMM models, the output directory of the search results, and the number od splits.")

from Bio import SearchIO
import os
import pandas as pd
import sys

class HMMRecord:
	def __init__(self, acc_numb, sampleType, sampleID, protType, bitscore, window, interval):
		self.acc_numb = acc_numb
		self.sampleType = sampleType
		self.sampleID = sampleID
		self.protType = protType
		self.bitscore = bitscore
		self.window = window
		self.interval = interval


"""
Function to parse HMM file into HMMRecord dict.
"""
def parseHMM(hmmPathFile, hmm_string_fmt, sampleType, sampleID, protType, window, interval):
    with open(hmmPathFile, 'r') as handle:
        results_dict = {}
        try:
            for record in SearchIO.parse(handle, hmm_string_fmt):
                hits = record.hits
                num_hits = len(hits)  #calculate how many hits per query
                if num_hits > 0:  #extract hits data
                    for i in range(0, num_hits):
                        hmm_name = hits[i].id
                        hmm_bitscore = hits[i].bitscore
                        if hmm_name not in results_dict: # add hits to results dictionary
                            hmmRec = HMMRecord(hmm_name, sampleType, sampleID, protType, hmm_bitscore, window, interval)
                            results_dict[hmm_name] = hmmRec
            handle.close()
        except ValueError:
            print ("ERROR: duplicated queryIDs in hmmer result file:", hmmPathFile)
    return results_dict


def createPandaDF(hmm_dict, outfile):
    outputDF = pd.DataFrame()
    if len(hmm_dict) != 0:  # check if counter dict is not empty to add to df
        df_rows = [(k, v.sampleType, v.sampleID, v.protType, v.bitscore, v.window, v.interval) for k, v in list(hmm_dict.items())]  # convert dictionary to list
        outputDF = pd.DataFrame(df_rows, columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window","interval"])
        sorted_outputDF = outputDF.sort_values(by=['HMMScore'],ascending=[False])  # sort in decending order by Hit.counts column
        sorted_outputDF.to_csv(outfile, index=False, sep='\t', header=False)
    else:
        outputDF_empty_columns = ["readID", "sampleType", "sampleID", "protType", "HMMScore", "window", "interval"]  # rename column names
        outputDF_empty = pd.DataFrame(columns=outputDF_empty_columns)
        outputDF_empty.to_csv(outfile, index=False, sep='\t', header=False)



if __name__ == '__main__':
    if len(sys.argv) > 3:
        protType = sys.argv[1]
        hmm_output_dir = sys.argv[2]
        identify_directory = sys.argv[3]
        os.makedirs(identify_directory, exist_ok=True)
        print(f"Merging all the HMMER output...")
        window = "30_10"
        process_ctr = 0
        for tbl_filename in os.listdir(hmm_output_dir):
            tbl_filepath = os.path.join(hmm_output_dir, tbl_filename)
            output_filepath = os.path.join(hmm_output_dir, os.path.splitext(tbl_filename)[0] + '.txt')
            if tbl_filename.endswith(".tbl") and not os.path.exists(output_filepath):
                sampleID = os.path.splitext(tbl_filename)[0].split('__')[0]
                interval = os.path.splitext(tbl_filename)[0].split('__')[-1]
                hmm_dict = parseHMM(tbl_filepath, "hmmer3-tab", 'ALL', sampleID, protType, window, interval)
                createPandaDF(hmm_dict, output_filepath)
            process_ctr = process_ctr + 1
            if process_ctr % 100 == 0:
                print('Processed: ' + str(process_ctr))

        allHMMResult = os.path.join(identify_directory, "CombinedHmmSearch.txt")
        found_hit_ctr = 0
        with open(allHMMResult, 'w') as outfile:
            for subdir, dirs, files in os.walk(hmm_output_dir):
                for file in files:
                    filePath = os.path.join(subdir, file)
                    if file.endswith(".txt") and os.path.getsize(filePath) > 0:
                        with open(filePath) as infile:
                            for line in infile:
                                outfile.write(line)
                                found_hit_ctr = found_hit_ctr + 1
        print('Found:' + str(found_hit_ctr))
    else:
        print("Not all parameters provided. Three parameters are required. The protein name, the directory with " + \
              "HMMER search results and the identify output directory.")

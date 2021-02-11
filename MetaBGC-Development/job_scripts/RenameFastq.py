from Bio import SeqIO
import sys

if __name__ == '__main__':
    fastq_file = sys.argv[1]
    name_suffix = sys.argv[2]
    output_file = sys.argv[3]
    records = list(SeqIO.parse(fastq_file, "fastq"))
    process_ctr = 0
    no_of_reads = len(records)
    for record in records:
        # read_name = record.id + record.description
        read_name = record.description
        read_name = read_name.split("length=")[0]
        read_name = read_name.strip().replace(" ","_")
        read_name = read_name + name_suffix
        record.id = read_name
        record.description = ""
        process_ctr = process_ctr + 1
        if (process_ctr % 50000) == 0:
            print("Processed " + str(process_ctr) + " of " + str(no_of_reads) + "...")
    count = SeqIO.write(records, output_file, "fastq")
    print("Processing complete. " + str(process_ctr) + " of " + str(no_of_reads) + " processed.")
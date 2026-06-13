import os
from Bio import SeqIO
import sys

if __name__ == '__main__':
    sample_name = sys.argv[1]
    scaffold_file = sys.argv[2]
    bin_dir = sys.argv[3]

    # Loop and load the binned scaffolds
    result = [os.path.join(dp, f)
              for dp, dn, filenames in os.walk(bin_dir)
              for f in filenames if (os.path.splitext(f)[1] == '.fna' or
                                     os.path.splitext(f)[1] == '.fa')]
    if len(result) > 0 and os.path.isfile(scaffold_file):
        binned_scaffolds = []
        for bin_file in result:
            for record in SeqIO.parse(bin_file, "fasta"):
                binned_scaffolds.append(record.id)

        # load thr original scaffolds
        unbinned_records = []
        for record in SeqIO.parse(scaffold_file, "fasta"):
            if record.id not in binned_scaffolds:
                unbinned_records.append(record)

        ext_str = os.path.splitext(result[0])[1]
        output_file = os.path.join(bin_dir, sample_name + '.unbinned' + ext_str)
        SeqIO.write(unbinned_records, output_file, "fasta")

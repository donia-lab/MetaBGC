import ftputil
import os, sys, os.path
from pathlib import Path

if __name__ == '__main__':
    acc_list_file = sys.argv[1]
    local_base_dir = sys.argv[2]
    host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
    with open(acc_list_file) as fp:
        for line in fp:
            if line.strip():
                acc_id = Path(line.strip()).resolve().stem
                asm_type = acc_id.split('_')[0]
                acc_code = acc_id.split('_')[1]
                asm_code = "_".join(acc_id.split('_')[2:])
                dirpath = '/genomes/all/' + asm_type + '/' + acc_code[0:3] + '/' + acc_code[3:6] + '/' + acc_code[6:9] + '/' + acc_id
                host.chdir(dirpath)
                filenames = host.listdir(dirpath)
                local_dirname = os.path.join(local_base_dir, acc_id)
                os.makedirs(local_dirname, 0o777, True)
                for filename in filenames:
                    if host.path.isfile(filename):
                        local_filename = os.path.join(local_dirname, filename)
                        #print("Downloading file " + os.path.join(local_dirname, filename))
                        host.download(filename, os.path.join(local_dirname, filename))
                print(line.strip() + " completed...")
    host.close()

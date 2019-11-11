#!/bin/bash

#SBATCH -N 1 # node count
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=30G
#SBATCH --qos=short
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited
export PATH="/tigress/DONIA/data/donia/abiswas/tools/hmmer-3.1b1/src:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/ncbi-blast-2.7.1+/bin:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/muscle3.8:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/cd-hit-v4.8.1-2019-0228:$PATH"
module load anaconda3
conda activate /tigress/DONIA/data/donia/abiswas/tools/py37 

DATA_PATH=/tigress/DONIA/abiswas/git/MetaBGC-V1.2/data/OxyN
NUCL_READ_PATH=${DATA_PATH}/reads
INSTALL_PATH=/tigress/DONIA/abiswas/git/MetaBGC-V1.2

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")

export PYTHONPATH="${INSTALL_PATH}"

python ${INSTALL_PATH}/metabgc.py identify --sphmm_directory ${DATA_PATH}/output/build/HiPer_spHMMs --prot_family_name Cyclase_OxyN --cohort_name OxyN --nucl_seq_directory ${DATA_PATH}/output/build/nucl_seq_dir --seq_fmt FASTA --pair_fmt interleaved --output_directory ${DATA_PATH}/output --cpu 20

python ${INSTALL_PATH}/metabgc.py quantify --identify_fasta ${DATA_PATH}/output/identified-biosynthetic-reads.fasta --prot_family_name Cyclase_OxyN --cohort_name OxyN --nucl_seq_directory ${DATA_PATH}/output/build/nucl_seq_dir --seq_fmt FASTA --pair_fmt interleaved --output_directory ${DATA_PATH}/output --cpu 1

python ${INSTALL_PATH}/metabgc.py cluster --table ${DATA_PATH}/output/unique-biosynthetic-reads-abundance-table.txt

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.



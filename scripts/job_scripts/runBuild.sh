#!/bin/bash

#SBATCH -N 1 # node count
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --mem=30G
#SBATCH --qos=short
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited
export PATH="/tigress/DONIA/data/donia/abiswas/tools/hmmer-3.1b1/src:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/ncbi-blast-2.7.1+/bin:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/muscle3.8:$PATH"
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

#rm -rf ${DATA_PATH}/output
#mkdir ${DATA_PATH}/output
python ${INSTALL_PATH}/metabgc.py build --prot_alignment ${DATA_PATH}/OxyN.fasta --prot_family_name Cyclase_OxyN --cohort_name OxyN --nucl_seq_directory ${NUCL_READ_PATH} --seq_fmt FASTQ --pair_fmt split --r1_file_suffix 1.fastq --r2_file_suffix 2.fastq --tp_genes_nucl ${DATA_PATH}/OxyN-True-Positive-genes.fasta --f1_thresh 0.25 --output_directory ${DATA_PATH}/output --cpu 20

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


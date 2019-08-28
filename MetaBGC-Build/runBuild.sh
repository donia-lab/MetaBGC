#!/bin/bash

#SBATCH -N 1 # node count
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH --mem=50G
#SBATCH --qos=short
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited
export PATH="/tigress/DONIA/data/donia/abiswas/tools/hmmer-3.1b1/src:$PATH"
export PATH="/tigress/DONIA/data/donia/abiswas/tools/ncbi-blast-2.7.1+/bin:$PATH"

module load anaconda3
conda activate /tigress/DONIA/data/donia/abiswas/tools/py37 

DATA_PATH=/tigress/DONIA/abiswas/git/MetaBGC-TIIPKS/benchmark_data/OxyN
NUCL_READ_PATH=/projects/DONIA/abiswas/git/MetaBGC-TIIPKS/benchmark_data/OxyN/reads
PROT_READ_PATH=/projects/DONIA/abiswas/git/MetaBGC-TIIPKS/benchmark_data/OxyN/reads
INSTALL_PATH=/tigress/DONIA/abiswas/git/MetaBGC-TIIPKS

cd ${INSTALL_PATH}

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")

export PYTHONPATH="${INSTALL_PATH}"

rm -rf ${DATA_PATH}/output

python ${INSTALL_PATH}/MetaBGC-Build/MetaBGC-Build.py --prot_alignment ${DATA_PATH}/OxyN.fasta --prot_family_name OxyN --cohort_name test --nucl_seq_directory ${NUCL_READ_PATH} --prot_seq_directory ${PROT_READ_PATH} --tp_genes ${DATA_PATH}/TPGenes/T2PKS-positive-genomes-cyclases-1-21-2019.fasta --gene_pos ${DATA_PATH}/Gene_Interval_Pos.txt --F1_Thresh 0.5 --output_directory ${DATA_PATH}/output --cpu 32

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


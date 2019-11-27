#!/bin/bash

#SBATCH -N 1 # node count
#SBATCH -A molbio
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

DATA_PATH=/projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs
NUCL_READ_PATH=${DATA_PATH}/reads
PROT_READ_PATH=${DATA_PATH}/reads_prot
INSTALL_PATH=/projects/DONIA/abiswas/git/MetaBGC/MetaBGC-Development

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")

export PYTHONPATH="${INSTALL_PATH}"

#rm -rf ${DATA_PATH}/output
#mkdir ${DATA_PATH}/output
python ${INSTALL_PATH}/runmetabgc.py build --prot_alignment ${DATA_PATH}/AcbK-homologs-for-HMM.fasta --prot_family_name Transport --cohort_name AbcK --nucl_seq_directory ${NUCL_READ_PATH} --prot_seq_directory ${PROT_READ_PATH} --seq_fmt FASTA --pair_fmt interleaved --tp_genes_nucl ${DATA_PATH}/TP_Genes.fasta --blastn_search_directory ${DATA_PATH}/output/build/blastn_result --f1_thresh 0.5 --output_directory ${DATA_PATH}/output --cpu 2

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


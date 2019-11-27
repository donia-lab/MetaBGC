#!/bin/bash

#SBATCH -N 1 # node count
#SBATCH -A molbio
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --mem=10G
#SBATCH --qos=short
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited
export PATH="/tigress/DONIA/data/donia/abiswas/tools/hmmer-3.1b1/src:$PATH"
OP_PATH=/projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/hmm_search/results
HMM_PATH=/projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/output/build/spHMMs
cd /projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/hmm_search/results

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")
echo Number of cores assigned is "$SLURM_CPUS_ON_NODE"
IFS=':' read -r -a FILE_NAMES <<< "$FILES"
FILENAME=${FILE_NAMES[$SLURM_ARRAY_TASK_ID]}
# make & move into new directory, and run!

f="$(basename -- $FILENAME)";
f=${f%.fasta};
echo "Processing: ${FILENAME}";

hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_0_30.tbl ${HMM_PATH}/Transport__30_10__0_30.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_100_130.tbl ${HMM_PATH}/Transport__30_10__100_130.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_10_40.tbl ${HMM_PATH}/Transport__30_10__10_40.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_110_140.tbl ${HMM_PATH}/Transport__30_10__110_140.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_120_150.tbl ${HMM_PATH}/Transport__30_10__120_150.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_130_160.tbl ${HMM_PATH}/Transport__30_10__130_160.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_140_170.tbl ${HMM_PATH}/Transport__30_10__140_170.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_150_180.tbl ${HMM_PATH}/Transport__30_10__150_180.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_160_190.tbl ${HMM_PATH}/Transport__30_10__160_190.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_170_200.tbl ${HMM_PATH}/Transport__30_10__170_200.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_180_210.tbl ${HMM_PATH}/Transport__30_10__180_210.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_190_220.tbl ${HMM_PATH}/Transport__30_10__190_220.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_200_230.tbl ${HMM_PATH}/Transport__30_10__200_230.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_20_50.tbl ${HMM_PATH}/Transport__30_10__20_50.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_210_240.tbl ${HMM_PATH}/Transport__30_10__210_240.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_220_250.tbl ${HMM_PATH}/Transport__30_10__220_250.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_230_260.tbl ${HMM_PATH}/Transport__30_10__230_260.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_240_270.tbl ${HMM_PATH}/Transport__30_10__240_270.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_250_280.tbl ${HMM_PATH}/Transport__30_10__250_280.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_260_290.tbl ${HMM_PATH}/Transport__30_10__260_290.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_270_300.tbl ${HMM_PATH}/Transport__30_10__270_300.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_280_310.tbl ${HMM_PATH}/Transport__30_10__280_310.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_290_320.tbl ${HMM_PATH}/Transport__30_10__290_320.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_30_60.tbl ${HMM_PATH}/Transport__30_10__30_60.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_40_70.tbl ${HMM_PATH}/Transport__30_10__40_70.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_50_80.tbl ${HMM_PATH}/Transport__30_10__50_80.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_60_90.tbl ${HMM_PATH}/Transport__30_10__60_90.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_70_100.tbl ${HMM_PATH}/Transport__30_10__70_100.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_80_110.tbl ${HMM_PATH}/Transport__30_10__80_110.hmm ${FILENAME} > /dev/null
hmmsearch --cpu 1 --F1 0.02 --F2 0.02 --F3 0.02 --tblout ${OP_PATH}/${f}_90_120.tbl ${HMM_PATH}/Transport__30_10__90_120.hmm ${FILENAME} > /dev/null

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


#!/bin/bash

#SBATCH -o PGAP_RUN
#SBATCH -A donia
###SBATCH -p donia
###SBATCH --nodelist della-h13n1
#SBATCH -N 1 # node count
#SBATCH -c 4
#SBATCH -t 48:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")
echo Number of cores assigned is "$SLURM_CPUS_ON_NODE"
IFS=':' read -r -a SAMPLE_LIST <<< "$ALL_SAMPLE_ARRAY"
SAMPLE_STRING=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}
echo "Running Command: ${SAMPLE_STRING}";

PFAM_NAME=AAC6__III-homologs
RUN_DIR=/scratch/gpfs/DONIA/abiswas/annotations/${PFAM_NAME}/run_folders/${SAMPLE_STRING}
LOG_FILE=/scratch/gpfs/DONIA/abiswas/annotations/${PFAM_NAME}/scripts/logs/${SAMPLE_STRING}.log
OUT_DIR=/scratch/gpfs/DONIA/abiswas/annotations/${PFAM_NAME}/annotations/${SAMPLE_STRING}_results

if [ ! -d "$OUT_DIR" ]; then

	/scratch/gpfs/DONIA/abiswas/tools/pgap_latest_old/pgap.py -n --no-internet --ignore-all-errors --no-self-update --container-path /scratch/gpfs/DONIA/abiswas/tools/pgap/pgap_2022-12-13.build6494.sif -o ${OUT_DIR} ${RUN_DIR}/input.yaml &> ${LOG_FILE}

fi 

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


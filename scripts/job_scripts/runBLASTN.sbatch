#!/bin/bash

#SBATCH -o BIN_SCF_BLAST
#SBATCH -a molbio
###SBATCH -p donia
###SBATCH --nodelist della-h13n1
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --mem=5G
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

ulimit -s unlimited
module load blast/2.10.0
export PATH="/tigress/DONIA/data/donia/abiswas/tools/bbmap:$PATH"
echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")
echo Number of cores assigned is "$SLURM_CPUS_ON_NODE"
IFS=':' read -r -a FILE_NAMES <<< "$FILES"
FILENAME=${FILE_NAMES[$SLURM_ARRAY_TASK_ID]}
echo "Searching DB: ${FILENAME}";
dbf="$(basename -- $FILENAME)";
dbf=${dbf%.fasta}

IFS=':' read -r -a BIN_NAMES <<< "$BINS"
BINNAME=${BIN_NAMES[$SLURM_ARRAY_TASK_ID]}
echo "Searching Query: ${BINNAME}";
f="$(basename -- $BINNAME)";
f=${f%.fasta};

OP_PATH=scf_blast_results
mkdir $OP_PATH
mkdir $OP_PATH/scf_match

makeblastdb -in ${FILENAME} -title ${dbf} -dbtype nucl -parse_seqids -out ${FILENAME}

#IFS='_' read -r -a array <<< "$f"

blastn -num_threads 1 -query ${BINNAME} -db $FILENAME -dust no -max_target_seqs 1 -perc_identity 95.0 -qcov_hsp_perc 50 -window_size 11 -outfmt "6 sseqid slen sstart send qseqid qlen qstart qend pident evalue" -out ${OP_PATH}/${f}_${dbf}.txt

cut -f1 ${OP_PATH}/${f}_${dbf}.txt > ${OP_PATH}/${f}_${dbf}.IDs.txt

filterbyname.sh in=$FILENAME out=${OP_PATH}/scf_match/${f}_${dbf}.fasta names=${OP_PATH}/${f}_${dbf}.IDs.txt include=t

#rm -rf ${OP_PATH}/${f}

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


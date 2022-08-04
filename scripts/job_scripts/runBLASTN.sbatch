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
export PATH="/tigress/DONIA/data/donia/abiswas/tools/ncbi-blast-2.9.0+/bin:$PATH"
Q_FILE=/projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/TP_Genes.fasta 
OP_PATH=/projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/search/results
cd ${OP_PATH}

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

#mkdir ${OP_PATH}/${f}
#makeblastdb -in ${FILENAME} -title ${f} -dbtype nucl -parse_seqids -out ${OP_PATH}/${f}/${f}

IFS='_' read -r -a array <<< "$f"

blastn -num_threads 1 -query ${Q_FILE} -db /projects/DONIA/data/synthetic_new/${array[2]}/${f}/${f}-raw-reads-fasta/${f}.fasta -perc_identity 90.0 -max_target_seqs 10000 -outfmt "6 sseqid slen sstart send qseqid qlen qstart qend pident evalue" -out ${OP_PATH}/${f}.txt

#rm -rf ${OP_PATH}/${f}

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.


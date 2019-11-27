#!/bin/bash

#SBATCH -p donia
#SBATCH --nodelist tiger-h19d1
#SBATCH -N 1 # node count
#SBATCH -c 10
#SBATCH -t 90:00:00
#SBATCH --mem=100000
#SBATCH --mail-type=end
#SBATCH --mail-user=fcamacho@princeton.edu
#SBATCH -D /tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/T2PKS-TP-reads

 
R CMD BATCH combine_BLAST-2PKS-results.R

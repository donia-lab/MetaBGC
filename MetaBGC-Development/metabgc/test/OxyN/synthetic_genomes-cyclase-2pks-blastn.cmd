#!/bin/bash

#SBATCH -p donia
#SBATCH --nodelist tiger-h19d1
#SBATCH -N 1 # node count
#SBATCH -c 56
#SBATCH -t 360:00:00
#SBATCH --mem=2000000
#SBATCH --mail-type=end
#SBATCH --mail-user=fcamacho@princeton.edu
#SBATCH -D /tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/T2PKS-TP-reads

module load blast
mgprun3.1 -f synthetic_genomes-cyclase-2pks-blastn.yaml --dbmode none

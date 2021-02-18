#!/bin/bash

#SBATCH -p donia
#SBATCH --nodelist=tiger-h19d1
#SBATCH -N 1 # node count
#SBATCH -c 10
#SBATCH -t 71:59:59
#SBATCH --mem=200000
#SBATCH --mail-type=end
#SBATCH --mail-user=shuow@princeton.edu

python synthesize_metagenomic_samples_PKS_random.py --indir1 /tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/background_genomes --indir2 /tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/T2PKS-Positive-genomes/T2PKS_pos-genomes-synthetic-metagenome -o /tigress/DONIA/data/synthetic_new/high -ss HS20 -l 100 -fl 400 -sd 20 -nr 51100000 -p 0.9 -ns 70 -t 10 -rs 412 -b synthetic_stool_high_S


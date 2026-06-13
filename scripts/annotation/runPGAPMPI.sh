#!/bin/bash
#SBATCH -J pgap_mpi
#SBATCH -o pgap_mpi_%j.out
#SBATCH -e pgap_mpi_%j.err
#SBATCH -A donia
#SBATCH -N 4                    # 1 controller rank + 3 worker ranks (scale -N as needed)
#SBATCH --ntasks-per-node 1     # 1 MPI rank per node; PGAP is internally multithreaded
#SBATCH -c 8                    # CPUs per task (PGAP uses ~4-8 threads)
#SBATCH --mem=20G               # Memory per node (PGAP default ~10G; 20G gives headroom)
#SBATCH -t 72:00:00
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=ab50@princeton.edu

# ---- Environment setup ----
module purge
module load openmpi               # or intel-mpi, depending on cluster config
module load singularity           # required for PGAP container
module load anaconda3

conda activate pgap_env           # environment with mpi4py and biopython installed

ulimit -s unlimited

# ---- Paths ----
FASTA=/scratch/gpfs/DONIA/abiswas/annotations/input_scaffolds.fasta
OUTDIR=/scratch/gpfs/DONIA/abiswas/annotations/AAC6__III-homologs/annotations
SUBMOL=/scratch/gpfs/DONIA/abiswas/annotations/AAC6__III-homologs/scripts/submol.yaml
PGAP_EXE=/scratch/gpfs/DONIA/abiswas/tools/pgap_latest_old/pgap.py
CONTAINER=/scratch/gpfs/DONIA/abiswas/tools/pgap/pgap_2022-12-13.build6494.sif

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---- Run ----
srun python "${SCRIPT_DIR}/run_pgap_mpi.py" \
    --fasta      "${FASTA}"     \
    --outdir     "${OUTDIR}"    \
    --submol     "${SUBMOL}"    \
    --pgap-exe   "${PGAP_EXE}"  \
    --container  "${CONTAINER}"

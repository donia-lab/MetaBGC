#!/bin/bash

#SBATCH --array=0-1012
#SBATCH --job-name=hmmer_search_array
#SBATCH --output=logs/hmmer_search_%A_%a.log
#SBATCH --error=logs/hmmer_search_error_%A_%a.log
#SBATCH -t 72:00:00
#SBATCH -N 1 # node count
#SBATCH -n 1
#SBATCH --mem=100M
#SBATCH --mail-type=end
#SBATCH --mail-user=ab50@princeton.edu

export PATH="/tigress/DONIA/data/donia/abiswas/tools/hmmer-3.1b1/src:$PATH"

# Create logs directory if not exist
mkdir -p logs

# Directory with command files
CMD_DIR="hmmersearch_splits"

# Get the file corresponding to this array task
CMD_FILE=$(ls ${CMD_DIR}/* | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")

echo "Running commands from: $CMD_FILE"
echo "Log for job ID: $SLURM_ARRAY_TASK_ID"

# Execute each command and capture output
while IFS= read -r cmd; do
    echo "Executing: $cmd"
    eval "$cmd"
done < "$CMD_FILE"


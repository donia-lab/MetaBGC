#!/bin/bash

declare -a FILE_NAME_ARR

IFS=$'\r\n' GLOBIGNORE='*' command eval 'FILE_NAME_ARR=($(cat pgap_all_samples))'

NUMFASTQ=${#FILE_NAME_ARR[@]}
echo "Number of runs: $NUMFASTQ"
ZBNUMFASTQ=$(($NUMFASTQ - 1))

FILE_NAMES=$( IFS=:; printf '%s' "${FILE_NAME_ARR[*]}" )
export ALL_SAMPLE_ARRAY=${FILE_NAMES}
# now submit to SLURM
if [ $ZBNUMFASTQ -ge 0 ]; then
   sbatch --array=0-$ZBNUMFASTQ runPGAP.sh
fi


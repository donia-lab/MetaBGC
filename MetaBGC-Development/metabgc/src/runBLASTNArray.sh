#!/bin/bash

declare -a FILE_NAME_ARR

declare -a BIN_NAME_ARR

IFS=$'\r\n' GLOBIGNORE='*' command eval 'FILE_NAME_PATHS=($(cat max_sample_assemblies.txt))'

IFS=$'\r\n' GLOBIGNORE='*' command eval 'BIN_FILE_PATHS=($(cat bin_paths.txt))'

for FILE in "${FILE_NAME_PATHS[@]}"
do 
   echo -e "$FILE"
   FILE_NAME_ARR+=("$FILE")
done

for BIN in "${BIN_FILE_PATHS[@]}"
do
   echo -e "Bin:$BIN"
   BIN_NAME_ARR+=("$BIN")
done


NUMFASTQ=${#FILE_NAME_ARR[@]}
echo "Number of Samples: $NUMFASTQ"
ZBNUMFASTQ=$(($NUMFASTQ - 1))

FILE_NAMES=$( IFS=:; printf '%s' "${FILE_NAME_ARR[*]}" )
export FILES=${FILE_NAMES}

BIN_NAMES=$( IFS=:; printf '%s' "${BIN_NAME_ARR[*]}" )
export BINS=${BIN_NAMES}

# now submit to SLURM
if [ $ZBNUMFASTQ -ge 0 ]; then
   sbatch --array=0-$ZBNUMFASTQ runBLASTN.sbatch 
fi



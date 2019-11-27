#!/bin/bash

declare -a FILE_NAME_ARR
for FILE in /projects/DONIA/abiswas/MetaBGCRuns/AcbK-homologs/search/reads/*.fasta ;
do 
   FILE_NAME_ARR+=("$FILE")
   echo -e $FILE
done

NUMFASTQ=${#FILE_NAME_ARR[@]}
echo "Number of Bins: $NUMFASTQ"
ZBNUMFASTQ=$(($NUMFASTQ - 1))

FILE_NAMES=$( IFS=:; printf '%s' "${FILE_NAME_ARR[*]}" )
export FILES=${FILE_NAMES}
# now submit to SLURM
if [ $ZBNUMFASTQ -ge 0 ]; then
   sbatch --array=0-$ZBNUMFASTQ runBLASTN.sbatch 
fi



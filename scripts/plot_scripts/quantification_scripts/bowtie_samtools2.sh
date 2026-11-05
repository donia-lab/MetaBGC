#!/bin/bash

if [ $# -ne 4 ]; then
    echo Usage: $0 query_file input_directory output_basename threads
    exit 1
fi

module load samtools/1.21

query_file=$1
input_directory=$2
output_basename=$3
threads=$4

echo args...
echo $query_file
echo $input_directory
echo $output_basename
echo $threads

samtools_tmpdir=/scratch/gpfs/$USER/samtools-sort-temp
mkdir -p $samtools_tmpdir

input_args=$(mgpInputArgs.py -d $input_directory -p bowtie)

bowtie2 \
    --end-to-end \
    --very-sensitive \
    --score-min L,-0.6,-0.3 \
    -p ${threads} \
    -x ${query_file} \
    $input_args \
    2> ${output_basename}.err | \
    samtools view -u -S -@ 3 | \
    samtools sort -l 0 -@ 3 -T $samtools_tmpdir \
	     > ${output_basename}.bam

samtools index -b ${output_basename}.bam

samtools idxstats ${output_basename}.bam > ${output_basename}.tsv

# rm ${output_basename}.bam
# rm ${output_basename}.bam.bai

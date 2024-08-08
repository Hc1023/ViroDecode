#!/bin/bash
# BAMCLIPPER
# Arguments:
#   $1  ${OutDir} 
#   $2  ${bedpe}
#   $3  ${thread}

samtools index $1/sort.bam

bamclipper=/workfolder/bamclipper/bamclipper.sh
bam=$1/sort.bam

${bamclipper} -b ${bam} -p $2 -n $3

# output: ${OutDir}/sort.primerclipped.bam

mv $1/sort.bam $1/sort.unclipped.bam
mv $1/sort.bam.bai $1/sort.unclipped.bam.bai
mv sort.primerclipped.bam $1/sort.bam
mv sort.primerclipped.bam.bai $1/sort.bam.bai
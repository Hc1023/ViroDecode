#!/bin/bash
# READS FILTERING
# Arguments:
#   $1  ${thread} 
#   $2  ${OutDir}

samtools view -@$1 -F12 -q30 -b $2/sort0.bam|\
    samtools sort -@$1 - > $2/sort.bam
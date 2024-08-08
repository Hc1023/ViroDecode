#!/bin/bash
# PILON
# Arguments:
#   $1  ${OutDir} 
#   $2  ${ref}

samtools index $1/sort.bam

pilon \
    --genome $2 --frags $1/sort.bam \
    --output $1/pilon_polished \
    --changes --mindepth 1 --fix all &> $1/pilon.log

#!/bin/bash
# ASSEMBLY
# Arguments:
#   $1  ${thread} 
#   $2  ${OutDir}

spades.py -t $1 -1 $2/r1_sub.fq \
    -2 $2/r2_sub.fq -o $2/spades_output_sub
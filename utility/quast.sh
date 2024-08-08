#!/bin/bash
# QUAST
# Arguments:
#   $1  $gff 
#   $2  $assemfa 
#   $3  ${thread} 
#   $4  ${OutDir} 
#   $5  ${myref}

if [ $1 ]; then
    quast $2 \
    -t $3 \
    -R $4/$5 \
    -G $1 \
    -o $4/resquast
else
    quast $2 \
    -t $3 \
    -R $4/$5 \
    -o $4/resquast
fi
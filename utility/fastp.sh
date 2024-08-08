#!/bin/bash
# FASTP: QUALITY CONTROL
# Arguments:
#   $1  ${InDir}/${fq1}
#   $2  ${InDir}/${fq2}
#   $3  ${thread}
#   $4  ${OutDir}/${fq1}
#   $5  ${OutDir}/${fq2}
#   $6  ${OutDir}

fastp -i $1 -I $2 -w $3 \
    -o $4 -O $5 \
    -h $6/fastp.html -j $6/fastp.json
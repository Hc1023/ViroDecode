#!/bin/bash
# BOWTIE2: INDEXING AND MAPPING
# Arguments:
#   $1  ${OutDir}/${myref}
#   $2  ${thread}
#   $3  $OutDir/$fq1
#   $4  $OutDir/$fq2
#   $5  ${OutDir}

bowtie2-build --threads $2 $1 $5/ref
bowtie2 -p $2 -x $5/ref \
    -1 $3 -2 $4 -S $5/out.sam
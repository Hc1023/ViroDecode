#!/bin/bash
# BWA MEM: INDEXING AND MAPPING
# Arguments:
#   $1  ${OutDir}/${myref}
#   $2  ${thread}
#   $3  $OutDir/$fq1
#   $4  $OutDir/$fq2
#   $5  ${OutDir}

### no thread option supported for bwa index
bwa index $1
### -v 1 output only error
bwa mem -t $2 -v 1 -M $1 $3 $4 > $5/out.sam
#!/bin/bash
# SNPEFF: SNP ANNOTATION
# Arguments:
#   $1  ${OutDir}  
#   $2  ${snpEff}
#   $3  ${genome_version}

mkdir -p $1/snpeff
java -Xmx8g -jar $2 $3 \
    $1/snp.vcf > $1/snpeff/snp.ann.vcf \
    -stats $1/snpeff/snpEff_summary.html
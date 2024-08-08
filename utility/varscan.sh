#!/bin/bash
# VARSCAN: SNP CALLING
# Arguments:
#   $1  ${OutDir}/${myref}  
#   $2  ${OutDir}
#   $3  ${varscan}
#   $4  ${min_reads2}
#   $5  ${min_var_freq}
#   $6  ${strand_filter}

samtools mpileup -f $1 $2/sort.bam |\
	java -jar $3 \
	mpileup2snp > $2/snp.vcf \
	--output-vcf 1 --min-reads2 $4 --min-var-freq $5 --strand-filter $6
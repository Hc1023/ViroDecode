#!/bin/bash
# QUALIMAP: MAPPING QUALITY
# Arguments:
#   $1  ${thread} 
#   $2  ${OutDir}

samtools view -@$1 -bS $2/out.sam |
    samtools sort -@$1 - > $2/sort0.bam

qualimap bamqc -nt $1 -bam $2/sort0.bam -outfile result.pdf

mv $2/sort0_stats $2/qualimap

results=$2/qualimap/genome_results.txt
logfile=$2/process.log
sed -n 20,21p ${results} >> ${logfile}
sed -n 94p ${results}  >> ${logfile}
sed -n 104p ${results} >> ${logfile}

len_ref=`tail -3 ${results}|head -1|cut -f 3`
mean_coverage=`tail -3 ${results}|head -1|cut -f 5`
echo "     Reference length: ${len_ref}" >> ${logfile}
echo "     Mean coverage: ${mean_coverage}" >> ${logfile}
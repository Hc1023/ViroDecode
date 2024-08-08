#!/bin/bash
# SUBSAMPLE
# Arguments:
#   $1  ${OutDir} 
#   $2  ${logfile} 
#   $3  ${coverage} 
#   $4  ${thread}

len_ref=`tail -3 $1/qualimap/genome_results.txt|head -1|cut -f 3`
mean_coverage=`tail -3 $1/qualimap/genome_results.txt|head -1|cut -f 5`
ratio=$(bc <<< "scale=6;$3/${mean_coverage}")
echo "Reference length: ${len_ref}" >> $2
echo "Mean coverage: ${mean_coverage}" >> $2
echo "Coverage goal: $3" >> $2
echo "Ratio: ${ratio}" >> $2

if (( $(echo "${ratio} < 1" |bc -l) ));then
    echo "Subsample: yes" >> $2
    sambamba view -h -t $4 -s $ratio -f bam \
        --subsampling-seed=100 $1/sort.bam -o $1/sub.bam   
else
    echo "Subsample: no" >> $2
    cp $1/sort.bam $1/sub.bam

fi

samtools fastq -1 $1/r1_sub.fq -2 $1/r2_sub.fq \
    -n $1/sub.bam
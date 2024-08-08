#!/bin/bash
# run directly on server

InDir=$1
if [ ! -d "$2/virodecode" ]; then mkdir -p "$2/virodecode"; fi
OutDir=$2/virodecode
fq1=$3
fq2=$4
ref=$5
thread=$6
align=$7
coverage=$8
gff=$9
genome_version=${10}
amplicon=${11}
bedpe=${12}
ass=${13}
min_reads2=${14}
min_var_freq=${15}
strand_filter=${16}
clean=${17}

cd ${OutDir}

# varscan=./VarScan.v2.3.9.jar
varscan=/workfolder/VarScan.v2.3.9.jar
# util=./utility
util=/workfolder/utility
# snpEff=./snpEff/snpEff.jar
cp -r /workfolder/snpEff ${OutDir}/workfolder
snpEff=${OutDir}/workfolder/snpEff.jar

myref=`basename ${ref}`
cp -n ${ref} ${OutDir}

touch ${OutDir}/process.log
logfile=${OutDir}/process.log

echo `date +%Y-%m-%d-%X`'===============BEGIN=================' > ${logfile}

## fastp
echo `date +%Y-%m-%d-%X`'=======FASTP: QUALITY CONTROL========' >> ${logfile}
$util/fastp.sh ${InDir}/${fq1} ${InDir}/${fq2} ${thread} \
    ${OutDir}/${fq1} ${OutDir}/${fq2} ${OutDir}

## alignment
echo `date +%Y-%m-%d-%X`'=====BOWTIE2/BWA-MEM: ALIGNMENT======' >> ${logfile}

if [ "$align" = "1" ]; then
    echo 'Alignment using BOWTIE2' >> ${logfile}
    align=bowtie2.sh
elif [ "$align" = "2" ]; then
    echo 'Alignment using BWA MEM' >> ${logfile}
    align=bwa-mem.sh
fi
$util/$align ${OutDir}/${myref} ${thread} $OutDir/$fq1 $OutDir/$fq2 ${OutDir}

## qualimap
echo `date +%Y-%m-%d-%X`'======QUALIMAP: MAPPING QUALITY======' >> ${logfile}
$util/qualimap.sh ${thread} ${OutDir}

## filtering
echo `date +%Y-%m-%d-%X`'===========READS FILTERING===========' >> ${logfile}
$util/filtering.sh ${thread} ${OutDir}

## bamclipper
if [ "${bedpe}" != "0" ]; then
echo `date +%Y-%m-%d-%X`'===========PRIMER REMOVING===========' >> ${logfile}
$util/bamclip.sh ${OutDir} ${bedpe} ${thread}
fi

## calling variants
echo `date +%Y-%m-%d-%X`'========VARSCAN: SNP CALLING=========' >> ${logfile}
$util/varscan.sh ${OutDir}/${myref} ${OutDir} ${varscan} ${min_reads2} ${min_var_freq} ${strand_filter}

## snpEff annotation
if [ "${genome_version}" != "0" ]; then
echo `date +%Y-%m-%d-%X`'=======SNPEFF: SNP ANNOTATION========' >> ${logfile}
$util/snpeff.sh ${OutDir} ${snpEff} ${genome_version}
fi

if [ "${ass}" = "0" ]; then
if [ "${clean}" = "1" ]; then
rm -rf ${OutDir}/${fq1} ${OutDir}/${fq2}
rm -rf ${OutDir}/*.fna* ${OutDir}/*.bt2
rm -rf ${OutDir}/*.json
rm -rf ${OutDir}/*.bam* ${OutDir}/*.sam
rm -rf ${OutDir}/workfolder
fi
echo `date +%Y-%m-%d-%X`'================END==================' >> ${logfile}
exit
fi

## subsample
echo `date +%Y-%m-%d-%X`'=============SUBSAMPLE===============' >> ${logfile}
$util/subsample.sh ${OutDir} ${logfile} ${coverage} ${thread}

if [ "$amplicon" = "1" ]; then
echo `date +%Y-%m-%d-%X`'==========PIlON ASSEMBLY=============' >> ${logfile}
echo "Pilon polish for amplicon..." >> ${logfile}
$util/pilon.sh ${OutDir} ${OutDir}/${myref}
assemfa=${OutDir}/pilon_polished.fasta
else
echo `date +%Y-%m-%d-%X`'=============ASSEMBLY================' >> ${logfile}
$util/spades.sh ${thread} ${OutDir}
assemfa=${OutDir}/spades_output_sub/scaffolds.fasta
fi

## assembly results quast (visualization)
echo `date +%Y-%m-%d-%X`'===============QUAST=================' >> ${logfile}
$util/quast.sh $gff $assemfa ${thread} ${OutDir} ${myref}

if [ "${clean}" = "1" ]; then
rm -rf ${OutDir}/${fq1} ${OutDir}/${fq2}
rm -rf ${OutDir}/*.fna* ${OutDir}/*.bt2
rm -rf ${OutDir}/*.json
rm -rf ${OutDir}/*.bam* ${OutDir}/*.sam ${OutDir}/*.fq
mkdir -p ${OutDir}/assembly
if [ "$amplicon" = "1" ]; then
mv ${OutDir}/pilon* ${OutDir}/assembly/
else
mv ${OutDir}/spades_output_sub/*.fasta ${OutDir}/assembly/
fi
rm -rf ${OutDir}/spades_output_sub
rm -rf ${OutDir}/workfolder
fi
echo `date +%Y-%m-%d-%X`'================END==================' >> ${logfile}

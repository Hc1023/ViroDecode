#!/bin/bash
#SBATCH --job-name=ass
#SBATCH -o output/%x_%A_%a.out
#SBATCH -e output/%x_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00
#SBATCH --array=1-2

source activate viro

fq=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /share/home/jianglab/huangsisi/usr/Lab_Analysis/Virus_Data/nova576-2021-04-19/sampid.txt)

InDir=/share/home/jianglab/huangsisi/usr/Virus/nova-2021-04-19
fq1=${fq}_R1.fastq.gz
fq2=${fq}_R2.fastq.gz
OutDir=/share/home/jianglab/huangsisi/usr/viro/21-04-19/220915assembly/${fq}
if [ ! -d $OutDir ]; then mkdir -p $OutDir; fi

ref=/share/home/jianglab/huangsisi/usr/viro/ref/NC_045512.2.fna

thread=${SLURM_CPUS_PER_TASK}
cd ${OutDir}

varscan=/share/home/jianglab/huangsisi/bin/VarScan.v2.3.9.jar
snpEff=/share/home/jianglab/huangsisi/bin/snpEff/snpEff.jar


myref=`basename ${ref}`
cp -n ${ref} ${OutDir}

touch ${OutDir}/process.log
logfile=${OutDir}/process.log

echo `date +%Y-%m-%d-%X`'===============BEGIN=================' > ${logfile}
## fastp
echo `date +%Y-%m-%d-%X`'=======FASTP: QUALITY CONTROL========' >> ${logfile}
fastp -i ${InDir}/${fq1} -I ${InDir}/${fq2} -w ${thread} \
    -o ${OutDir}/${fq1} -O ${OutDir}/${fq2} \
    -h ${OutDir}/fastp.html -j ${OutDir}/fastp.json

## alignment
echo `date +%Y-%m-%d-%X`'=====BOWTIE2/BWA-MEM: ALIGNMENT======' >> ${logfile}

echo 'Alignment using BOWTIE2' >> ${logfile}
bowtie2-build --threads ${thread} ${OutDir}/${myref} ${OutDir}/ref
bowtie2 -p ${thread} -x ${OutDir}/ref \
    -1 $OutDir/$fq1 -2 $OutDir/$fq2 -S ${OutDir}/out.sam

## filtering
echo `date +%Y-%m-%d-%X`'===========READS FILTERING===========' >> ${logfile}
samtools view -@${thread} -F12 -q30 -b ${OutDir}/sort0.bam|\
    samtools sort -@${thread} - > ${OutDir}/sort.bam

echo `date +%Y-%m-%d-%X`'===========PRIMER REMOVING===========' >> ${logfile}
## bamclipper
samtools index ${OutDir}/sort.bam

### Dir of bamclipper
bamclipper=/share/home/jianglab/huangsisi/bin/bamclipper/bamclipper.sh
bam=${OutDir}/sort.bam

${bamclipper} -b ${bam} -p ${bedpe} -n ${thread}

samtools fastq -1 ${OutDir}/r1.fq -2 ${OutDir}/r2.fq \
    -n ${OutDir}/sort.primerclipped.bam

echo `date +%Y-%m-%d-%X`'=============ASSEMBLY================' >> ${logfile}
spades.py -t ${thread} -1 ${OutDir}/r1.fq \
    -2 ${OutDir}/r2.fq -o ${OutDir}/spades_output

echo `date +%Y-%m-%d-%X`'================END==================' >> ${logfile}


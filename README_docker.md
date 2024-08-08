# ViroDecode

There is a plethora of research on genetic mutations detection and whole-genome assembly during the dynamic spread of the SARS-CoV-2 pandemic. ViroDecode is a containerized integrated bioinformatic tool for DECODING virome, providing relevant analysis reports. 

- Tool basic functions: 
    - Variant calling and mutation frequency
    - Genome assembly
- Additions:
    - Reads quality
    - Alignment evaluation
    - Assembly quality evaluation and visualization
    - SNP annotation

For genome assembly, difficulties in dealing with ultra-deep sequencing data (with some exceeding 2 million) call for a relevant analysis pipeline. High coverage sequencing data for small genomes, such as virus, bacteria, or bacterial artificial nucleic acid clones, is very common. We adopt the idea of taking subsets of the excessive sequencing data and combine some polishing steps to improve workflows. 

Although ViroDecode is particularly designed for analysing coronavirus, it can similarly be used in other genome sequencing studies, especially for those ultra-deep sequencing data.

In general, we integrate the functions of the following software, and design a robust analysing pipeline.
- fastp 0.20.1
- bowtie2 2.3.5.1
- bwa 0.7.17-r1188
- samtools 1.9
- qualimap 2.2.2-dev
- VarScan 2.3.9
- snpEff 5.0e
- SPAdes 3.13.0
- Pilon 1.24
- QUAST 5.0.2

However, as a Docker image, ViroDecode is user-friendly: no environment building, no dependencies installing, ‘no’ coding…

## Input and Output

With the help of [run_virodecode_docker.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_docker.py) or [run_virodecode_singularity.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_singularity.py), ViroDecode takes as input paired-end reads in FASTQ format. Other options listed as follows.

The complete output include the following files

- `snp.vcf`: variants and the frequency of the mutation.
- `fastp.html`: reads quality report.
- `qualimap/`: results of `qualimap`, reads mapping evaluation.
- `resquast/`: results of `quast`, with visualization of the alignment by `icarus`.
- `snpeff/`: results of `snpEff`, genetic variant annotation and functional effect.
- `assembly/`: the genome assembly fasta file.
- `process.log`: log file.

## Usage

We have created two helper scripts for Docker and Singularity, respectively. To download, simply clone the ViroDecode repository from GitHub. From the commandline type:

```sh
git clone https://github.com/Hc1023/ViroDecode.git
```

### Using Docker (Require root privilege)

Run ViroDecode via the Docker image. For your convenience, we have created a script, [run_virodecode_docker.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_docker.py), which handles Docker-related issues such as mounting files for you and is cross-platform compatible. 

1. Install Docker and Python on your machine.
2. Pull the Docker image.
```sh
sudo docker pull sisih/virodecode:1
```
3. Run the python script.
```sh
python run_virodecode_docker.py -i sisih/virodecode:1 -1 fq1 -2 fq2
```

### Using Singularity (Recommended for HPC)

Singularity can also be used with Docker images. For security reasons, many high performance computing resources disallow Docker usage. Instead, many of these platforms support Singularity, which is very similar to Docker. Thus, if you are using a compute resource that supports Singularity, you can use the [run_virodecode_singularity.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_singularity.py) script.

1. Load Singularity and Python on your computer cluster. This must be done each time you log in (or you can add it to your `.bashrc` for convenience). 
2. Import the Docker image into a Singularity Image.
```sh
singularity pull docker://sisih/virodecode:1
```
3. Run the python script.
```sh
python run_virodecode_singularity.py -i virodecode_1.sif -1 fq1 -2 fq2
```

### Command line options

```{sh}
usage: run_virodecode_singularity.py [-h] -i IMAGE -1 READ1 -2 READ2
                                     [-o OUTDIR] [-r REF] [-t THREAD]
                                     [-a ALIGN] [-c COVERAGE] [-g GFF]
                                     [-G GENOME_VERSION] [-p AMPLICON]
                                     [-A ASSEMBLY]

optional arguments:
  -h, --help            show this help message and exit
  -i IMAGE, --image IMAGE
                        Image path
  -1 READ1, --read1 READ1
                        Provide read1 file path
  -2 READ2, --read2 READ2
                        Provide read2 file path
  -o OUTDIR, --outDir OUTDIR
                        Provide path to a directory to place the output
  -r REF, --ref REF     Reference genome
  -t THREAD, --thread THREAD
                        Thread
  -a ALIGN, --align ALIGN
                        Alignment tool 1: BOWTIE2 2: BWA MEM
  -c COVERAGE, --coverage COVERAGE
                        Subsample reads if the mean coverage is higher than
                        xxx
  -g GFF, --gff GFF     Gff file, for Icarus visualization
  -G GENOME_VERSION, --genome_version GENOME_VERSION
                        Genome version (e.g. NC_045512.2) for SNP annotation
                        by snpEff
  -p AMPLICON, --amplicon AMPLICON
                        Set 1 if amplicon sequencing
  -A ASSEMBLY, --assembly ASSEMBLY
                        Set 1 if needing assembly results (fasta)
```

If you are analysing coronavirus, the required arguments are `-i` (image path), `-1` and `-2` (paired-end reads). Otherwise, you should provide the reference genome for `-r`. However, `gff` file and `genome_version` options are optional.

### Example

```sh
python ${DIR_OF_RUN_VIRODECODE_SINGULARITY.PY} \
    -i ${DIR_OF_VIRODECODE_1.SIF} \
    -1 ${DIR_OF_READ1} -2 ${DIR_OF_READ2}
```

- `process.log`

```sh
2022-05-20-18:49:06===============BEGIN=================
2022-05-20-18:49:06=======FASTP: QUALITY CONTROL========
2022-05-20-18:49:08=====BOWTIE2/BWA-MEM: ALIGNMENT======
Alignment using BOWTIE2
2022-05-20-18:49:23======QUALIMAP: MAPPING QUALITY======
     number of reads = 174,020
     number of mapped reads = 142,794 (82.06%)
     There is a 99.22% of reference with a coverageData >= 20X
     There is a 99.22% of reference with a coverageData >= 30X
     Reference length: 29903
     Mean coverage: 700.6409724776778
2022-05-20-18:49:32===========READS FILTERING===========
2022-05-20-18:49:33========VARSCAN: SNP CALLING=========
2022-05-20-18:49:55=======SNPEFF: SNP ANNOTATION========
2022-05-20-18:50:20=============SUBSAMPLE===============
Reference length: 29903
Mean coverage: 700.6409724776778
Coverage goal: 5000
Ratio: 7.136322
Subsample: no
2022-05-20-18:50:21==========PIlON ASSEMBLY=============
Pilon polish for amplicon...
2022-05-20-18:50:22===============QUAST=================
2022-05-20-18:50:24================END==================
```

- Output

```sh
$ tree virodecode/
├── assembly
│   ├── before_rr.fasta
│   ├── contigs.fasta
│   ├── pilon.log
│   └── scaffolds.fasta
├── fastp.html
├── process.log
├── qualimap
│   ├── genome_results.txt
│   └── result.pdf
├── resquast
│   ├── aligned_stats
│   │   ├── cumulative_plot.pdf
│   │   ├── NAx_plot.pdf
│   │   └── NGAx_plot.pdf
│   ├── basic_stats
│   │   ├── coverage_histogram.pdf
│   │   ├── cumulative_plot.pdf
│   │   ├── GC_content_plot.pdf
│   │   ├── gc.icarus.txt
│   │   ├── NGx_plot.pdf
│   │   ├── Nx_plot.pdf
│   │   ├── scaffolds_coverage_histogram.pdf
│   │   └── scaffolds_GC_content_plot.pdf
│   ├── contigs_reports
│   │   ├── all_alignments_scaffolds.tsv
│   │   ├── contigs_report_scaffolds.mis_contigs.info
│   │   ├── contigs_report_scaffolds.stderr
│   │   ├── contigs_report_scaffolds.stdout
│   │   ├── contigs_report_scaffolds.unaligned.info
│   │   ├── minimap_output
│   │   │   ├── scaffolds.coords
│   │   │   ├── scaffolds.coords.filtered
│   │   │   ├── scaffolds.coords_tmp
│   │   │   ├── scaffolds.sf
│   │   │   ├── scaffolds.unaligned
│   │   │   └── scaffolds.used_snps.gz
│   │   ├── misassemblies_frcurve_plot.pdf
│   │   ├── misassemblies_plot.pdf
│   │   ├── misassemblies_report.tex
│   │   ├── misassemblies_report.tsv
│   │   ├── misassemblies_report.txt
│   │   ├── scaffolds.mis_contigs.fa
│   │   ├── transposed_report_misassemblies.tex
│   │   ├── transposed_report_misassemblies.tsv
│   │   ├── transposed_report_misassemblies.txt
│   │   ├── unaligned_report.tex
│   │   ├── unaligned_report.tsv
│   │   └── unaligned_report.txt
│   ├── genome_stats
│   │   ├── features_cumulative_plot.pdf
│   │   ├── features_frcurve_plot.pdf
│   │   ├── genome_info.txt
│   │   ├── scaffolds_gaps.txt
│   │   └── scaffolds_genomic_features_gene.txt
│   ├── icarus.html
│   ├── icarus_viewers
│   │   ├── alignment_viewer.html
│   │   └── contig_size_viewer.html
│   ├── quast.log
│   ├── report.html
│   ├── report.pdf
│   ├── report.tex
│   ├── report.tsv
│   ├── report.txt
│   ├── transposed_report.tex
│   ├── transposed_report.tsv
│   └── transposed_report.txt
├── snpeff
│   ├── snp.ann.vcf
│   ├── snpEff_summary.genes.txt
│   └── snpEff_summary.html
└── snp.vcf
```

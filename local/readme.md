# ViroDecode

The workflow was developed within the miniconda3:4.8.2 Docker environment (https://hub.docker.com/r/continuumio/miniconda3, accessed on January 10, 2023). 

ViroDecode was composed of two main functions: 
- SNP calling
- Genome assembly

The workflow was optimized and standardized for ultra-deep sequencing and amplicon sequencing data, with a focus on research on the SARS-CoV-2 virus. The filtered bam files from the SNP calling workflow were used for genome assembly. The alignment tools bwa-mem (0.7.17-r1188) (Li, 2013) and bowtie2 (2.3.5.1) (Langmead and Salzberg, 2012) were both available for use. For efficient genome assembly by SPAdes (3.13.0) (Bankevich et al., 2012), a downsampling step was added for ultra-high coverage sequencing data using sambamba (0.7.1) (Tarasov et al., 2015). For amplicon sequencing, assembly was conducted using pilon (1.24) (Walker et al., 2014). Relevant analysis reports and visualization modules were also provided, including qualimap (2.2.2-dev) (Okonechnikov et al., 2016) to analyze alignment metrics, snpEff (5.0e) (Cingolani et al., 2012) to annotate SNPs, and QUAST (5.0.2) (Gurevich et al., 2013) to inspect and visualize assembly results.

## Input and Output

With the help of [run_virodecode_docker.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_docker.py) or [run_virodecode_singularity.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_singularity.py), ViroDecode takes as input paired-end reads in FASTQ format. Other options are listed as follows.

The complete output include the following files

- `snp.vcf`: variants and the frequency of mutations
- `assembly/`: genome assembly results
- `fastp.html`: reads quality report
- `qualimap/`: results of `qualimap`, reads mapping evaluation
- `snpeff/`: results of `snpEff`, genetic variant annotation and functional effect
- `resquast/`: results of `quast`, with visualization of the alignment by `icarus`
- `process.log`: log file

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
python3 run_virodecode_docker.py -i sisih/virodecode:1 -1 fq1 -2 fq2
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
python3 run_virodecode_singularity.py -i virodecode_1.sif -1 fq1 -2 fq2
```

### Command line options

```{sh}
usage: run_virodecode_singularity.py [-h] -i IMAGE -1 READ1 -2 READ2
                                     [-o OUTDIR] [-r REF] [-t THREAD]
                                     [-a ALIGN] [-c COVERAGE] [-g GFF]
                                     [-G GENOME_VERSION] [-p AMPLICON]
                                     [-b BEDPE] [-A ASSEMBLY]
                                     [--min_reads2 MIN_READS2]
                                     [--min_var_freq MIN_VAR_FREQ]
                                     [--strand_filter STRAND_FILTER]
                                     [--clean CLEAN]

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
                        Alignment tool 1: BOWTIE2 (default) 2: BWA MEM
  -c COVERAGE, --coverage COVERAGE
                        Subsample reads if the mean coverage is higher than
                        xxx
  -g GFF, --gff GFF     Gff file, for Icarus visualization
  -G GENOME_VERSION, --genome_version GENOME_VERSION
                        Genome version (e.g. NC_045512.2) for SNP annotation
                        by snpEff
  -p AMPLICON, --amplicon AMPLICON
                        Set to 1 if amplicon sequencing
  -b BEDPE, --bedpe BEDPE
                        BEDPE file of primer pair locations
  -A ASSEMBLY, --assembly ASSEMBLY
                        Set to 1 if needing assembly results (fasta)
  --min_reads2 MIN_READS2
                        VarScan parameter: Was supported by at least the
                        minimum number of supporting reads
  --min_var_freq MIN_VAR_FREQ
                        VarScan parameter: Meets the minimum allele frequency
                        threshold
  --strand_filter STRAND_FILTER
                        VarScan parameter: Passes a basic strand-bias filter
                        if set to 1
  --clean CLEAN         Set to 0 if needing complete processing files
```

If you are analysing coronavirus, the required arguments are `-i` (image path), `-1` and `-2` (paired-end reads). Otherwise, you should provide the reference genome for `-r`. However, `gff` file and `genome_version` options are optional.

### Example

- Call variants

```sh
python3 ${DIR_OF_RUN_VIRODECODE_SINGULARITY.PY} \
    -i ${DIR_OF_VIRODECODE_1.SIF} \
    -1 ${DIR_OF_READ1} -2 ${DIR_OF_READ2}
```

- Call variants from amplicon sequencing data

```sh
python3 ${DIR_OF_RUN_VIRODECODE_SINGULARITY.PY} \
    -p 1 -b ${DIR_OF_BEDPE} \
    -i ${DIR_OF_VIRODECODE_1.SIF} \
    -1 ${DIR_OF_READ1} -2 ${DIR_OF_READ2}
```

- Call variants from amplicon sequencing data, assemble, and reserve all processing files

```sh
python3 ${DIR_OF_RUN_VIRODECODE_SINGULARITY.PY} \
    -A 1 \
    -p 1 -b ${DIR_OF_BEDPE} \
    -i ${DIR_OF_VIRODECODE_1.SIF} \
    -1 ${DIR_OF_READ1} -2 ${DIR_OF_READ2} \
    -- clean 0
```

- Output

```sh
$ tree -L 3
.
└── virodecode
    ├── assembly
    │   ├── pilon.log
    │   ├── pilon_polished.changes
    │   └── pilon_polished.fasta
    ├── fastp.html
    ├── process.log
    ├── qualimap
    │   ├── genome_results.txt
    │   ├── result.pdf
    │   └── sort0_stats
    ├── resquast
    │   ├── aligned_stats
    │   ├── basic_stats
    │   ├── contigs_reports
    │   ├── genome_stats
    │   ├── icarus.html
    │   ├── icarus_viewers
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

- `process.log`

```sh
2022-06-07-17:38:54===============BEGIN=================
2022-06-07-17:38:54=======FASTP: QUALITY CONTROL========
2022-06-07-17:40:16=====BOWTIE2/BWA-MEM: ALIGNMENT======
Alignment using BOWTIE2
2022-06-07-17:46:22======QUALIMAP: MAPPING QUALITY======
     number of reads = 21,105,690
     number of mapped reads = 19,724,086 (93.45%)
     There is a 97.36% of reference with a coverageData >= 20X
     There is a 97.34% of reference with a coverageData >= 30X
     Reference length: 29903
     Mean coverage: 80190.87329030532
2022-06-07-18:04:12===========READS FILTERING===========
2022-06-07-18:04:21========VARSCAN: SNP CALLING=========
2022-06-07-18:19:57=======SNPEFF: SNP ANNOTATION========
2022-06-07-18:20:29=============SUBSAMPLE===============
Reference length: 29903
Mean coverage: 80190.87329030532
Coverage goal: 5000
Ratio: .062351
Subsample: yes
2022-06-07-18:20:28=============ASSEMBLY================
2022-06-07-18:21:28===============QUAST=================
2022-06-07-18:22:25================END==================
```

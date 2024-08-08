FROM continuumio/miniconda3:4.8.2    

RUN conda install -y python=2
RUN conda install -y -c bioconda samtools=1.9
RUN conda install -y -c bioconda fastp=0.20.1
RUN conda install -y -c bioconda bowtie2=2.3.5.1
RUN conda install -y -c bioconda bwa=0.7.17
RUN conda install -y -c conda-forge tbb=2020.2
RUN conda install -y -c bioconda spades=3.13.0

RUN conda install -y -c bioconda qualimap=2.2.2
RUN conda install -y -c conda-forge bc
RUN conda install -y -c bioconda sambamba
RUN conda install -y -c bioconda pilon=1.24
RUN conda install -y -c bioconda/label/cf201901 quast
RUN conda install -y -c conda-forge parallel

RUN mkdir -p /workfolder
RUN mkdir -p /ref

COPY ./main.sh /workfolder/

COPY ./VarScan.v2.3.9.jar /workfolder/
COPY ./snpEff /workfolder/snpEff
COPY ./bamclipper /workfolder/bamclipper
COPY ./utility /workfolder/utility

COPY ./NC_045512.2.fna /ref/
COPY ./sequence.gff3 /ref/

CMD [ "bash", "/workfolder/main.sh" ]
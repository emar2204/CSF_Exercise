FROM mambaorg/micromamba:1.5.8

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.12 \
    fastqc=0.12.1 \
    fastp=1.1.0 \
    bowtie2 \
    samtools \
    metaphlan \
    humann \
    diamond \
 && micromamba clean -a -y

WORKDIR /metagenomic_exercise
COPY pipeline.sh ./pipeline.sh

CMD ["bash", "./pipeline.sh"]

#!/usr/bin/env bash
set -ue
if [ ! -f chrM.fa ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrM.fa.gz
    gunzip chrM.fa.gz
fi
if [ ! -f chrM.ncbiRefSeq.gtf ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
    gzip -cfd ce11.ncbiRefSeq.gtf.gz | grep '^chrM\s' >chrM.ncbiRefSeq.gtf
fi

if [ ! -f isoform_depth.tsv ]; then
    python -m yasim generate_gene_depth \
        -g chrM.ncbiRefSeq.gtf \
        -o gene_depth.tsv \
        -d 20
    python -m yasim generate_isoform_depth \
        -g chrM.ncbiRefSeq.gtf \
        -d gene_depth.tsv \
        -o isoform_depth.tsv \
        --alpha 4
fi
if [ ! -f isoform_low_depth.tsv ]; then
    python -m yasim generate_gene_depth \
        -g chrM.ncbiRefSeq.gtf \
        -o gene_low_depth.tsv \
        -d 5
    python -m yasim generate_isoform_depth \
        -g chrM.ncbiRefSeq.gtf \
        -d gene_low_depth.tsv \
        -o isoform_low_depth.tsv \
        --alpha 4
fi
if [ ! -f chrm_trans.fa ]; then
    python -m labw_utils.bioutils transcribe \
        -g chrM.ncbiRefSeq.gtf \
        -f chrM.fa \
        -o chrm_trans.fa
fi
python -m labw_utils.bioutils describe_gtf chrM.ncbiRefSeq.gtf
if [ ! -f chrm_pbsim3_errhmm.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M errhmm \
        -e pbsim \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_errhmm \
        -j 40
    python -m labw_utils.bioutils describe_fastq --input chrm_pbsim3_errhmm.fq --input_fmt fastq
fi
if [ ! -f chrm_pbsim3_qshmm.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e pbsim \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_qshmm \
        -j 40
    python -m labw_utils.bioutils describe_fastq --input chrm_pbsim3_qshmm.fq --input_fmt fastq
fi
if [ ! -f chrm_pbsim3_clr.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e pbsim \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_clr \
        --ccs_pass 1 \
        -j 40
    python -m labw_utils.bioutils describe_fastq --input chrm_pbsim3_clr.fq --input_fmt fastq
fi
if [ ! -f chrm_pbsim3_ccs.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e pbsim \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_ccs \
        --ccs_pass 10 \
        -j 40
    python -m labw_utils.bioutils describe_fastq --input chrm_pbsim3_ccs.fq --input_fmt fastq
fi

#!/usr/bin/env bash

if [ ! -f ce11.ncbiRefSeq.chr1.gtf ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
    gunzip ce11.ncbiRefSeq.gtf.gz
    grep -i '^chrI\s' <ce11.ncbiRefSeq.gtf >ce11.ncbiRefSeq.chr1.gtf
fi
if [ ! -f ce11.chr1.fa ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz
    zcat chrI.fa.gz >ce11.chr1.fa
fi

python -m yasim generate_as_events \
    -f ce11.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq.chr1.as.gtf \
    -c 5
python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_gene_depth.tsv \
    -d 5
python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -d ce11_gene_depth.tsv \
    -o ce11_isoform_depth.tsv
python -m labw_utils.bioutils transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_trans_as.fa
python -m yasim pbsim3 \
    -F ce11_trans_as.fa.d \
    --hmm_method errhmm \
    --hmm_model RSII \
    -j 40 \
    -d ce11_isoform_depth.tsv \
    -o pbsim3_mode

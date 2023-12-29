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
if [ ! -f chrm_trans.fa ]; then
    python -m labw_utils.bioutils transcribe \
        -g chrM.ncbiRefSeq.gtf \
        -f chrM.fa \
        -o chrm_trans.fa
fi
if [ ! -f isoform_low_depth.tsv ]; then
    python -m yasim generate_gene_depth \
        -g chrM.ncbiRefSeq.gtf \
        -o gene_low_depth.tsv \
        -d 5 \
        --low_cutoff 1
    python -m yasim generate_isoform_depth \
        -g chrM.ncbiRefSeq.gtf \
        -d gene_low_depth.tsv \
        -o isoform_low_depth.tsv \
        --alpha 4 \
        --low_cutoff 1
fi

if [ ! -f chrm_ccs_isoseq.fq ]; then
    python -m yasim pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_ccs_isoseq \
        -j 40 \
        --ccs_pass 10 \
        --preserve_intermediate_files
fi
if [ ! -f chrm_ccs_isoseq.ccs.bam ]; then
    python -m yasim_scripts merge_pbccs \
        --out chrm_ccs_isoseq.ccs.bam \
        --input_bam_glob 'chrm_ccs_isoseq.d/*/tmp*.ccs.bam'
    pbindex chrm_ccs_isoseq.ccs.bam
    samtools index chrm_ccs_isoseq.ccs.bam
fi
if [ ! -f chrm_ccs_isoseq.transcripts.xml.hq.bam ]; then
    isoseq3 cluster \
        chrm_ccs_isoseq.ccs.bam \
        chrm_ccs_isoseq.transcripts.xml \
        --log-level INFO \
        --num-threads 40
fi
if [ ! -f chrm_ccs_isoseq.aln.bam ]; then
    pbmm2 align \
        --preset ISOSEQ \
        --sort \
        --log-level INFO \
        chrm_ccs_isoseq.transcripts.xml.hq.bam \
        chrM.fa \
        chrm_ccs_isoseq.aln.bam
fi
if [ ! -f chrm_ccs_isoseq.collapse.gff ]; then
    isoseq3 collapse \
        --do-not-collapse-extra-5exons \
        --log-level INFO \
        chrm_ccs_isoseq.aln.bam \
        chrm_ccs_isoseq.ccs.bam \
        chrm_ccs_isoseq.collapse.gff
fi

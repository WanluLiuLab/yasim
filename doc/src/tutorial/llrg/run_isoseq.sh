#!/usr/bin/env bash
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

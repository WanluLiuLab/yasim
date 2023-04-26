#!/usr/bin/env bash
set -ue
if [ ! -f chrM.fa ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz &>/dev/null
    tar xzvf chromFa.tar.gz
    rm -f chrI.fa chrII.fa chrIII.fa chrIV.fa
fi
if [ ! -f chrM.ncbiRefSeq.gtf ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz &>/dev/null
    gzip -cfd ce11.ncbiRefSeq.gtf.gz | grep '^chrM\s' >chrM.ncbiRefSeq.gtf
fi
function grep_pbar() {
    { grep -v '1%' || true; } |
        { grep -v '26%' || true; } |
        { grep -v '51%' || true; } |
        { grep -v '76%' || true; }
}
if [ ! -f isoform_depth.tsv ]; then
    python -m yasim generate_gene_depth \
        -g chrM.ncbiRefSeq.gtf \
        -o gene_depth.tsv \
        -d 60
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
if [ ! -f chrm_pbsim3_wgs.fq ]; then
    python -m yasim pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy wgs \
        -F chrm_trans.fa.d \
        -d isoform_depth.tsv \
        -o chrm_pbsim3_wgs \
        -j 40
fi
if [ ! -f chrm_pbsim3_trans.fq ]; then
    python -m yasim pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_depth.tsv \
        -o chrm_pbsim3_trans \
        -j 40
fi
python -m labw_utils.bioutils describe_fastq chrm_pbsim3_trans.fq chrm_pbsim3_wgs.fq
python -m labw_utils.bioutils describe_gtf chrM.ncbiRefSeq.gtf
if [ ! -f chrm_pbsim3_errhmm.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M errhmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_errhmm \
        -j 40
    python -m labw_utils.bioutils describe_fastq chrm_pbsim3_errhmm.fq
fi
if [ ! -f chrm_pbsim3_qshmm.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_qshmm \
        -j 40
    python -m labw_utils.bioutils describe_fastq chrm_pbsim3_qshmm.fq
fi
if [ ! -f chrm_pbsim3_clr.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_clr \
        --ccs_pass 1 \
        -j 40
    python -m labw_utils.bioutils describe_fastq chrm_pbsim3_clr.fq
fi
if [ ! -f chrm_pbsim3_ccs.fq ]; then
    python -m yasim pbsim3 \
        -m RSII \
        -M qshmm \
        -e /home/yuzj/bin/pbsim3 \
        --strategy trans \
        -F chrm_trans.fa.d \
        -d isoform_low_depth.tsv \
        -o chrm_pbsim3_ccs \
        --ccs_pass 10 \
        -j 40
    python -m labw_utils.bioutils describe_fastq chrm_pbsim3_ccs.fq
fi

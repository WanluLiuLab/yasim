#!/usr/bin/env bash
set -ueo pipefail
cd "$(readlink -f "$(dirname "${0}")")" || exit 1
. config.sh

[ ! -f "${FASTQ_NAME}".GENE.bam ] && minimap2 -t "${THREAD_NUM}" -a -x splice "${GENE_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_NAME}".GENE.bam

[ ! -f "${FASTQ_NAME}".TRANS.bam ] && minimap2 -t "${THREAD_NUM}" -a "${TRANSCRIPT_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_NAME}".TRANS.bam

[ ! -f "${FASTQ_NAME}".GENE.bam.bai ] && samtools index "${FASTQ_NAME}".GENE.bam
[ ! -f "${FASTQ_NAME}".TRANS.bam.bai ] && samtools index "${FASTQ_NAME}".TRANS.bam
[ ! -f "${FASTQ_NAME}".TRANS.bam.depth.tsv ] && samtools depth "${FASTQ_NAME}".TRANS.bam > "${FASTQ_NAME}".TRANS.bam.depth.tsv

[ ! -f "${FASTQ_NAME}".transcript.depth.tsv ] &&  \
    Rscript R/transform_depth_results.R \
    --libfile R/lib.R \
    --input "${FASTQ_NAME}".TRANS.bam.depth.tsv \
    --fa_stats "${TRANSCRIPT_FA_STATS}" \
    --output "${FASTQ_NAME}".transcript.depth.tsv
#!/usr/bin/env bash
set -ue

GENE_REFERENCE="${1:-/dev/null}"
TRANSCRIPT_REFERENCE="${2:-/dev/null}"
FASTQ_NAME="${3:-/dev/null}"
THREAD_NUM="${4:-8}"

minimap2 -t "${THREAD_NUM}" -a -x splice "${GENE_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_NAME}".GENE.bam

minimap2 -t "${THREAD_NUM}" -a "${TRANSCRIPT_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_NAME}".TRANS.bam

samtools index"${FASTQ_NAME}".GENE.bam
samtools index"${FASTQ_NAME}".TRANS.bam
samtools depth "${FASTQ_NAME}".TRANS.bam | \
xz -9vcf -T"${THREAD_NUM}" > "${FASTQ_NAME}".TRANS.bam.depth.tsv

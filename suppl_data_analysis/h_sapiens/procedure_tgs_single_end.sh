#!/usr/bin/env bash
set -ueo pipefail

cd "$(readlink -f "$(dirname "${0}")")" || exit 1
. shlib/libstr.sh
. shlib/libisopt.sh
. shlib/libdo.sh
. shlib/libprivate.sh


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
    --fa_stats "${TRANSCRIPT_REFERENCE}.stats" \
    --output "${FASTQ_NAME}".transcript.depth.tsv
[ ! -f "${FASTQ_NAME}".stringtie_unguided.gtf ] && \
stringtie "${FASTQ_NAME}".GENE.bam -o "${FASTQ_NAME}".stringtie_unguided.gtf

[ ! -f "${FASTQ_NAME}".stringtie_guided.gtf ] && \
stringtie "${FASTQ_NAME}".GENE.bam \
-G "${GENE_GTF}" \
-o "${FASTQ_NAME}".stringtie_guided.gtf

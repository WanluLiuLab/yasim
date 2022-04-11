#!/usr/bin/env bash
set -ueo pipefail

cd "$(readlink -f "$(dirname "${0}")")" || exit 1
. shlib/libstr.sh
. shlib/libisopt.sh
. shlib/libdo.sh
. shlib/libprivate.sh

if [ -f "${FASTQ_BASE_NAME}".fq ]; then
    FASTQ_NAME="${FASTQ_BASE_NAME}".fq
elif [ -f "${FASTQ_BASE_NAME}".fq.gz ]; then
    FASTQ_NAME="${FASTQ_BASE_NAME}".fq.gz
elif [ -f "${FASTQ_BASE_NAME}".fastq ]; then
    FASTQ_NAME="${FASTQ_BASE_NAME}".fastq
elif [ -f "${FASTQ_BASE_NAME}".fastq.gz ]; then
    FASTQ_NAME="${FASTQ_BASE_NAME}".fastq.gz
else
    errh "No FASTQ found!"
fi

[ ! -f "${FASTQ_BASE_NAME}".GENE.bam ] && minimap2 -t "${THREAD_NUM}" -a -x splice "${GENE_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".GENE.bam

[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam ] && minimap2 -t "${THREAD_NUM}" -a "${TRANSCRIPT_REFERENCE}" "${FASTQ_NAME}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".TRANS.bam

[ ! -f "${FASTQ_BASE_NAME}".GENE.bam.bai ] && samtools index "${FASTQ_BASE_NAME}".GENE.bam
[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam.bai ] && samtools index "${FASTQ_BASE_NAME}".TRANS.bam
[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv ] && samtools depth "${FASTQ_BASE_NAME}".TRANS.bam > "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv

[ ! -f "${FASTQ_BASE_NAME}".transcript.depth.tsv ] &&  \
    Rscript R/transform_depth_results.R \
    --libfile R/lib.R \
    --input "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv \
    --fa_stats "${TRANSCRIPT_REFERENCE}.stats" \
    --output "${FASTQ_BASE_NAME}".transcript.depth.tsv
[ ! -f "${FASTQ_BASE_NAME}".stringtie_unguided.gtf ] && \
stringtie "${FASTQ_BASE_NAME}".GENE.bam -o "${FASTQ_BASE_NAME}".stringtie_unguided.gtf

[ ! -f "${FASTQ_BASE_NAME}".stringtie_guided.gtf ] && \
stringtie "${FASTQ_BASE_NAME}".GENE.bam \
-G "${GENE_GTF}" \
-o "${FASTQ_BASE_NAME}".stringtie_guided.gtf

#!/usr/bin/env bash
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate yasim_c_elegans_as_depth_analysis
set -ueo pipefail

SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}"/shlib/libstr.sh
. "${SHDIR}"/shlib/libisopt.sh
. "${SHDIR}"/shlib/libdo.sh
. "${SHDIR}"/shlib/libprivate.sh

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

. "${SHDIR}"/shlib/libpost_alignment_analysis.sh

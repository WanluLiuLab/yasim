#!/usr/bin/env bash
set -ueo pipefail

cd "$(readlink -f "$(dirname "${0}")")" || exit 1
. shlib/libstr.sh
. shlib/libisopt.sh
. shlib/libdo.sh
. shlib/libprivate.sh

if [ -f "${FASTQ_BASE_NAME}.fq" ] ; then
    FASTQ_NAME="${FASTQ_BASE_NAME}.fq"
    IS_GZ=0
elif [ -f "${FASTQ_BASE_NAME}.fq.gz" ] ; then
    FASTQ_NAME="${FASTQ_BASE_NAME}.fq.gz"
    IS_GZ=1
elif [ -f "${FASTQ_BASE_NAME}.fastq" ] ; then
    FASTQ_NAME="${FASTQ_BASE_NAME}.fastq"
    IS_GZ=0
elif [ -f "${FASTQ_BASE_NAME}.fastq.gz" ] ; then
    FASTQ_NAME="${FASTQ_BASE_NAME}.fastq.gz"
    IS_GZ=1
else
    errh "No FASTQ found!"
fi

if [ -d "${GENE_REFERENCE}_STAR_idx" ]; then
    GENE_REFERENCE="${GENE_REFERENCE}_STAR_idx"
else
    errh "STAR index not found at ${GENE_REFERENCE}_STAR_idx"
fi

if [ ${IS_GZ} -eq 0 ]; then
    [ ! -f "${FASTQ_BASE_NAME}".GENE.bam ] && STAR --runThreadN "${THREAD_NUM}" --genomeDir "${GENE_REFERENCE}" \
    --readFilesIn "${FASTQ_NAME}" --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outFileNamePrefix "${FASTQ_BASE_NAME}" | \
    samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".GENE.bam
else
    [ ! -f "${FASTQ_BASE_NAME}".GENE.bam ] && STAR --runThreadN "${THREAD_NUM}" --genomeDir "${GENE_REFERENCE}" --readFilesCommand zcat \
    --readFilesIn "${FASTQ_NAME}" --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outFileNamePrefix "${FASTQ_BASE_NAME}" | \
    samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".GENE.bam
fi

[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam ] && bwa mem -t 12 "${GENE_REFERENCE}" "${FASTQ_NAME}" | \
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

#!/usr/bin/env bash
set -ueo pipefail

cd "$(readlink -f "$(dirname "${0}")")" || exit 1
. shlib/libstr.sh
. shlib/libisopt.sh
. shlib/libdo.sh
. shlib/libprivate.sh

if [ -f "${FASTQ_BASE_NAME}_1.fq" ] ; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fq"
    IS_GZ = 0
elif [ -f "${FASTQ_BASE_NAME}_1.fq.gz" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fq.gz"
    IS_GZ = 1
elif [ -f "${FASTQ_BASE_NAME}_1.fastq" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fastq"
    IS_GZ = 0
elif [ -f "${FASTQ_BASE_NAME}_1.fastq.gz" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fastq.gz"
    IS_GZ = 1
else
    errh "No FASTQ found!"
fi

if [ -f "${FASTQ_BASE_NAME}_2.fq.gz" ]; then
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fq.gz"
    IS_GZ = 1
elif [ -f "${FASTQ_BASE_NAME}_2.fastq" ]; then
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fastq"
    IS_GZ = 0
elif [ -f "${FASTQ_BASE_NAME}_2.fastq.gz" ]; then
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fastq.gz"
    IS_GZ = 1
elif [ -f "${FASTQ_BASE_NAME}_2.fq" ] ; then
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fq"
    IS_GZ = 0
else
    errh "No FASTQ found!"
fi

if [ -d "${GENE_REFERENCE}_STAR_idx" ]; then
    GENE_REFERENCE = "${GENE_REFERENCE}_STAR_idx"
fi

if [ ${IS_GZ} -eq 0 ]; then
    [ ! -f "${FASTQ_BASE_NAME}".GENE.bam ] && STAR --runThreadN "${THREAD_NUM}" --genomeDir "${GENE_REFERENCE}" \
    --readFilesIn "${FASTQ_NAME_1}" "${FASTQ_NAME_2}" --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outFileNamePrefix "${FASTQ_BASE_NAME}"
    | \ samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".GENE.bam
fi

if [ ${IS_GZ} -eq 1 ]; then
    [ ! -f "${FASTQ_BASE_NAME}".GENE.bam ] && STAR --runThreadN "${THREAD_NUM}" --genomeDir "${GENE_REFERENCE}" --readFilesCommand zcat \
    --readFilesIn "${FASTQ_NAME_1}" "${FASTQ_NAME_2}" --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outFileNamePrefix "${FASTQ_BASE_NAME}"
    | \ samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".GENE.bam
fi

[ ! -f "${FASTQ_NAME}".TRANS.bam ] && bwa mem -t 12 "${GENE_REFERENCE}" "${FASTQ_NAME_1}" "${FASTQ_NAME_2}" | \
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

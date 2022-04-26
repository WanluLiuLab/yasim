#!/usr/bin/env bash
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate yasim_c_elegans_as_depth_analysis
set -ueo pipefail

SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}"/shlib/libstr.sh
. "${SHDIR}"/shlib/libisopt.sh
. "${SHDIR}"/shlib/libdo.sh
. "${SHDIR}"/shlib/libprivate.sh

if [ -f "${FASTQ_BASE_NAME}_1.fq" ] && [ -f "${FASTQ_BASE_NAME}_2.fq" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fq"
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fq"
    IS_GZ=0
elif [ -f "${FASTQ_BASE_NAME}_1.fq.gz" ] && [ -f "${FASTQ_BASE_NAME}_2.fq.gz" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fq.gz"
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fq.gz"
    IS_GZ=1
elif [ -f "${FASTQ_BASE_NAME}_1.fastq" ] && [ -f "${FASTQ_BASE_NAME}_2.fastq" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fastq"
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fastq"
    IS_GZ=0
elif [ -f "${FASTQ_BASE_NAME}_1.fastq.gz" ]  && [ -f "${FASTQ_BASE_NAME}_2.fastq.gz" ]; then
    FASTQ_NAME_1="${FASTQ_BASE_NAME}_1.fastq.gz"
    FASTQ_NAME_2="${FASTQ_BASE_NAME}_2.fastq.gz"
    IS_GZ=1
else
    errh "No FASTQ found!"
fi

if [ -d "${GENE_REFERENCE}_STAR_idx" ]; then
    STAR_REFERENCE="${GENE_REFERENCE}_STAR_idx"
else
    errh "STAR index not found at ${GENE_REFERENCE}_STAR_idx"
fi

if [ -d "${TRANSCRIPT_REFERENCE}_bwa_idx" ]; then
    BWA_REFERENCE="${TRANSCRIPT_REFERENCE}_bwa_idx"
else
    errh "BWA index not found at ${TRANSCRIPT_REFERENCE}_bwa_idx"
fi

if [ ! -f "${FASTQ_BASE_NAME}".GENE.bam ];then
    if [ ${IS_GZ} -eq 0 ]; then
        readFilesCommand="cat"
    else
        readFilesCommand="zcat"
    fi
    STAR --runThreadN "${THREAD_NUM}" \
    --genomeDir "${STAR_REFERENCE}" \
    --readFilesCommand "${readFilesCommand}" \
    --readFilesIn "${FASTQ_NAME_1}" "${FASTQ_NAME_2}" \
    --outSAMtype BAM Unsorted \
    --outSAMattributes All \
    --outFileNamePrefix "${FASTQ_BASE_NAME}_STAR/"
    samtools sort "${FASTQ_BASE_NAME}_STAR/Aligned.out.bam" -@ "${THREAD_NUM}"  -o "${FASTQ_BASE_NAME}".GENE.bam
    rm -rf "${FASTQ_BASE_NAME}_STAR/"
fi


[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam ] && bwa mem -t "${THREAD_NUM}" "${BWA_REFERENCE}/bwa" "${FASTQ_NAME_1}" "${FASTQ_NAME_2}" | \
samtools sort - -@ "${THREAD_NUM}" -o "${FASTQ_BASE_NAME}".TRANS.bam

. "${SHDIR}"/shlib/libpost_alignment_analysis.sh
. "${SHDIR}"/shlib/libspladder_analysis_ngs.sh

#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

rm -rf dataset_sc_ngs
STAR --genomeDir "${STAR_INDEX}" \
    --runThreadN "${THREADS}" \
    --readFilesIn "${FASTQ_BASENAME}"_1.fq "${FASTQ_BASENAME}"_2.fq  \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix dataset_sc_ngs/
mv dataset_sc_ngs/Aligned.sortedByCoord.out.bam aligned_ngs.bam


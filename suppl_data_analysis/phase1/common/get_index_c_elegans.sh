#!/usr/bin/env bash
set -vueo pipefail

GENE_REFERENCE="ce11.fa"
GENE_GTF="ce11.ncbiRefSeq.gtf"
TRANSCRIPT_REFERENCE="ce11.trans.fa"

[ ! -f "${GENE_REFERENCE}" ] &&
    curl -L https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz | gzip -cdf - >"${GENE_REFERENCE}"
[ ! -f "${GENE_GTF}" ] &&
    curl -L https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz | gzip -cdf - >"${GENE_GTF}"

[ ! -f "${TRANSCRIPT_REFERENCE}" ] &&
    bedtools getfasta -nameOnly -s -fi "${GENE_REFERENCE}" -bed "${GENE_GTF}" >"${TRANSCRIPT_REFERENCE}"

STAR_REFERENCE="${GENE_REFERENCE}_STAR_idx"
[ ! -d "${STAR_REFERENCE}" ] && STAR --runThreadN 40 \
    --runMode genomeGenerate \
    --genomeDir "${STAR_REFERENCE}" \
    --genomeFastaFiles "${GENE_REFERENCE}" \
    --sjdbGTFfile "${GENE_GTF}"

BWA_REFERENCE="${TRANSCRIPT_REFERENCE}_bwa_idx"

if [ ! -d "${BWA_REFERENCE}" ]; then
    mkdir -p "${BWA_REFERENCE}"
    bwa index -p "${BWA_REFERENCE}/bwa" "${TRANSCRIPT_REFERENCE}" || rm -drf "${BWA_REFERENCE}"
fi

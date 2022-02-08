#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

rm "${FASTQ_BASENAME}".d -rf
python -m yasim dwgsim -F "${SEL_CDNA_FASTA}" -o "${FASTQ_BASENAME}" -d "${DEPTH_TSV}"

rm transcripts_quant/ -rf ; salmon quant -i "${SALMON_INDEX}" -l IU -1 "${FASTQ_BASENAME}"_1.fq -2 "${FASTQ_BASENAME}"_2.fq -o transcripts_quant

Rscript R/test_salmon.R

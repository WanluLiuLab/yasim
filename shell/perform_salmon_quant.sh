#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

rm transcripts_quant/ -rf ; salmon quant -i "${SALMON_INDEX}" --threads "${THREADS}" -l IU -1 "${FASTQ_BASENAME}"_1.fq -2 "${FASTQ_BASENAME}"_2.fq -o transcripts_quant

 Rscript R/test_salmon.R  --salmon_quant_sf transcripts_quant/quant.sf --yasim_tsv hg38.chr1.gene.depth.tsv --output yasim_to_salmon_quant

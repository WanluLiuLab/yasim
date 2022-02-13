#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

stringtie aligned_ngs.bam -G "${REFERENCE_GTF}" -o dataset_as_ngs.gtf -p "${THREADS}"

python src/scripts/parse_stringtie_into_tsv.py -g dataset_as_ngs.gtf -o stringtie_quant.tsv

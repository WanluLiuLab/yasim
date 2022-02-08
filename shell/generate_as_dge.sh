#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

[ -f "${SEL_GTF}" ] || python -m yasim alternative_splicing -f "${REFERENCE_FASTA}" -g "${REFERENCE_GTF}" -o "${SEL_GTF}"
[ -f "${SEL_CDNA_FASTA}" ] || python -m yasim transcribe -g "${SEL_GTF}" -o "${SEL_CDNA_FASTA}" -f "${REFERENCE_FASTA}"
[ -f "${DEPTH_TSV}" ] || python -m yasim dge -g "${SEL_GTF}" -o "${DEPTH_TSV}" -d 100 -l 20

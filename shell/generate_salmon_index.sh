#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

python -m yasim transcribe -g "${REFERENCE_GTF}" -o hg38.chr1.reference_transcripts.fa -f "${REFERENCE_FASTA}"
salmon index --transcripts hg38.chr1.reference_transcripts.fa --index "${SALMON_INDEX}" -p "${THREADS}"

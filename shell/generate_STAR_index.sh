#!/usr/bin/env bash
SHDIR="$(readlink -f "$(dirname "${0}")")"
. "${SHDIR}/conf.sh"

STAR --runThreadN "${THREADS}" \
--runMode genomeGenerate \
--genomeDir "${STAR_INDEX}" \
--genomeFastaFiles "${REFERENCE_FASTA}" \
--sjdbGTFfile "${REFERENCE_GTF}"

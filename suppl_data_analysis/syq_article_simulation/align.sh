#!/usr/bin/env bash
set -ue
echo "ALN ${1} -> ${2} START" >&2
FQ="${1}"
FA="${2}"
PARALLEL=66

minimap2 -ax splice --MD -t "${PARALLEL}" "${FA}" "${FQ}" > "${FQ}".sam
samtools sort -@ "${PARALLEL}" "${FQ}".sam -o "${FQ}".sam
samtools view -bS -@ "${PARALLEL}" "${FQ}".sam > "${FQ}".sam.bam
samtools sort -@ "${PARALLEL}" "${FQ}".sam.bam -o "${FQ}".sam.bam
samtools index -@ "${PARALLEL}" "${FQ}".sam.bam
echo "ALN ${1} -> ${2} FIN" >&2

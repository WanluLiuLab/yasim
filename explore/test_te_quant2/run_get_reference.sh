#!/usr/bin/env bash
# shellcheck disable=SC2002
set -ue
mkdir -p ref aln
# Get Reference
cd ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
cat ce11.ncbiRefSeq.gtf.gz | pigz -d >ce11.ncbiRefSeq.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
cat ce11.fa.gz | pigz -d >ce11.fa
samtools faidx ce11.fa

wget https://dfam.org/releases/current/families/Dfam_curatedonly.h5.gz






#!/usr/bin/env bash
set -ue
mkdir -p ref aln
# Get Reference
cd ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
cat ce11.ncbiRefSeq.gtf.gz | pigz -d >ce11.ncbiRefSeq.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
cat ce11.fa.gz | pigz -d >ce11.fa
samtools faidx ce11.fa
cd ..

perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir ref/ce11.rmsk_rmblast.d \
    -gff -small -s \
    ce11.fa




#!/usr/bin/env bash
perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir ref/ce11.rmsk_rmblast.d \
    -gff -small -s \
    ce11.fa

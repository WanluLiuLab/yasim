#!/usr/bin/env bash
set -e
cd "$(readlink -f "$(dirname "${0}")")"
cd ..

perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir ref/ce11.rmsk_rmblast.d \
    -gff -small -s \
    ref/ce11.fa

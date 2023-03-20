#!/usr/bin/env bash
set -uex
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz

gunzip ./*.gz

# Generate AS events with different complexity level (defaults 2)
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq.gtf

function generate_as_events(){
    python -m yasim generate_as_events \
        -f ce11.fa -g ce11.ncbiRefSeq.gtf \
        -o ce11_as_"${1}".gtf \
        --complexity "${1}"
    python -m labw_utils.bioutils describe_gtf ce11_as_"${1}".gtf
    python -m yasim transcribe \
        -g ce11_as_"${1}".gtf \
        -f ce11.fa \
        -o ce11_trans_"${1}".fa
}

generate_as_events 2

# Install LLRGs
mkdir -p bin
cd src || exit 1
bash run.sh
cd .. || exit 1

# Read accuracy
cd diff_accuracy_zipf4 || exit 1
bash run.sh
cd .. || exit 1

# Genome Complexity
cd diff_genome_complexity_zipf4 || exit 1
bash run.sh
cd .. || exit 1

# Read Completeness
cd diff_read_compl_zipf4 || exit 1
bash run.sh
cd .. || exit 1

# Read depth
cd diff_read_depth_ng_zipf33 || exit 1
bash run.sh
cd .. || exit 1

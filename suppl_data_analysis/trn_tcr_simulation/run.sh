#!/usr/bin/env bash
set -ue

axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz

gunzip ./*.gz

export PYTHONPATH="../../src:../../deps/labw_utils/src:${PYTHONPATH}"

# Target: Cell 1000

python -m yasim_scripts exclude_tcr_from_gtf hg38.ncbiRefSeq.gtf
python -m yasim generate_gene_depth \
    -g hg38.ncbiRefSeq.gtf.tcr_filtered.gtf \
    -d 10 \
    -o notcr_gene_depth.tsv
python -m yasim generate_isoform_depth \
    -g hg38.ncbiRefSeq.gtf.tcr_filtered.gtf \
    -d notcr_gene_depth.tsv \
    -o notcr_isoform_depth.tsv
python -m yasim generate_isoform_replicates \
    -d notcr_isoform_depth.tsv \
    -n 100

python -m yasim transcribe \
    -f hg38.fa \
    -g hg38.ncbiRefSeq.gtf.tcr_filtered.gtf \
    -o hg38_tcr_filtered_trans.fa

python -m yasim dwgsim \
    -F hg38_tcr_filtered_trans.fa.d \
    -o ref \
    -d notcr_isoform_depth.tsv \
    -e $(which dwgsim) \
    -j 20

python create_tcr_cache.py
python rearrange_tcr.py

art_illumina -f 10 --in sim.fa -ss HS25 -l 150 --out sim

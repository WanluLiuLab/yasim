#!/usr/bin/env bash
set -ue

axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz

gunzip ./*.gz

export PYTHONPATH="../../src:../../deps/labw_utils/src:${PYTHONPATH}"

grep -e '^chr7\s' -e '^chr14\s' hg38.ncbiRefSeq.gtf > hg38.ncbiRefSeq_chr7_14.gtf

# TCR
python generate_tcr_depth.py
python -m yasim_sc generate_barcode -n 10 -o barcode.txt
python create_tcr_cache.py
python rearrange_tcr.py
python -m labw_utils.bioutils split_fasta sim_tcr.fa
python -m yasim art \
    -F sim_tcr.fa.d \
    -o sim_tcr  \
    -d tcr_depth.tsv \
    -e art_illumina\
    -j 20

# Gene
python -m yasim generate_gene_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d 5 \
    -o notcr_gene_depth.tsv
python -m yasim generate_isoform_depth \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -d notcr_gene_depth.tsv \
    -o notcr_isoform_depth.tsv
python -m yasim_sc generate_barcoded_isoform_replicates \
    -d notcr_isoform_depth.tsv \
    -b barcode.txt

python -m yasim transcribe \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_subsampled.gtf \
    -o notcr_trans.fa

for fn in notcr_isoform_depth.tsv.d/*.tsv; do
    python -m yasim art \
        -F notcr_trans.fa.d \
        -o "${fn}"_sim \
        -d "${fn}" \
        -e art_illumina \
        -j 20
done

mkdir -p trust4_result

run-trust4 -u sim_tcr.fq \
    -f trust4_index/bcrtcr.fa \
    -t 40 \
    --ref trust4_index/IMGT+C.fa \
    --od trust4_result

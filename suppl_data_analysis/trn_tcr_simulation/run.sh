#!/usr/bin/env bash
set -ue

axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz

gunzip ./*.gz

export PYTHONPATH="../../src:../../deps/labw_utils/src:${PYTHONPATH}"

# Target: Cell 1000
python -m yasim_sc generate_barcode -n 10 -o barcode.txt
python create_tcr_cache.py
python rearrange_tcr.py
python sequencing_art_illumina.py

python -m yasim_scripts exclude_tcr_from_gtf hg38.ncbiRefSeq.gtf
cat hg38.ncbiRefSeq.gtf | grep '^chr1\s' > hg38.ncbiRefSeq_chr1.gtf
python -m yasim generate_gene_depth \
    -g hg38.ncbiRefSeq_chr1.gtf \
    -d 5 \
    -o notcr_gene_depth.tsv
python -m yasim generate_isoform_depth \
    -g hg38.ncbiRefSeq_chr1.gtf \
    -d notcr_gene_depth.tsv \
    -o notcr_isoform_depth.tsv
python -m yasim_sc generate_barcoded_isoform_replicates \
    -d notcr_isoform_depth.tsv \
    -b barcode.txt

python -m yasim transcribe \
    -f hg38.fa \
    -g hg38.ncbiRefSeq_chr1.gtf \
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

run-trust4 -u sim.fq \
    -f trust4_index/bcrtcr.fa \
    -t 40 \
    --ref trust4_index/IMGT+C.fa \
    --od trust4_result

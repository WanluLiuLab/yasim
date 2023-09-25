#!/usr/bin/env bash
set -e
conda activate yasim_dev

mkdir -p sim
python -m yasim generate_gene_te_fusion \
    -g ref/ce11.ncbiRefSeq.gtf \
    -f ref/ce11.fa \
    --tedb ref/teidx.pkl.xz \
    -o sim/ce11_denovo_test \
    -d 150 \
    -n 100000
sbatch <slrum/run_aln_gtts.sh

python -m yasim generate_gene_te_fusion \
    -g ref/ce11.ncbiRefSeq.gtf \
    -f ref/ce11.fa \
    --tedb ref/teidx.pkl.xz \
    -o sim/ce11_denovo_test2 \
    -d 150 \
    -n 10

cat sim/ce11_denovo_test_split/*.nhmmer_original.out > aln/nhmmer_original.out
python run_convert_aln.py
xz -9 -T5 -f -vv all.cmp.tsv

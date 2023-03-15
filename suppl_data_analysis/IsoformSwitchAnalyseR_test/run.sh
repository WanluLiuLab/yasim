#!/usr/bin/env bash

# Use Ensembl C. Elegans MT as reference
python -m yasim generate_gene_depth -g fake.gtf -o fake_gene_depth.tsv
python -m yasim generate_isoform_depth -g fake.gtf -d fake_gene_depth.tsv -o fake_isoform_depth_1.tsv --alpha 4
python -m yasim generate_isoform_depth -g fake.gtf -d fake_gene_depth.tsv -o fake_isoform_depth_2.tsv --alpha 4
python -m labw_utils.bioutils transcribe -g fake.gtf -f fake.fa -o fake_trans.fa
Rscript perform.R

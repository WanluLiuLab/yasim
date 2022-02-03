#!/usr/bin/env bash
set -ue
# PYTHONPATH=src python -m yasim alternative_splicing -f /home/yuzj/ref/Ensembl/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa -g /home/yuzj/ref/Ensembl/Homo_sapiens.GRCh38.105_chr1.gtf -o Homo_sapiens.GRCh38.105_chr1_sel.gtf
# PYTHONPATH=src python -m yasim transcribe -g Homo_sapiens.GRCh38.105_chr1_sel.gtf -o "hg38.chr1.transcripts.fa" -f /home/yuzj/ref/Ensembl/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa
samtools faidx hg38.chr1.transcripts.fa
PYTHONPATH=src python -m yasim dge -g Homo_sapiens.GRCh38.105_chr1_sel.gtf -o "hg38.chr1.gene.depth.tsv" -d 1000 -l 100
PYTHONPATH=src python -m yasim dwgsim -F "hg38.chr1.transcripts.fa" -o "hg38.chr1.transcripts.fq" -d "hg38.chr1.gene.depth.tsv"

rm transcripts_quant/ -rf ; salmon quant -i /home/yuzj/ref/Ensembl/salmon_index -l IU -1 hg38.chr1.transcripts.fq_1.fq.gz -2 hg38.chr1.transcripts.fq_2.fq.gz -o transcripts_quant



#!/usr/bin/env bash

REF=/home/tgs/Refdata/ensembl/fa/Drosophila_melanogaster.BDGP6.32.dna.toplevel.trans.fa
for fn in drosophila/*/*.fastq; do
    if [ ! -e "${fn}".trans.depth.filtered.parquet ]; then
        minimap2 -t 20 -a $REF "${fn}" 2> /dev/null | \
            samtools sort --write-index -@ 20 -o "${fn}".trans.bam &>> /dev/null
        samtools depth -@ 20 "${fn}".trans.bam | xz -9 -T20 > "${fn}".trans.depth.tsv.xz
        echo "${fn}"
    fi
    Rscript depth_aln_to_transcriptome.R \
        --input "${fn}".trans.depth.tsv.xz \
        --trans_stats "${REF}.stats" \
        --output "${fn}".trans.depth.filtered.parquet
    rm -f "${fn}".trans.bam
done

REF=/home/tgs/Refdata/ensembl/fa/Mus_musculus.GRCm39.dna.primary_assembly.trans.fa
for fn in mouse/*/*.fastq; do
    if [ ! -e "${fn}".trans.depth.filtered.parquet ]; then
        minimap2 -t 20 -a $REF "${fn}" 2> /dev/null | \
            samtools sort --write-index -@ 20 -o "${fn}".trans.bam &>> /dev/null
        samtools depth -@ 20 "${fn}".trans.bam | xz -9 -T20 > "${fn}".trans.depth.tsv.xz
        echo "${fn}"
    fi
    Rscript depth_aln_to_transcriptome.R \
        --input "${fn}".trans.depth.tsv.xz \
        --trans_stats "${REF}.stats" \
        --output "${fn}".trans.depth.filtered.parquet
    rm -f "${fn}".trans.bam
done
REF=/home/tgs/Refdata/ensembl/fa/Homo_sapiens.GRCh38.dna.primary_assembly.trans.fa
for fn in human/*/*.fastq; do
    if [ ! -e "${fn}".trans.depth.filtered.parquet ]; then
        minimap2 -t 20 -a $REF "${fn}" 2> /dev/null | \
            samtools sort --write-index -@ 20 -o "${fn}".trans.bam &>> /dev/null
        samtools depth -@ 20 "${fn}".trans.bam | xz -9 -T20 > "${fn}".trans.depth.tsv.xz
        echo "${fn}"
    fi
    Rscript depth_aln_to_transcriptome.R \
        --input "${fn}".trans.depth.tsv.xz \
        --trans_stats "${REF}.stats" \
        --output "${fn}".trans.depth.filtered.parquet
    rm -f "${fn}".trans.bam
done

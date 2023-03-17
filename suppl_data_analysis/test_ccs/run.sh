#!/usr/bin/env bash
set -uex
export NUM_THREADS=40

axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz


grep -i '^chrI\s' ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" > ce11.chr1.fa
samtools faidx ce11.chr1.fa

python -m labw_utils.bioutils transcribe \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -f ce11.chr1.fa \
    -o ce11_trans.chr1.fa
python -m yasim generate_as_events \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq_as.chr1.gtf \
    -c 2
python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_trans_as.chr1.fa
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq_as.chr1.gtf
python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_gene_depth.chr1.tsv \
    -d 20
python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_isoform_depth.chr1.tsv \
    -d ce11_gene_depth.chr1.tsv \
    --alpha 4
python -m yasim pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -F ce11_trans_as.chr1.fa.d \
    -d ce11_isoform_depth.chr1.tsv \
    -o ce11_ccs \
    -j 80 \
    --ccs_pass 10

# pbmerge -o out.ccs.bam ce11_ccs.d/*/tmp.ccs.bam
python merge.py out.ccs.bam ce11_ccs.d/*/tmp.ccs.bam
pbindex out.ccs.bam
samtools index out.ccs.bam

isoseq3 cluster \
    out.ccs.bam \
    out.transcripts.xml \
    --log-level INFO \
    --log-file isoseq3_cluster.log \
    --num-threads 40

pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --log-file aln.log \
    --log-level INFO \
    out.transcripts.xml.hq.bam \
    ce11.chr1.fa \
    aln.bam

isoseq3 collapse \
    --do-not-collapse-extra-5exons \
    --log-level INFO \
    --log-file isoseq3_collapse.log \
    aln.bam out.ccs.bam collapse.gff
Rscript get_express_gtf.R -s ce11_ccs.fq.stats -g ce11.ncbiRefSeq_as.chr1.gtf -o ce11.ncbiRefSeq_as.chr1.expressed.gtf
gffcompare collapse.gff -r  ce11.ncbiRefSeq_as.chr1.expressed.gtf -o gffcmp
pigeon sort ce11.ncbiRefSeq.chr1.gtf -o ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf
pigeon index ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf

pigeon classify collapse.gff ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf ce11.chr1.fa -o classified
pigeon filter classified_classification.txt --isoforms collapse.gff

#!/usr/bin/env bash
set -ue
mkdir -p CE11_GENOME_STAR_INDEX/
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir CE11_GENOME_STAR_INDEX --genomeFastaFiles ce11.fa --sjdbGTFfile ce11.ncbiRefSeq.gtf
minimap2 -d ce11.fa.mmi ce11.fa -t 50

minimap2 -x splice -a -t 50 ce11.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz |
    samtools sort -@ 50 -o NANO_GENE.bam # 85.55% mapping rate
minimap2 -x splice -a -t 50 ce11.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz |
    samtools sort -@ 50 -o PACB_GENE.bam # 97.27% mapping rate

STAR --genomeDir CE11_GENOME_STAR_INDEX \
    --runThreadN 50 \
    --readFilesCommand zcat \
    --readFilesIn \
    SRR5123644_1_c.fastq.gz,SRR5123648_1_c.fastq.gz,SRR5123649_1_c.fastq.gz \
    SRR5123644_2_c.fastq.gz,SRR5123648_2_c.fastq.gz,SRR5123649_2_c.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ILLM_GENE

for fn in *_GENE.bam; do samtools index "${fn}" & done
wait

spladder build -b ILLM_GENE.bam -o TEST_ILLM/ -a ce11.ncbiRefSeq.gtf --parallel 50 --ignore-mismatche

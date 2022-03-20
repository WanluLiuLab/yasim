#!/usr/bin/env bash
set -ue
mkdir -p CE11_TRANSCRIPTOME_BWA_INDEX
bwa index ce11.reference_transcripts.fa -p CE11_TRANSCRIPTOME_BWA_INDEX/ce11
minimap2 -d ce11.reference_transcripts.fa.mmi ce11.reference_transcripts.fa -t 50

minimap2 -a -t 50 ce11.reference_transcripts.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz | \
samtools sort -@ 50 -o NANO_TRANS.bam # 84.24% mapping rate
minimap2 -a -t 50 ce11.reference_transcripts.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz | \
samtools sort -@ 50 -o PACB_TRANS.bam # 96.65% mapping rate

bwa mem -t 64 CE11_TRANSCRIPTOME_BWA_INDEX/ce11 \
<(zcat SRR5123644_1_c.fastq.gz SRR5123648_1_c.fastq.gz SRR5123649_1_c.fastq.gz) \
<(zcat SRR5123644_2_c.fastq.gz SRR5123648_2_c.fastq.gz SRR5123649_2_c.fastq.gz) |\
samtools sort -@ 50 -o ILLM_TRANS.bam # 78.37% mapping rate

for fn in *_TRANS.bam; do samtools index "${fn}"; samtools depth "${fn}" > "${fn}".depth.tsv & done
wait

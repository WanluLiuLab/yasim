#!/usr/bin/env bash
cat << EOF | while read line; do [ ! -f "$(basename "${line}")" ] && axel -n 20 "${line}"; done
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/001/ERR3245471/ERR3245471.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/000/ERR3245470/ERR3245470.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8568877/SRR8568877_subreads.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/008/SRR8568878/SRR8568878_subreads.fastq.gz
EOF

bwa index ce11.reference_transcripts.fa -p CE11_TRANSCRIPTOME_BWA_INDEX
minimap2 -d ce11.reference_transcripts.fa.mmi ce11.reference_transcripts.fa -t 16

minimap2 -x splice -a -t 16 ce11.reference_transcripts.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz | \
samtools sort -@ 16 -o nanopore_transcriptome_1.bam # 84.24% mapping rate
minimap2 -x splice -a -t 16 ce11.reference_transcripts.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz | \
samtools sort -@ 16 -o pacbio_transcriptome_1.bam # 96.65% mapping rate


# SRR5123648_1.fastq.gz SRR5123649_1.fastq.gz
# SRR5123648_2.fastq.gz SRR5123649_2.fastq.gz
bwa mem -t 16 CE11_TRANSCRIPTOME_BWA_INDEX <(zcat SRR5123644_1.fastq.gz) <(zcat SRR5123644_2.fastq.gz) |\
samtools sort -@ 16 -o ngs_transcript.bam # 70.35%, 74.77%, 88.43% respectively

for fn in *.bam; do samtools index $fn; samtools depth $fn > $fn.depth.tsv; done

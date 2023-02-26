#!/usr/bin/env bash
set -ue
cat << EOF | while read -r line; do [ ! -f "$(basename "${line}")" ] && axel -n 20 "${line}"; done
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
https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
EOF

gzip -dk ce11.fa.gz ce11.ncbiRefSeq.gtf.gz
samtools faidx ce11.fa
gffread ce11.ncbiRefSeq.gtf -g ce11.fa -w ce11.reference_transcripts.fa

fastqc --extract SRR512364?_?.fastq.gz -t 50

cutadapt -j 10 -u 10 -u -20 -o SRR5123644_1_c.fastq.gz SRR5123644_1.fastq.gz &
cutadapt -j 10 -u 10 -u -20 -o SRR5123644_2_c.fastq.gz SRR5123644_2.fastq.gz &
cutadapt -j 10 -u 10 -u -10 -o SRR5123648_1_c.fastq.gz SRR5123648_1.fastq.gz &
cutadapt -j 10 -u 10 -u -10 -o SRR5123648_2_c.fastq.gz SRR5123648_2.fastq.gz &
cutadapt -j 10 -u 10 -u -30 -o SRR5123649_1_c.fastq.gz SRR5123649_1.fastq.gz &
cutadapt -j 10 -u 10 -u -30 -o SRR5123649_2_c.fastq.gz SRR5123649_2.fastq.gz &
wait

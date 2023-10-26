#!/usr/bin/env bash
set -ue

mkdir -p ref
cd ref
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes
cd ..

mkdir -p raw_data
prefetch SRR9286043 SRR9286044 --output-directory raw_data
cd raw_data
fasterq-dump --outdir . --threads 32 --progress --split-files SRR9286043/SRR9286043.sra
fasterq-dump --outdir . --threads 32 --progress --split-files SRR9286044/SRR9286044.sra
rm -fr SRR9286044 SRR9286043
fastqc -o . -t 36 SRR9286043_1.fastq SRR9286043_2.fastq SRR9286044_1.fastq SRR9286044_2.fastq
cd ..
for acc in SRR9286044 SRR9286043; do
    fastp \
        --in1 raw_data/"${acc}"_1.fastq \
        --out1 fastp_tmp/"${acc}"_1.fastq \
        --in2 raw_data/"${acc}"_2.fastq \
        --out2 fastp_tmp/"${acc}"_2.fastq \
        --thread 16 \
        --json fastp_tmp/"${acc}".fastp.json \
        --html fastp_tmp/"${acc}".fastp.html \
        --verbose \
        --trim_front1 10 \
        --trim_front2 10 &
done
wait
bwa-mem2 index ref/ce11.fa
mkdir -p aln
for acc in SRR9286044 SRR9286043; do
        bwa-mem2 mem \
            -t 32 -S -P \
            ref/ce11.fa \
            fastp_tmp/"${acc}"_1.fastq \
            fastp_tmp/"${acc}"_2.fastq | \
        samtools view -h -@ 36 -o aln/"${acc}".bam
    pairtools parse \
        -o aln/"${acc}".pair.gz \
        -c ref/ce11.chrom.sizes \
        --drop-sam \
        --drop-seq \
        --output-stats aln/"${acc}".pair.gz.stats \
        --assembly ce11 \
        --no-flip \
        --add-columns mapq \
        --walks-policy mask \
        aln/"${acc}".bam
    pairtools sort \
        --nproc 36 \
        -o aln/"${acc}".sorted.pair.gz \
        aln/"${acc}".pair.gz
    pairtools dedup \
        --max-mismatch 3 \
        --mark-dups \
        --output \
            >( pairtools split \
                --output-pairs aln/"${acc}".nodups.pairs.gz \
                --output-sam aln/"${acc}".nodups.bam \
             ) \
        --output-stats aln/"${acc}".dedup.stats \
        aln/"${acc}".pairs.gz
    pairtools select \
        "mapq1>0 and mapq2>0" \
        aln/"${acc}".nodups.pairs.gz \
        -o aln/"${acc}".nodups.UU.pairs.gz

done

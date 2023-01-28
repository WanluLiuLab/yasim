#!/usr/bin/env bash
export NUM_THREADS=40

mkdir -p sim
rm -rf sim/*

pbsim3 \
    --strategy wgs \
    --method qshmm \
    --qshmm ../../src/yasim/llrg_adapter/pbsim3_dist/QSHMM-ONT.model \
    --prefix sim/ce11 \
    --genome ce11.fa  &> pbsim.log

cat sim/*.fastq | pigz -9  > simulated.fastq.gz

mkdir -p lastdb

lastdb -v -P"${NUM_THREADS}" lastdb/ce11 ce11.fa
zcat simulated.fastq.gz | last-train -v -Qfastx -P"${NUM_THREADS}" lastdb/ce11 /dev/stdin > simulated.train
zcat simulated.fastq.gz | lastal -v -P"${NUM_THREADS}" -p simulated.train -Qfastx lastdb/ce11 /dev/stdin | pigz -9 > simulated.maf.gz

cat sim/*.maf | pigz -9 > simulated_gt.maf.gz

#!/usr/bin/env bash
export NUM_THREADS=40

mkdir -p sim
rm -rf sim/*
python filter_pbsim_reference.py < ce11_trans.fa | head -n 1000 > ce11_trans.filtered.fa

pbsim --prefix sim/sim \
    --data-type CLR \
    --depth 50 \
    --model_qc ../../src/yasim/llrg_adapter/pbsim_dist/model_qc_clr \
    ce11_trans.filtered.fa  &> pbsim.log

cat sim/*.fastq > simulated.fastq

mkdir -p lastdb

lastdb -P"${NUM_THREADS}" -uRY4 lastdb/ce11 ce11.fa
last-train -P"${NUM_THREADS}" -Q0 lastdb/ce11 simulated.fastq > simulated.train
lastal -P"${NUM_THREADS}" --split -p simulated.train lastdb/ce11 simulated.fastq > simulated.maf

lastdb -P"${NUM_THREADS}" -uRY4 lastdb/ce11_trans.filtered ce11_trans.filtered.fa
last-train -P"${NUM_THREADS}" -Q0 lastdb/ce11_trans.filtered simulated.fastq > simulated_trans.train
lastal -P"${NUM_THREADS}" -p simulated_trans.train lastdb/ce11_trans.filtered simulated.fastq > simulated_trans.maf



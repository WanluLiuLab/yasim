#!/usr/bin/env bash
set -ue
python filter_pbsim_reference.py < ce11_trans.fa | head -n 4000 > ce11_trans.filtered.fa

function perform_pbsim() {
    mkdir -p pbsim_"${1}"
    pbsim \
        --prefix pbsim_"${1}"/ce11 \
        --data-type CLR \
        --depth 0.5 \
        --length-mean "${1}" \
        --model_qc ../../src/yasim/llrg_adapter/pbsim_dist/model_qc_clr \
        ce11_trans.filtered.fa &> pbsim_"${1}".log
    cat pbsim_"${1}"/ce11_*.fastq > pbsim_"${1}".fastq
    rm -rf pbsim_"${1}"
    python -m labw_utils.bioutils describe_fastq pbsim_"${1}".fastq
}

function perform_pbsim2() {
    mkdir -p pbsim2_"${1}"
    pbsim2 \
        --prefix pbsim2_"${1}"/ce11 \
        --depth 0.5 \
        --length-mean "${1}" \
        --hmm_model ../../src/yasim/llrg_adapter/pbsim2_dist/R103.model \
        ce11_trans.filtered.fa &> pbsim2_"${1}".log
    cat pbsim2_"${1}"/ce11_*.fastq > pbsim2_"${1}".fastq
    rm -rf pbsim2_"${1}"
    python -m labw_utils.bioutils describe_fastq pbsim2_"${1}".fastq
}

for i in 400 800 1200 1600 2000 2400 2800 3200; do
    perform_pbsim "${i}" &
    perform_pbsim2 "${i}" &
done
wait

#!/usr/bin/env bash
set -ue
ln -s ../diff_genome_complexity_zipf4/ce11_as_3.gtf .
ln -s ../diff_genome_complexity_zipf4/ce11_as_3.gtf.*.tsv .
ln -s ../diff_genome_complexity_zipf4/ce11_trans_3.fa .
ln -s ../diff_genome_complexity_zipf4/ce11_trans_3.fa.d .

function generate_depth() {
    LOG_FILE_NAME="yasim_generate_gene_depth_${1}.log" \
        python -m yasim generate_gene_depth \
        -g ce11_as_3.gtf \
        -o ce11_as_3_gene_depth_"${1}".tsv \
        -d "${1}" \
        --low_cutoff 0
    LOG_FILE_NAME="yasim_generate_isoform_depth_${1}.log" \
        python -m yasim generate_isoform_depth \
        -g ce11_as_3.gtf \
        -d ce11_as_3_gene_depth_"${1}".tsv \
        -o ce11_as_3_isoform_depth_"${1}".tsv \
        --alpha 3 \
        --low_cutoff 0
}

function perform_housekeeping() {
    rm -rf "${1}".d && \
    touch "${1}".finished || return 1
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ce11_trans_3.fa.d \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_RSII_CCS_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim3_RSII_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ce11_trans_3.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CLR_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim3_SEQUEL_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ce11_trans_3.fa.d \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CCS_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim3_SEQUEL_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ce11_trans_3.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim2_"${2}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${2}" \
        -F ce11_trans_3.fa.d \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_simulation() {
    perform_pbsim3_RSII_CLR_simulation "${1}" &
    perform_pbsim3_RSII_CCS_simulation "${1}" &
    perform_pbsim3_SEQUEL_CLR_simulation "${1}" &
    perform_pbsim3_SEQUEL_CCS_simulation "${1}" &
    for pbsim2_mode in R94 R103; do
        perform_pbsim2_simulation "${1}" "${pbsim2_mode}" &
    done
}

 generate_depth 100

 for mean_depth in 10 25 40 55 70; do
     python downscale_depth.py \
         ce11_as_3_isoform_depth_100.tsv \
         ce11_as_3_isoform_depth_"${mean_depth}".tsv \
         "0.${mean_depth}" &
 done
 wait

for mean_depth in 10 25 40 55 70; do
    perform_simulation ce11_as_3_isoform_depth_"${mean_depth}"
done
wait


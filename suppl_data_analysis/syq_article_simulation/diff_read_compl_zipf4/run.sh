#!/usr/bin/env bash
set -uex
LOG_FILE_NAME="yasim_generate_gene_depth.log" \
    python -m yasim generate_gene_depth \
    -g ../ce11_as_2.gtf \
    -o ce11_as_2_gene_depth_20.tsv \
    -d 20
LOG_FILE_NAME="yasim_generate_isoform_depth.log" \
    python -m yasim generate_isoform_depth \
    -g ../ce11_as_2.gtf \
    -d ce11_as_2_gene_depth_20.tsv \
    -o ce11_as_2_isoform_depth_20.tsv \
    --alpha 4

function perform_housekeeping() {
    rm -rf "${1}".d &&
        touch "${1}".finished || return 1
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_"${1}"_"${2}"_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --truncate_ratio_3p "${1}" \
        --truncate_ratio_5p "${2}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_RSII_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_"${1}"_"${2}"_pbsim3_RSII_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --truncate_ratio_3p "${1}" \
        --truncate_ratio_5p "${2}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_SEQUEL_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_"${1}"_"${2}"_pbsim3_SEQUEL_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --truncate_ratio_3p "${1}" \
        --truncate_ratio_5p "${2}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_SEQUEL_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_"${1}"_"${2}"_pbsim3_SEQUEL_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --truncate_ratio_3p "${1}" \
        --truncate_ratio_5p "${2}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_"${1}"_"${2}"_pbsim2_"${3}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${3}" \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --truncate_ratio_3p "${1}" \
        --truncate_ratio_5p "${2}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_simulation() {
    perform_pbsim3_RSII_CLR_simulation "${1}" "${2}" &
    perform_pbsim3_RSII_CCS_simulation "${1}" "${2}" &
    perform_pbsim3_SEQUEL_CCS_simulation "${1}" "${2}" &
    perform_pbsim3_SEQUEL_CLR_simulation "${1}" "${2}" &
    for pbsim2_mode in R94 R103; do
        perform_pbsim2_simulation "${1}" "${2}" "${pbsim2_mode}" &
    done
    wait
}

perform_simulation 0.0 0.4
perform_simulation 0.2 0.2
perform_simulation 0.4 0.0

perform_simulation 0.0 0.2
perform_simulation 0.1 0.1
perform_simulation 0.2 0.0

perform_simulation 0.0 0.0

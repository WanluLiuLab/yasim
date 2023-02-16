#!/usr/bin/env bash
set -uex

function generate_depth() {
    LOG_FILE_NAME="yasim_generate_gene_depth_${1}.log" \
        python -m yasim generate_gene_depth \
        -g ../ce11_as_2.gtf.gz \
        -o ce11_as_2_gene_depth_"${1}".tsv.xz \
        -d "${1}"
    LOG_FILE_NAME="yasim_generate_isoform_depth_${1}.log" \
        python -m yasim generate_isoform_depth \
        -g ../ce11_as_2.gtf.gz \
        -d ce11_as_2_gene_depth_"${1}".tsv.xz \
        -o ce11_as_2_isoform_depth_"${1}".tsv.xz
}

function perform_housekeeping() {
    gzip -9 "${1}".fq &&
        cat "${1}".d/*/*.maf | gzip -9 >"${1}".maf.gz &&
        rm -rf "${1}".d
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME="${2}"_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_RSII_CCS_simulation() {
    OUTPUT_BASENAME="${2}"_pbsim3_RSII_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CLR_simulation() {
    OUTPUT_BASENAME="${2}"_pbsim3_SEQUEL_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CCS_simulation() {
    OUTPUT_BASENAME="${2}"_pbsim3_SEQUEL_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME="${2}"_pbsim2_"${3}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${3}" \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 \
        --accuracy-mean "${1}" || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_simulation() {
    perform_pbsim3_RSII_CLR_simulation "${1}" "${2}" &
    perform_pbsim3_RSII_CCS_simulation "${1}" "${2}" &
    perform_pbsim3_SEQUEL_CLR_simulation "${1}" "${2}" &
    perform_pbsim3_SEQUEL_CCS_simulation "${1}" "${2}" &
    for pbsim2_mode in R94 R103; do
        perform_pbsim2_simulation "${1}" "${2}" "${pbsim2_mode}" &
    done
    wait
}

for mean_depth in 10 25 40 55 70 85 100; do
    generate_depth "${mean_depth}" &
done
wait

for mean_depth in 10 25 40 55 70 85 100; do
    perform_simulation ce11_as_2_isoform_depth_"${mean_depth}" ce11_as_2_isoform_depth_"${mean_depth}"
done

xz ./*.stats -vv -9 -T0

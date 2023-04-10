#!/usr/bin/env bash
set -ue

function generate_as_events() {
    LOG_FILE_NAME="yasim_generate_as_events_${1}.log" \
        python -m yasim generate_as_events \
        -f ../ce11.fa \
        -g ../ce11.ncbiRefSeq.gtf \
        -o ce11_as_"${1}".gtf \
        --complexity "${1}"
    python -m labw_utils.bioutils describe_gtf ce11_as_"${1}".gtf

    LOG_FILE_NAME="yasim_transcribe_${1}.log" \
        python -m yasim transcribe \
        -g ce11_as_"${1}".gtf \
        -f ../ce11.fa \
        -o ce11_trans_"${1}".fa
    LOG_FILE_NAME="yasim_generate_gene_depth_${1}.log" \
        python -m yasim generate_gene_depth \
        -g ce11_as_"${1}".gtf \
        -o ce11_as_"${1}"_gene_depth_20.tsv \
        -d 20
    LOG_FILE_NAME="yasim_generate_isoform_depth_${1}.log" \
        python -m yasim generate_isoform_depth \
        -g ce11_as_"${1}".gtf \
        -d ce11_as_"${1}"_gene_depth_20.tsv \
        -o ce11_as_"${1}"_isoform_depth_20.tsv \
        --alpha 4
}

function perform_housekeeping() {
    rm -rf "${1}".d &&
        touch "${1}".finished || return 1
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ce11_trans_"${1}".fa.d \
        -d ce11_as_"${1}"_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_RSII_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim3_RSII_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ce11_trans_"${1}".fa.d \
        --ccs_pass 10 \
        -d ce11_as_"${1}"_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim3_SEQUEL_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ce11_trans_"${1}".fa.d \
        -d ce11_as_"${1}"_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim3_SEQUEL_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim3_SEQUEL_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ce11_trans_"${1}".fa.d \
        --ccs_pass 10 \
        -d ce11_as_"${1}"_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim2_"${2}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${2}" \
        -F ce11_trans_"${1}".fa.d \
        -d ce11_as_"${1}"_isoform_depth_20.tsv \
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
    wait
}

for gene_complexity_level in 1 3 5 7 9; do
    generate_as_events "${gene_complexity_level}" &
done
wait

for gene_complexity_level in 1 3 5 7 9; do
    perform_simulation "${gene_complexity_level}"
done

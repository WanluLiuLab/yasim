#!/usr/bin/env bash
set -uex

function generate_as_events(){
    # LOG_FILE_NAME="yasim_generate_as_events_${1}.log" \
    #     python -m yasim generate_as_events \
    #     -f ../ce11.fa.gz \
    #     -g ../ce11.ncbiRefSeq.gtf.gz \
    #     -o ce11_as_"${1}".gtf.gz \
    #     --complexity "${1}"
    python -m labw_utils.bioutils describe_gtf ce11_as_"${1}".gtf.gz
    xz -9 -vv -f ce11_as_"${1}".gtf.gz.*.tsv
    LOG_FILE_NAME="yasim_transcribe_${1}.log" \
        python -m yasim transcribe \
        -g ce11_as_"${1}".gtf.gz \
        -f ../ce11.fa.gz \
        -o ce11_trans_"${1}".fa
    LOG_FILE_NAME="yasim_generate_gene_depth_${1}.log" \
        python -m yasim generate_gene_depth \
        -g ce11_as_"${1}".gtf.gz \
        -o ce11_as_"${1}"_gene_depth_20.tsv.xz \
        -d 20
    LOG_FILE_NAME="yasim_generate_isoform_depth_${1}.log" \
        python -m yasim generate_isoform_depth \
        -g ce11_as_"${1}".gtf.gz \
        -d ce11_as_"${1}"_gene_depth_20.tsv.xz \
        -o ce11_as_"${1}"_isoform_depth_20.tsv.xz \
        --alpha 4
}

function perform_housekeeping() {
    gzip -9f "${1}".fq &&
        cat "${1}".d/*/*.maf | gzip -9f >"${1}".maf.gz &&
        rm -rf "${1}".d
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_"${1}"_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ce11_trans_"${1}".fa.d \
        -d ce11_as_"${1}"_isoform_depth_20.tsv.xz \
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
        -d ce11_as_"${1}"_isoform_depth_20.tsv.xz \
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
        -d ce11_as_"${1}"_isoform_depth_20.tsv.xz \
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
        -d ce11_as_"${1}"_isoform_depth_20.tsv.xz \
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
        -d ce11_as_"${1}"_isoform_depth_20.tsv.xz \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_simulation() {
    perform_pbsim3_RSII_CLR_simulation "${1}"  &
    perform_pbsim3_RSII_CCS_simulation "${1}" &
    perform_pbsim3_SEQUEL_CLR_simulation "${1}" &
    perform_pbsim3_SEQUEL_CCS_simulation "${1}" &
    for pbsim2_mode in R103; do
        perform_pbsim2_simulation "${1}" "${pbsim2_mode}" &
    done
    wait
}

# for gene_complexity_level in 9; do
#     generate_as_events "${gene_complexity_level}" &
# done
# wait

for gene_complexity_level in 1 3 5 7; do
    perform_simulation "${gene_complexity_level}"
done

xz -vv -9 -T0 -f ./*.stats

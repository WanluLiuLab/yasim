#!/usr/bin/env bash
# shellcheck disable=SC2001
# shellcheck disable=SC2016
set -uex
LOG_FILE_NAME="yasim_generate_gene_depth.log" \
    python -m yasim generate_gene_depth \
    -g ../ce11_as_2.gtf \
    -o ce11_as_2_gene_depth_20.tsv \
    -d 20 \
    --low_cutoff 1.0 \
    --high_cutoff_ratio 200
LOG_FILE_NAME="yasim_generate_isoform_depth.log" \
    python -m yasim generate_isoform_depth \
    -g ../ce11_as_2.gtf \
    -d ce11_as_2_gene_depth_20.tsv \
    -o ce11_as_2_isoform_depth_20.tsv \
    --alpha 4 \
    --low_cutoff 1.0 \
    --high_cutoff_ratio 200

function perform_housekeeping() {
    touch "${1}".finished || return 1
}

function perform_pbsim3_RSII_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_0.0_0.0_pbsim3_RSII_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        --strategy wgs \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 60 \
        --not_perform_assemble || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_RSII_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_0.0_0.0_pbsim3_RSII_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        --strategy wgs \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 60 \
        --not_perform_assemble || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_SEQUEL_CCS_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_0.0_0.0_pbsim3_SEQUEL_CCS
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        --strategy wgs \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 60 \
        --not_perform_assemble || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_pbsim3_SEQUEL_CLR_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_0.0_0.0_pbsim3_SEQUEL_CLR
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        --strategy wgs \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 60 \
        --not_perform_assemble || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME=ce11_as_2_rcompl_0.0_0.0_pbsim2_"${1}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${1}" \
        -F ../ce11_trans_2.fa.d \
        -d ce11_as_2_isoform_depth_20.tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 60 \
        --not_perform_assemble || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}

function perform_simulation() {
    perform_pbsim3_RSII_CLR_simulation &
    perform_pbsim3_RSII_CCS_simulation &
    perform_pbsim3_SEQUEL_CCS_simulation &
    perform_pbsim3_SEQUEL_CLR_simulation &
    for pbsim2_mode in R94 R103; do
        perform_pbsim2_simulation "${pbsim2_mode}" &
    done
    wait
}

perform_simulation

perform_assemble() {
    for fn in ce11_as_2_rcompl_0.0_0.0_*.d; do
        python -m yasim assemble \
            -F ../ce11_trans_2.fa.d \
            -d ce11_as_2_isoform_depth_20.tsv \
            -i "${fn}" \
            -o "$(echo "${fn/0\.0_0\.0/${1}_${2}}" | sed 's;\.d;;')" \
            --simulator_name "$(echo "${fn}" | sed -E 's;ce11_as_2_rcompl_0\.0_0\.0_(.*)\.d;\1;')" \
            --truncate_ratio_3p "${1}" \
            --truncate_ratio_5p "${2}" &
    done
}

perform_assemble 0.0 0.4
perform_assemble 0.2 0.2
perform_assemble 0.4 0.0

perform_assemble 0.0 0.2
perform_assemble 0.1 0.1
perform_assemble 0.2 0.0
perform_assemble 0.0 0.0

wait

for fn in *.fq; do
    {
        python -m labw_utils.bioutils describe_fastq "${fn}" &>/dev/null &&
            echo "FIN ${fn}" ||
            echo "ERR ${fn}"
    } &
done
wait
rm ./*.parquet
Rscript plot_fastq.R

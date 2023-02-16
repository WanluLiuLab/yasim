#!/usr/bin/env bash
function generate_depth(){
    python -m yasim generate_gene_depth \
        -g ../ce11_as_2.gtf.gz \
        -o ce11_as_2_gene_depth_"${1}".tsv.xz \
        -d "${1}"
    python -m yasim generate_isoform_depth \
        -g ../ce11_as_2.gtf.gz \
        -d ce11_as_2_gene_depth_"${1}".tsv.xz \
        -o ce11_as_2_isoform_depth_"${1}".tsv.xz
}

function perform_simulation(){
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_RSII_CLR \
        -j 40 &
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_RSII_CCS \
        -j 40 &
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_SEQUEL_CLR \
        -j 40 &
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_SEQUEL_CCS \
        -j 40 &
    for pbsim2_mode in R94 R103; do
        python -m yasim pbsim2 \
            -e ../bin/pbsim2 \
            -m "${pbsim2_mode}" \
            -F ../ce11_trans_2.fa.d \
            -d "${1}".tsv.xz \
            -o "${2}"_pbsim2_"${pbsim2_mode}" \
            -j 40 &
    done
    wait
}
cd diff_read_depth || exit 1

for mean_depth in 10 25 40 55 70 85 100; do
    # generate_depth "${mean_depth}" &
    perform_simulation ce11_as_2_isoform_depth_"${mean_depth}" ce11_as_2_isoform_depth_"${mean_depth}"
done

wait

cd .. || exit 1

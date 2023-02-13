#!/usr/bin/env bash
set -uex
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz

# Generate AS events with different complexity level (defaults 2)
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq.gtf.gz
xz -9 -T0 ce11.ncbiRefSeq.gtf.gz.*.tsv -vv

function generate_as_events(){
    python -m yasim generate_as_events \
        -f ce11.fa.gz -g ce11.ncbiRefSeq.gtf.gz \
        -o ce11_as_"${1}".gtf.gz \
        --complexity "${1}"
    python -m labw_utils.bioutils describe_gtf ce11_as_"${1}".gtf.gz
    python -m yasim transcribe \
        -g ce11_as_"${1}".gtf.gz \
        -f ce11.fa.gz \
        -o ce11_trans_"${1}".fa
    xz -9 -T0 ce11_as_"${1}".gtf.gz.*.tsv -vv
}

function perform_simulation(){
        python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_RSII_CLR \
        -j 40
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m RSII \
        -M qshmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_RSII_CCS \
        -j 40
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_SEQUEL_CLR \
        -j 40
    python -m yasim pbsim3 \
        -e ../bin/pbsim3 \
        -m SEQUEL \
        -M errhmm \
        -F ../ce11_trans_2.fa.d \
        --ccs_pass 10 \
        -d "${1}".tsv.xz \
        -o "${2}"_pbsim3_SEQUEL_CCS \
        -j 40
    for pbsim2_mode in R94 R103; do
        python -m yasim pbsim2 \
            -e ../bin/pbsim2 \
            -m "${pbsim2_mode}" \
            -F ../ce11_trans_2.fa.d \
            -d "${1}".tsv.xz \
            -o "${2}"_pbsim2_"${pbsim2_mode}" \
            -j 40
    done
}

for gene_complexity_level in 1 2 3 5 7 9; do
    generate_as_events "${gene_complexity_level}"
done

mkdir -p src bin
cd src || exit 1
git clone http://github.com/yukiteruono/pbsim3
cd pbsim3 || exit 1
git checkout v3.0.0
./configure
make -j20
ln -s "$(pwd)"/src/pbsim ../../bin/pbsim3
cd .. || exit 1
git clone http://github.com/yukiteruono/pbsim2
cd pbsim2 || exit 1
git checkout eeb5a19420534a0f672c81db2670117e62a9ee38
./configure
make -j20
ln -s "$(pwd)"/src/pbsim ../../bin/pbsim2
cd .. || exit 1
cd .. || exit 1

# Read depth
mkdir -p diff_read_depth
cd diff_read_depth || exit 1

for mean_depth in 10 25 40 55 70 85 100; do
    python -m yasim generate_gene_depth \
        -g ../ce11_as_2.gtf.gz \
        -o ce11_as_2_gene_depth_"${mean_depth}".tsv.xz \
        -d "${mean_depth}"
    python -m yasim generate_isoform_depth \
        -g ../ce11_as_2.gtf.gz \
        -d ce11_as_2_gene_depth_"${mean_depth}".tsv.xz \
        -o ce11_as_2_isoform_depth_"${mean_depth}".tsv.xz
    perform_simulation ce11_as_2_isoform_depth_"${mean_depth}" ce11_as_2_isoform_depth_"${mean_depth}" &
done
wait


cd .. || exit 1

for accuracy in 0.6 0.7 0.8 0.9 1.0; do
    pass
done

for five_prime_break in ; do
    pass
done

# NIpG
# RLEN
# ERROR RATE
# RC


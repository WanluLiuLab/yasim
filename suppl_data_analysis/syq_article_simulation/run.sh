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
    xz -9 -vv ce11_as_"${1}".gtf.gz.*.tsv
}

for gene_complexity_level in 1 2 3 5 7 9; do
    generate_as_events "${gene_complexity_level}" &
done
wait

# Install LLRGs

mkdir -p src bin
cd src || exit 1
git clone http://github.com/yukiteruono/pbsim3
cd pbsim3 || exit 1
git checkout f3e5f1cf3d0e8346b5e4598ac238b2b570b223e8
./configure
make -j20
ln -sf "$(pwd)"/src/pbsim ../../bin/pbsim3
cd .. || exit 1
git clone http://github.com/yukiteruono/pbsim2
cd pbsim2 || exit 1
git checkout eeb5a19420534a0f672c81db2670117e62a9ee38
./configure
make -j20
ln -sf "$(pwd)"/src/pbsim ../../bin/pbsim2
cd .. || exit 1
cd .. || exit 1

# Read depth
cd diff_read_depth || exit 1
bash run.sh
cd .. || exit 1

# Read accuracy
cd diff_accuracy || exit 1
bash run.sh
cd .. || exit 1

# Genome Complexity
cd diff_genome_complexity || exit 1
bash run.sh
cd .. || exit 1



for five_prime_break in ; do
    pass
done

# NIpG
# RLEN
# RC


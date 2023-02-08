#!/usr/bin/env bash
set -ue

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
gunzip ./*.gz

grep '^chr1\s' hg38.ncbiRefSeq.gtf > hg38.ncbiRefSeq.chr1.gtf
head hg38.fa -n "$(($(cat -n hg38.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  hg38.chr1.fa
python -m labw_utils.bioutils transcribe -g hg38.ncbiRefSeq.chr1.gtf -f hg38.chr1.fa -o hg38_trans.fa

shuf < hg38_trans.fa | python -m yasim_scripts filter_pbsim_reference | head -n 4000 > hg38_trans.filtered.fa

function perform_pbsim2_length() {
    mkdir -p pbsim2_length_"${1}"
    pbsim2 \
        --prefix pbsim2_length_"${1}"/ce11 \
        --depth 0.5 \
        --length-mean "${1}" \
        --hmm_model ../../src/yasim/llrg_adapter/pbsim2_dist/R103.model \
        hg38_trans.filtered.fa &> pbsim2_length_"${1}".log
    cat pbsim2_length_"${1}"/ce11_*.fastq > pbsim2_length_"${1}".fastq
    rm -rf pbsim2_length_"${1}"
    python -m labw_utils.bioutils describe_fastq pbsim2_length_"${1}".fastq
}

for i in 400 1200 2000 2800; do
    perform_pbsim2_length "${i}" &
done
wait

function perform_pbsim2_accuracy() {
    mkdir -p pbsim2_accuracy_"${1}"
    pbsim2 \
        --prefix pbsim2_accuracy_"${1}"/ce11 \
        --depth 0.5 \
        --accuracy-mean "${1}" \
        --hmm_model ../../src/yasim/llrg_adapter/pbsim2_dist/R103.model \
        hg38_trans.filtered.fa &> pbsim2_accuracy_"${1}".log
    cat pbsim2_accuracy_"${1}"/ce11_*.fastq > pbsim2_accuracy_"${1}".fastq
    rm -rf pbsim2_accuracy_"${1}"
    python -m labw_utils.bioutils describe_fastq pbsim2_accuracy_"${1}".fastq
}


for i in 0.60 0.70 0.80 0.90 1.0; do
    perform_pbsim2_accuracy "${i}" &
done
wait

mkdir -p lastdb
lastdb -v -P40 lastdb/hg38_trans.filtered hg38_trans.filtered.fa

for fn in pbsim2_accuracy_*.fastq; do
    echo "${fn}"
    last-train -Qfastx -P40 lastdb/hg38_trans.filtered "${fn}" > "${fn}"_trans.train
    lastal -P40 -Qfastx -m1 -l1  \
    -p "${fn}"_trans.train \
    lastdb/hg38_trans.filtered "${fn}" \
    > "${fn}"_trans.maf
done

function perform_pbsim3_multi_pass_parallel(){
    mkdir -p pbsim3_multi_pass_"${1}"
    pbsim3 \
    --method errhmm \
    --errhmm ../../src/yasim/llrg_adapter/pbsim3_dist/ERRHMM-SEQUEL.model \
    --strategy wgs \
    --genome hg38_trans.filtered.fa \
    --depth 0.5 \
    --prefix pbsim3_multi_pass_"${1}"/ce11  \
    --pass-num "${1}" &> pbsim3_multi_pass_"${1}".log
}

function perform_pbsim3_multi_pass() {
    samtools merge -f \
        -@ 40 \
        pbsim3_multi_pass_"${1}"/subreads.bam \
        pbsim3_multi_pass_"${1}"/ce11_*.sam
    ccs \
        --report-json pbsim3_multi_pass_"${1}"/ccs.report.json \
        --report-file pbsim3_multi_pass_"${1}"/ccs.report.txt \
        --log-level INFO \
        --log-file pbsim3_multi_pass_"${1}"/ccs.log \
        pbsim3_multi_pass_"${1}"/subreads.bam pbsim3_multi_pass_"${1}"/ccs.bam
    samtools fastq pbsim3_multi_pass_"${1}"/ccs.bam > pbsim3_multi_pass_"${1}".fastq
    rm -rf pbsim3_multi_pass_"${1}"
    python -m labw_utils.bioutils describe_fastq pbsim3_multi_pass_"${1}".fastq
}


for i in 1 10 20 30 40 50; do
    perform_pbsim3_multi_pass_parallel "${i}" &
done
wait

cat pbsim3_multi_pass_1"${1}"/ce11_*.fastq > pbsim3_multi_pass_1"${1}".fastq
rm -rf pbsim3_multi_pass_1"${1}"
python -m labw_utils.bioutils describe_fastq pbsim3_multi_pass_1"${1}".fastq

for i in 10 20 30 40 50; do
    perform_pbsim3_multi_pass "${i}"
done
wait

for fn in pbsim3_multi_pass*.fastq; do
    echo "${fn}"
    last-train -Qfastx -P40 lastdb/hg38_trans.filtered "${fn}" > "${fn}"_trans.train
    lastal -P40 -Qfastx -m1 -l1  \
    -p "${fn}"_trans.train \
    lastdb/hg38_trans.filtered "${fn}" \
    > "${fn}"_trans.maf
done

printf "FILENAME\tINSERTION\tDELETION\tMATCH\tSUBSTITUTION\n" > all_last_mapq.tsv
for fn in *.maf; do
    python -m extract_quality_from_maf "${fn}" >> all_last_mapq.tsv
done

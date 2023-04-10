#!/usr/bin/env bash
set -uex
export NUM_THREADS=40

axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz

grep -i '^chrI\s' <ce11.ncbiRefSeq.gtf >ce11.ncbiRefSeq.chr1.gtf
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1) - 1))" >ce11.chr1.fa
rm -f ce11.ncbiRefSeq.gtf ce11.fa

python -m labw_utils.bioutils transcribe \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -f ce11.chr1.fa \
    -o ce11_trans.chr1.fa

mkdir -p sim
rm -rf sim/* pbsim.log
for fn in ce11_trans.chr1.fa.d/*.fa; do
    pbsim2 \
        --hmm_model ../../src/yasim/llrg_adapter/pbsim2_dist/R103.model \
        --prefix sim/ce11_"$(basename "${fn}")" \
        --depth 5 \
        "${fn}" 2>&1 | tee -a pbsim.log || true
done

cat sim/*.fastq | pigz -9 >simulated.fastq.gz
cat sim/*.maf | pigz -9 >simulated_gt.maf.gz
rm sim -rf

mkdir -p lastdb
lastdb -v -P"${NUM_THREADS}" lastdb/ce11 ce11.chr1.fa
lastdb -v -P"${NUM_THREADS}" lastdb/ce11_t ce11_trans.chr1.fa
last-train -v \
    -Qfastx \
    -P"${NUM_THREADS}" \
    lastdb/ce11 simulated.fastq.gz >simulated.train
last-train -v \
    -Qfastx \
    -P"${NUM_THREADS}" \
    lastdb/ce11_t simulated.fastq.gz >simulated_t.train
lastal -v \
    -Qfastx \
    -P"${NUM_THREADS}" \
    -p simulated.train \
    -m100 \
    -j7 \
    --splice \
    lastdb/ce11 simulated.fastq.gz |
    last-map-probs /dev/stdin |
    pigz -9 >simulated.maf.gz

lastal -v \
    -Qfastx \
    -P"${NUM_THREADS}" \
    -p simulated_t.train \
    -m100 \
    -j7 \
    lastdb/ce11_t simulated.fastq.gz |
    last-map-probs /dev/stdin |
    pigz -9 >simulated_t.maf.gz

for fn in *.maf.gz; do
    echo "${fn}"
    python -m yasim_scripts extract_quality_from_maf "${fn}"
done
python analyse_pbsim_log.py

# simulated_gt.maf.gz {'I': '7.65%', 'D': '5.64%',  'M': '85.87%', 'S': '0.85%'}
# simulated.maf.gz    {'I': '5.35%', 'D': '3.38%',  'M': '89.55%', 'S': '1.72%'}
# simulated_t.maf.gz  {'I': '5.43%', 'D': '3.46%',  'M': '89.33%', 'S': '1.78%'}
# pbsim.log           {'I': '8.23%', 'D': '16.54%', 'M': '84.67%', 'S': '0.91%'}

#!/usr/bin/env bash
set -ue -o pipefail
# mkdir -p lastdb
# lastdb -P240 lastdb/ce11_trans_2 ../ce11_trans_2.fa

for fn in *CCS.fq; do
    echo "${fn}"
    last-train -Qfastx -P240 lastdb/ce11_trans_2 "${fn}" > "${fn}"_trans.train
    lastal -P240 -Qfastx -m100 -j7 \
    -p "${fn}"_trans.train \
    lastdb/ce11_trans_2 "${fn}" |\
    last-map-probs /dev/stdin > "${fn}"_trans.maf
done

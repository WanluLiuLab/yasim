#!/usr/bin/env bash
# shellcheck disable=SC2034
# shellcheck disable=SC2155

export NUM_THREADS=$(cat /proc/cpuinfo | grep '^processor' | wc -l)
export MAX_MEM=4
export PYTHONPATH="${SHDIR}/../src:${SHDIR}/../test:${PYTHONPATH:-}"
export MKL_DEBUG_CPU_TYPE=5
export OMP_NUM_THREADS="${NUM_THREADS}"
export MKL_NUM_THREADS="${NUM_THREADS}"
export OPENBLAS_NUM_THREADS="${NUM_THREADS}"
export BLIS_NUM_THREADS="${NUM_THREADS}"
export USE_PYTHON="$(which python3 || which python)"
export PYSPARK_DRIVER_PYTHON="${USE_PYTHON}"

for fn in real_stats/*.fastq.sam.stats.d/pileup_stat.tsv.gz; do
    spark-submit \
        --master local["${NUM_THREADS}"] \
        --driver-memory "${MAX_MEM}"G \
        --executor-cores "${NUM_THREADS}" \
        --total-executor-cores "${NUM_THREADS}" \
        --executor-memory "${MAX_MEM}"G \
        --deploy-mode client \
        sample_sam_statistics.py "${fn}" &
done
wait

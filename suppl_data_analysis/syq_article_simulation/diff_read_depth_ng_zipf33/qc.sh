#!/usr/bin/env bash
for fn in *.fq.gz; do
   { python -m labw_utils.bioutils describe_fastq "${fn}" &> /dev/null || echo "${fn}"; } &
done
wait

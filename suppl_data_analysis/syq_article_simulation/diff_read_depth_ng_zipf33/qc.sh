#!/usr/bin/env bash
for fn in *.fq; do
   {
       python -m labw_utils.bioutils describe_fastq "${fn}" &> /dev/null \
        && echo "FIN ${fn}" \
        || echo "ERR ${fn}"
    } &
done
wait

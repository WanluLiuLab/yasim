#!/usr/bin/env bash
set -uex

for fn in *.maf.gz; do
    python -m yasim_scripts extract_quality_from_maf "${fn}" > "${fn}".mapq.tsv  &
done
wait

printf "FILENAME\tINSERTION\tDELETION\tMATCH\tSUBSTITUTION\n" > all_last_mapq.tsv
cat ./*.maf.gz.mapq.tsv >> all_last_mapq.tsv
rm -fr ./*.maf.gz.mapq.tsv

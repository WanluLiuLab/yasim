#!/usr/bin/env bash
set -ue
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz

tar xzvf chromFa.tar.gz
gunzip ce11.ncbiRefSeq.gtf.gz
rm -f chromFa.tar.gz

for chrname in I II III IV V X M; do
    grep '^chr'"${chrname}"'\s' ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr"${chrname}".gtf
done


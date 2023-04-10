#!/usr/bin/env bash
set -ue
{
if [ ! -f ce11.ncbiRefSeq.chr1.gtf ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
    gunzip ce11.ncbiRefSeq.gtf.gz
    grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf
fi
if [ ! -f ce11.chr1.fa ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
    gunzip ce11.fa.gz
    head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa
fi
} &>> preparation.log
{
if [ ! -f ce11.ncbiRefSeq.chr1.as.gtf ]; then
python -m yasim generate_as_events \
        -f ce11.fa \
        -g ce11.ncbiRefSeq.chr1.gtf \
        -o ce11.ncbiRefSeq.chr1.as.gtf \
        -c 5
fi
} 2>&1 | grep -v 'inferred from feature transcript' || true >> generate_as_events.log


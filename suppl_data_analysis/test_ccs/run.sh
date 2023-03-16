#!/usr/bin/env bash
set -uex
export NUM_THREADS=40

axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz


grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" > ce11.chr1.fa
rm -f ce11.ncbiRefSeq.gtf ce11.fa
samtools faidx ce11.chr1.fa

python -m labw_utils.bioutils transcribe \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -f ce11.chr1.fa \
    -o ce11_trans.chr1.fa
python -m yasim generate_as_events \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq_as.chr1.gtf \
    -c 2
python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_trans_as.chr1.fa
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq_as.chr1.gtf
python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_gene_depth.chr1.tsv \
    -d 20
python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_isoform_depth.chr1.tsv \
    -d ce11_gene_depth.chr1.tsv \
    --alpha 4
python -m yasim pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -F ce11_trans_as.chr1.fa.d \
    -d ce11_isoform_depth.chr1.tsv \
    -o ce11_ccs \
    -j 40 \
    --ccs_pass 20

python -m yasim_scripts generate_missing_subreads_xml ce11_ccs.d/
/opt/yuzj/smrttools-11.1.0.166339/smrtcmds/bin/python3 \
    merge.py out.subreads.xml ce11_ccs.d/*/tmp.subreads.bam.xml

# FIXME: Error in following step
#~/opt/smrttools-11.1.0.166339/smrtcmds/bin/dataset consolidate \
#    out.subreads.xml \
#    out_consolidated.subreads.bam \
#    out_consolidated.subreads.xml
# Requires CCS 6.0.0. 6.5.0 would FAIL.
ccs \
    --report-json ccs.json \
    --report-file ccs.txt \
    --log-level INFO \
    --log-file ccs.log \
    --num-threads 40 \
    out_consolidated.subreads.xml \
    out.ccs.xml

isoseq3 cluster \
    out.ccs.xml \
    out.transcripts.xml \
    --log-level INFO \
    --log-file isoseq3_cluster.log \
    --num-threads 40

pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --log-file aln.log \
    --log-level INFO \
    out.transcripts.xml.hq.bam \
    ce11.chr1.fa \
    aln.bam

isoseq3 collapse \
    --do-not-collapse-extra-5exons \
    --log-level INFO \
    --log-file isoseq3_collapse.log \
    aln.bam out.ccs.xml collapse.gff
gffcompare collapse.gff -r ce11.ncbiRefSeq_as.chr1.gtf -o gffcmp
pigeon sort ce11.ncbiRefSeq.chr1.gtf -o ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf
pigeon index ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf

pigeon classify collapse.gff ce11.ncbiRefSeq.chr1.pigeon_sorted.gtf ce11.chr1.fa -o classified
pigeon filter classified_classification.txt --isoforms collapse.gff

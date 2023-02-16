#!/usr/bin/env bash
set -ue
python -m yasim transcribe \
    -g ../ce11.ncbiRefSeq.chrM.gtf \
    -f ../chrM.fa \
    -o chrM.trans.fa

python -m yasim generate_gene_depth \
    -g ../ce11.ncbiRefSeq.chrM.gtf \
    -o depth_gene.tsv \
    --mu 10

python -m yasim generate_isoform_depth \
    -g ../ce11.ncbiRefSeq.chrM.gtf \
    -o depth.tsv \
    -d depth_gene.tsv

python -m yasim pbsim3 \
    -e /mnt/volume2/TGS/phase_2_simulation/yasim-dev/suppl_data_analysis/syq_article_simulation/bin/pbsim3 \
    -m RSII \
    -M qshmm \
    -F chrM.trans.fa.d \
    -d depth.tsv \
    -o pbsim3_RSII_CLR \
    -j 40 \
    --accuracy-mean 0.4
python -m yasim pbsim3 \
    -e /mnt/volume2/TGS/phase_2_simulation/yasim-dev/suppl_data_analysis/syq_article_simulation/bin/pbsim3 \
    -m RSII \
    -M qshmm \
    -F chrM.trans.fa.d \
    --ccs_pass 10 \
    -d depth.tsv \
    -o pbsim3_RSII_CCS \
    -j 40
python -m yasim pbsim3 \
    -e /mnt/volume2/TGS/phase_2_simulation/yasim-dev/suppl_data_analysis/syq_article_simulation/bin/pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -F chrM.trans.fa.d \
    -d depth.tsv \
    -o pbsim3_SEQUEL_CLR \
    -j 40
python -m yasim pbsim3 \
    -e /mnt/volume2/TGS/phase_2_simulation/yasim-dev/suppl_data_analysis/syq_article_simulation/bin/pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -F chrM.trans.fa.d \
    --ccs_pass 10 \
    -d depth.tsv \
    -o pbsim3_SEQUEL_CCS \
    -j 40
for pbsim2_mode in R94 R103; do
    python -m yasim pbsim2 \
        -e /mnt/volume2/TGS/phase_2_simulation/yasim-dev/suppl_data_analysis/syq_article_simulation/bin/pbsim2 \
        -m "${pbsim2_mode}" \
        -F chrM.trans.fa.d \
        -d depth.tsv \
        -o pbsim2_"${pbsim2_mode}" \
        -j 40
done
for fn in *.fq; do
    python -m labw_utils.bioutils describe_fastq "${fn}" &
done
wait

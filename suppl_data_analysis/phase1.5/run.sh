#!/usr/bin/env bash
set -uex
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz

cat ce11.ncbiRefSeq.gtf | grep -i '^chrI\s' > ce11.ncbiRefSeq.chr1.gtf
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa

rm -f ce11.ncbiRefSeq.gtf ce11.fa

if ! mamba env list | grep ^yasim-pbsim1 &>> /dev/null; then
    mamba create -y -n yasim-pbsim1 -c bioconda pbsim=1.0.3
fi
if ! mamba env list | grep ^yasim-pbsim2 &>> /dev/null; then
    mamba create -y -n yasim-pbsim2 -c bioconda pbsim2=2.0.1
fi
if ! mamba env list | grep ^yasim-badread &>> /dev/null; then
    mamba create -y -n yasim-badread -c bioconda badread=0.2.0 python-edlib
fi

PYTHONPATH=../../src python -m yasim generate_as_events \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq.chr1_as.gtf
PYTHONPATH=../../src python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq.chr1_as.gtf \
    -o ce11_gene_depth.chr1.tsv \
    -d 10.0
PYTHONPATH=../../src python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq.chr1_as.gtf \
    -o ce11_depth_diu1.chr1.tsv \
    -d ce11_gene_depth.chr1.tsv
PYTHONPATH=../../src python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq.chr1_as.gtf \
    -o ce11_depth_diu2.chr1.tsv \
    -d ce11_gene_depth.chr1.tsv
PYTHONPATH=../../src python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1_as.gtf \
    -o ce11_trans.chr1_as.fa
python -m labw_utils.bioutils get_gtf_statistics \
    -g ce11.ncbiRefSeq.chr1_as.gtf \
    -o ce11.ncbiRefSeq.chr1_as.gtf

# PBSIM 1, 2 and 3 is lightweighted so can be run parallely

for diu_grp in diu1 diu2; do
    PYTHONPATH=../../src python -m yasim pbsim \
        -d ce11_depth_"${diu_grp}".chr1.tsv \
        -F ce11_trans.chr1_as.fa.d \
        -e pbsim \
        -o ce11_pbsim_clr_"${diu_grp}" &
    PYTHONPATH=../../src python -m yasim pbsim \
        -d ce11_depth_"${diu_grp}".chr1.tsv \
        -F ce11_trans.chr1_as.fa.d \
        -e pbsim \
        --ccs \
        -o ce11_pbsim_ccs_"${diu_grp}" &

    for pbsim2_models in R95 P6C4 R94 P5C3 R103 P4C2; do
        PYTHONPATH=../../src python -m yasim pbsim2 \
        -d ce11_depth_"${diu_grp}".chr1.tsv \
        -F ce11_trans.chr1_as.fa.d \
        -e pbsim2 \
        -o ce11_pbsim2_"${pbsim2_models}"_"${diu_grp}" -m "${pbsim2_models}" &
    done

    for pbsim3_models in RSII ONT; do
        PYTHONPATH=../../src python \
        -m yasim pbsim3 \
        -d ce11_depth_"${diu_grp}".chr1.tsv \
        -F ce11_trans.chr1_as.fa.d \
        -e pbsim3 \
        -o ce11_pbsim3_"${pbsim3_models}"_"${diu_grp}" \
        -m "${pbsim3_models}" &
    done
done
wait

python -m labw_utils.bioutils describe_fastq ./*.fq

PYTHONPATH=../../src python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11_trans.chr1.fa
python -m labw_utils.bioutils get_gtf_statistics \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq.chr1.gtf

for fn in *.fq; do
    minimap2 -x splice -a -t 50 ce11.chr1.fa "${fn}" > "${fn}".sam
    samtools sort "${fn}".sam -@ 50 -o "${fn}".bam
    samtools index "${fn}".bam
    rm -f "${fn}".sam
    minimap2 -a -t 50 ce11_trans.chr1.fa "${fn}" > "${fn/.fq/_trans.fq}".sam
    samtools sort "${fn/.fq/_trans.fq}".sam -@ 50 -o "${fn/.fq/_trans.fq}".bam
    samtools index "${fn/.fq/_trans.fq}".bam
    rm -f "${fn/.fq/_trans.fq}".sam
done

for fn in *.fq.bam; do
    python -m labw_utils.bioutils describe_sam "${fn}"
    Rscript R/transform_depth_results.R "${fn}.stats.d"/pileup_stat.tsv.gz "${fn}.stats.d"/pileup_stat.merged.tsv
done

find . | grep .fq.bam$ | grep -v trans | while read -r fn; do
    featureCounts -O -M --primary --ignoreDup -a ce11.ncbiRefSeq.chr1.gtf -g transcript_id -L -o "${fn}".fc.tsv "${fn}"
    stringtie -L -G ce11.ncbiRefSeq.chr1.gtf -o "${fn}".stringtie.gtf -p 50 "${fn}"
    python -m labw_utils.bioutils get_gtf_statistics -g "${fn}".stringtie.gtf -o "${fn}".stringtie.gtf
done

Rscript ./plot_fastq.R
Rscript ./plot_nipg.R
Rscript ./plot_sam.R
Rscript ./plot_gep.R


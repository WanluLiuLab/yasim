#!/usr/bin/env bash
set -uex
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz

if ! mamba env list | grep ^yasim-pbsim1 &>>/dev/null; then
    mamba create -y -n yasim-pbsim1 -c bioconda pbsim=1.0.3
fi
if ! mamba env list | grep ^yasim-pbsim2 &>>/dev/null; then
    mamba create -y -n yasim-pbsim2 -c bioconda pbsim2=2.0.1
fi
if ! mamba env list | grep ^yasim-badread &>>/dev/null; then
    mamba create -y -n yasim-badread -c bioconda badread=0.2.0 python-edlib
fi

PYTHONPATH=../../src python -m yasim generate_as_events -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq_as.gtf
PYTHONPATH=../../src python -m yasim generate_depth -g ce11.ncbiRefSeq_as.gtf -o ce11_depth.tsv -d 20.0
PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq_as.gtf -o ce11_trans_as.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq_as.gtf -o ce11.ncbiRefSeq_as.gtf

# PBSIM 1, 2 and 3 is lightweighted so can be run parallely
PYTHONPATH=../../src python -m yasim pbsim -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim.sh" -o ce11_pbsim_clr &
PYTHONPATH=../../src python -m yasim pbsim -d ce11_depth.tsv -F ce11_trans_as.fa.d --ccs -e "wrappers/pbsim.sh" -o ce11_pbsim_ccs &
for pbsim2_models in R95 P6C4 R94 P5C3 R103 P4C2; do
    PYTHONPATH=../../src python -m yasim pbsim2 -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim2.sh" -j 40 -o ce11_pbsim2_"${pbsim2_models}" -m "${pbsim2_models}" &
done
for pbsim3_models in RSII ONT; do
    PYTHONPATH=../../src python -m yasim pbsim3 -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "pbsim3" -j 40 -o ce11_pbsim3_"${pbsim3_models}" -m "${pbsim3_models}" &
done
wait

# Badread consumes lot more resources and should be run sequentially
for badread_models in nanopore2018 pacbio2016 nanopore2020 verybad verynice; do
    PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_"${badread_models}" -m "${badread_models}"
done

python -m labw_utils.bioutils describe_fastq ./*.fq

for fn in *.fq.stats.d; do
    Rscript R/plot_describe_fastq.R "${fn}"
done

PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11_trans.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq.gtf

for fn in *.fq; do
    minimap2 -x splice -a -t 50 ce11.fa "${fn}" >"${fn}".sam
    samtools sort "${fn}".sam -@ 50 -o "${fn}".bam
    samtools index "${fn}".bam
    rm -f "${fn}".sam
    minimap2 -a -t 50 ce11_trans.fa "${fn}" >"${fn/.fq/_trans.fq}".sam
    samtools sort "${fn/.fq/_trans.fq}".sam -@ 50 -o "${fn/.fq/_trans.fq}".bam
    samtools index "${fn/.fq/_trans.fq}".bam
    rm -f "${fn/.fq/_trans.fq}".sam
done

for fn in *.fq.bam; do
    python -m labw_utils.bioutils describe_sam "${fn}"
    Rscript R/plot_describe_sam.R "${fn}.stats.d"
done

find . | grep .fq.bam$ | grep -v trans | while read -r fn; do
    featureCounts -O -M --primary --ignoreDup -a ce11.ncbiRefSeq.gtf -g transcript_id -L -o "${fn}".fc.tsv "${fn}"
    stringtie -L -G ce11.ncbiRefSeq.gtf -o "${fn}".stringtie.gtf -p 50 "${fn}"
    python -m labw_utils.bioutils get_gtf_statistics -g "${fn}".stringtie.gtf -o "${fn}".stringtie.gtf
done

find . | grep .fq.bam$ | grep trans | while read -r fn; do
    printf "REFERENCE_NAME\tPOS\tNUM_READS\n" >"${fn}".depth.tsv
    samtools depth -aa "${fn}" >>"${fn}".depth.tsv
    Rscript R/transform_depth_results.R "${fn}".depth.tsv "${fn}".depth.mean.tsv
done

Rscript ./plot_fastq.R
Rscript ./plot_nipg.R
Rscript ./plot_sam.R
Rscript ./plot_gep.R

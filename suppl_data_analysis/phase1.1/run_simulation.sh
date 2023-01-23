#!/usr/bin/env bash
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz

if ! mamba env list | grep ^yasim-pbsim1 &>> /dev/null; then
    mamba create -y -n yasim-pbsim1 -c bioconda pbsim=1.0.3
fi
if ! mamba env list | grep ^yasim-pbsim2 &>> /dev/null; then
    mamba create -y -n yasim-pbsim2 -c bioconda pbsim2=2.0.1
fi
if ! mamba env list | grep ^yasim-badread &>> /dev/null; then
    mamba create -y -n yasim-badread -c bioconda badread=0.2.0 python-edlib
fi

PYTHONPATH=../../src python -m yasim generate_as_events -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq_as.gtf
PYTHONPATH=../../src python -m yasim generate_depth -g ce11.ncbiRefSeq_as.gtf -o ce11_depth.tsv -d 20
PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq_as.gtf -o ce11_trans_as.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq_as.gtf -o ce11.ncbiRefSeq_as.gtf
PYTHONPATH=../../src python -m yasim pbsim -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim.sh" -j 40 -o ce11_pbsim_clr
PYTHONPATH=../../src python -m yasim pbsim2 -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim2.sh" -j 40 -o ce11_pbsim2_r94 -m R94
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_nanopore2018 -m nanopore2018
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_pacbio2016 -m pacbio2016
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_nanopore2020 -m nanopore2020

PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11_trans.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq.gtf

STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ce11_STAR_index --genomeFastaFiles ce11.fa --sjdbGTFfile ce11.ncbiRefSeq.gtf
# !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=100286401, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 12
mkdir -p ce11_trans_bwa_index/
bwa index ce11_trans.fa -p ce11_trans_bwa_index/ce11



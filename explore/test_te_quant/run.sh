#!/usr/bin/env bash
set -ue
mkdir -p ref aln 
# Get Reference
cd ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
cat ce11.ncbiRefSeq.gtf.gz | pigz -d >ce11.ncbiRefSeq.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
cat ce11.fa.gz | pigz -d >ce11.fa
samtools faidx ce11.fa

mkdir -p rmsk_wd
cd rmsk_wd
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromOut.tar.gz
tar xzvf chromOut.tar.gz
cat */*.out > ce11.out
perl \
    ../../../../deps/RepeatMasker/util/rmOutToGFF3.pl \
    ce11.out > ce11.out.gff
cd ../../../
python rmsk_filter.py
cd -
bedtools getfasta \
    -fi ../ce11.fa \
    -fo ../ce11.te_loci.fa \
    -bed ce11.out.filtered.gtf \
    -name \
    -s
cd ..

wget https://dfam.org/releases/current/families/Dfam_curatedonly.h5.gz
python -m yasim generate_te_index \
    -i Dfam_curatedonly.h5.gz \
    -o teidx.pkl.xz \
    -t 6239
cd ..
sbatch <run_aln_build_index_slrum.sh

# Simulating data
mkdir -p sim
cd sim
python -m yasim generate_gene_te_fusion \
    -g ../ref/ce11.ncbiRefSeq.gtf \
    -f ../ref/ce11.fa \
    --tedb ../ref/teidx.pkl.xz \
    -o ce11_denovo_test \
    -d 150 \
    -n 100000
samtools faidx ce11_denovo_test.fa
python -m labw_utils.bioutils split_fasta ce11_denovo_test.fa
cd ..
python plot_n_frags.py

sbatch <run_art_slrum.sh
sbatch <run_pbsim_slrum.sh

seqkit fq2fa sim/ce11_denovo_test_pbsim.fq >sim/ce11_denovo_test_pbsim.fa

# Running Aligner
rm -fr aln
mkdir -p aln
sbatch <run_sim_aln_slrum.sh

# Run RepeatMasker
python prepare_for_rmsk.py sim/ce11_denovo_test_pbsim
python prepare_for_rmsk.py sim/ce11_denovo_test
sbatch <run_sim_rmsk_slrum.sh
python rmsk_cmp.py
python plot_rmsk.py

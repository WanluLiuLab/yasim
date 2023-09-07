#!/usr/bin/env bash
#SBATCH --job-name=run_rnaspades
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-spades

rnaspades.py \
    --threads 40 \
    --memory 100 \
    -o assmb/rnaspades_ngs_de_novo \
    -1 sim/ce11_denovo_test_art_1.fq \
    -2 sim/ce11_denovo_test_art_2.fq &>> \
    assmb/rnaspades_ngs_de_novo.log
# TGS only not supported.
# rnaspades.py \
#     --threads 40 \
#     --memory 100 \
#     -o assmb/rnaspades_tgs_de_novo \
#     --pacbio sim/ce11_denovo_test_pbsim.fa &>> \
#     assmb/rnaspades_tgs_de_novo.log
rnaspades.py \
    --threads 40 \
    --memory 100 \
    -o assmb/rnaspades_ngs_tgs_de_novo \
    -1 sim/ce11_denovo_test_art_1.fq \
    -2 sim/ce11_denovo_test_art_2.fq \
    --pacbio sim/ce11_denovo_test_pbsim.fa &>> \
    assmb/rnaspades_ngs_tgs_de_novo.log

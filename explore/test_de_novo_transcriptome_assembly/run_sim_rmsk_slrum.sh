#!/usr/bin/env bash
#SBATCH --job-name=run_sim_rmsk
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=36:00:00
set -e

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"

conda activate yasim_dev

perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir aln/ce11_denovo_test_rmsk_prep.rmsk.d \
    -gff -small -s \
    sim/ce11_denovo_test_rmsk_prep.fa

perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir aln/ce11_denovo_test_pbsim_rmsk_prep.rmsk.d \
    -gff -small -s \
    sim/ce11_denovo_test_pbsim_rmsk_prep.fa

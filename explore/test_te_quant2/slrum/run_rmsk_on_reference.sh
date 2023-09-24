#!/usr/bin/env bash
#SBATCH --job-name=rum_rmsk_ce11
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
    -dir ref/ce11.rmsk_rmblast.d \
    -gff -small -s \
    ref/ce11.fa

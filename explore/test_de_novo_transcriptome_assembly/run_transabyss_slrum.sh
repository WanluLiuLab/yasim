#!/usr/bin/env bash
#SBATCH --job-name=run_transabyss
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-transabyss

transabyss \
    --threads 40 \
    --outdir assmb/transabyss_ngs_de_novo \
    --pe sim/ce11_denovo_test_art_1.fq sim/ce11_denovo_test_art_2.fq &> \
    assmb/transabyss_ngs_de_novo.log
transabyss \
    --threads 40 \
    --outdir assmb/transabyss_tgs_de_novo \
    --se sim/ce11_denovo_test_pbsim.fa &> \
    assmb/transabyss_tgs_de_novo.log
transabyss \
    --threads 40 \
    --outdir assmb/transabyss_ngs_tgs_de_novo \
    --se sim/ce11_denovo_test_pbsim.fa \
    --pe sim/ce11_denovo_test_art_1.fq sim/ce11_denovo_test_art_2.fq &> \
    assmb/transabyss_ngs_tgs_de_novo.log

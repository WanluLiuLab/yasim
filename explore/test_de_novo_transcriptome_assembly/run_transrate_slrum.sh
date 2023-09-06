#!/usr/bin/env bash
#SBATCH --job-name=run_transrate
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-transrate

FILES="$(echo assmb_final/*.fa| sed 's; ;,;g')"
echo "${FILES}" | sed 's;,;\n;g'
rm -fr assmb_final/transrate

transrate \
    --assembly="${FILES}" \
    --reference=sim/ce11_denovo_test.fa \
    --output=assmb_final/transrate \
    --threads=40

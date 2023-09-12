#!/usr/bin/env bash
#SBATCH --job-name=run_pbsim
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim_dev
CURDIR="$(pwd)"
cd ../../ && . bashrc && cd "${CURDIR}" || exit 1

"${CONDA_PREFIX}/bin/python" -m yasim pbsim \
    -F ref/WBcel235.trans.fa.d \
    -j 40 \
    -d sim/WBcel235.isoform_depth.tsv \
    -o sim/WBcel235_pbsim \
    --ccs
rm -fr sim/WBcel235_pbsim.d

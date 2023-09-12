#!/usr/bin/env bash
#SBATCH --job-name=run_art
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

"${CONDA_PREFIX}/bin/python" -m yasim art \
    -F ref/WBcel235.trans.fa.d \
    -j 40 \
    -d sim/WBcel235.isoform_depth.tsv \
    -o sim/WBcel235_art \
    --is_pair_end \
    --pair_end_fragment_length_std 50 \
    --pair_end_fragment_length_mean 200 \
    --sequencer_name HS25 \
    --read_length 150
rm -fr sim/WBcel235_art.d

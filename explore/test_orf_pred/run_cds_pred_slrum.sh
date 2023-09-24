#!/usr/bin/env bash
#SBATCH --job-name=run_cds_pred
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00
set -e

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim_dev

# cd assmb_transcripts_final
# for fn in *.fa; do
#     ../../../deps/GMST/gmst.pl \
#         --output ../assmb_final/"$(basename "${fn}")".gff \
#         --faa \
#         --format GFF \
#         "${fn}" &
# done
# wait
# cd -

plass nuclassemble \
    sim/WBcel235_pbsim.fq \
    assmb_final/plass_tgs_de_novo.faa \
    assmb/tmp \
    --threads 40

plass nuclassemble \
    sim/WBcel235_art_1.fq sim/WBcel235_art_2.fq \
    assmb_final/plass_ngs_de_novo.faa \
    assmb/tmp \
    --threads 40

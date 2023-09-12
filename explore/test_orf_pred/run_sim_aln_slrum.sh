#!/usr/bin/env bash
#SBATCH --job-name=run_aln_build_index
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00
set -e

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"

conda activate yasim_dev

STAR \
    --runThreadN 40 \
    --genomeDir ref/WBcel235_STAR_idx \
    --readFilesIn sim/WBcel235_art_1.fq sim/WBcel235_art_2.fq \
    --outSAMtype BAM Unsorted \
    --outSAMattributes All \
    --outFileNamePrefix aln/WBcel235_art_STAR/ \
    --twopassMode Basic
samtools sort \
    -@ 40 \
    -o aln/WBcel235_art.bam \
    aln/WBcel235_art_STAR/Aligned.out.bam
samtools index aln/WBcel235_art.bam

minimap2 \
    -x splice \
    -a \
    -t 40 \
    ref/WBcel235.genome.fa \
    sim/WBcel235_pbsim.fq >aln/WBcel235_pbsim.sam
samtools sort \
    -@ 40 \
    -o aln/WBcel235_pbsim.bam \
    aln/WBcel235_pbsim.sam
samtools index aln/WBcel235_pbsim.bam
rm -fr aln/WBcel235_pbsim.sam

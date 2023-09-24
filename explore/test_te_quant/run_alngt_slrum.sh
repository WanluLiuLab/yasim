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

minimap2 -t 40 \
    -a ref/ce11.te_loci.fa sim/ce11_denovo_test.fa \
    >aln/ce11_denovo_test.minimap2.te_loci.sam
samtools sort -@40 \
    aln/ce11_denovo_test.minimap2.te_loci.sam \
    -o aln/ce11_denovo_test.minimap2.te_loci.bam
samtools index aln/ce11_denovo_test.minimap2.te_loci.bam

minimap2 -t 40 \
    -a ref/ce11.fa sim/ce11_denovo_test.fa \
    >aln/ce11_denovo_test.minimap2.ref_genome.sam
samtools sort -@40 \
    aln/ce11_denovo_test.minimap2.ref_genome.sam \
    -o aln/ce11_denovo_test.minimap2.ref_genome.bam
samtools index aln/ce11_denovo_test.minimap2.ref_genome.bam
mkdir -p aln/ce11_denovo_test.minimap2.ref_genome.fc_assign.d
featureCounts -L -O -M --primary --ignoreDup \
    -a ref/rmsk_wd/ce11.out.gff \
    -g Target \
    -t dispersed_repeat \
    -o aln/ce11_denovo_test.minimap2.ref_genome.fc.tsv \
    aln/ce11_denovo_test.minimap2.ref_genome.bam \
    -R CORE \
    --Rpath aln/ce11_denovo_test.minimap2.ref_genome.fc_assign.d

pblat \
    -t=dna \
    -q=dna \
    -threads=40 \
    -out=blast9 \
    ref/ce11.te_loci.2bit \
    sim/ce11_denovo_test.fa \
    aln/ce11_denovo_test.blat.te_loci.blast9

minimap2 -t 40 \
    -a ref/teidx.pkl.xz.fa sim/ce11_denovo_test.fa \
    >aln/ce11_denovo_test.minimap2.te_consensus.sam
samtools sort -@40 \
    aln/ce11_denovo_test.minimap2.te_consensus.sam \
    -o aln/ce11_denovo_test.minimap2.te_consensus.bam
samtools index aln/ce11_denovo_test.minimap2.te_consensus.bam

pblat \
    -t=dna \
    -q=dna \
    -threads=40 \
    -out=blast9 \
    ref/teidx.2bit \
    sim/ce11_denovo_test.fa \
    aln/ce11_denovo_test.blat.te_consensus.blast9

perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine rmblast \
    -parallel 40 \
    -dir aln/ce11_denovo_test_rmsk_prep.rmsk_rmblast.d \
    -gff -small -s \
    sim/ce11_denovo_test_rmsk_prep.fa
perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine hmmer \
    -parallel 40 \
    -dir aln/ce11_denovo_test_rmsk_prep.rmsk_hmmer.d \
    -gff -small -s \
    sim/ce11_denovo_test_rmsk_prep.fa

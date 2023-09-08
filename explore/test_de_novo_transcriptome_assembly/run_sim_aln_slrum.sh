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

# STAR \
#     --runThreadN 40 \
#     --genomeDir ref/ce11_STAR_idx \
#     --readFilesIn sim/ce11_denovo_test_art_1.fq sim/ce11_denovo_test_art_2.fq \
#     --outSAMtype BAM Unsorted \
#     --outSAMattributes All \
#     --outFileNamePrefix aln/ce11_denovo_test_art_STAR/ \
#     --twopassMode Basic
# samtools sort \
#     -@ 40 \
#     -o aln/ce11_denovo_test_art.bam \
#     aln/ce11_denovo_test_art_STAR/Aligned.out.bam
# samtools index aln/ce11_denovo_test_art.bam

minimap2 \
    -x splice \
    -a \
    -t 40 \
    ref/ce11.fa \
    sim/ce11_denovo_test_pbsim.fq >aln/ce11_denovo_test_pbsim.sam
samtools sort \
    -@ 40 \
    -o aln/ce11_denovo_test_pbsim.bam \
    aln/ce11_denovo_test_pbsim.sam
samtools index aln/ce11_denovo_test_pbsim.bam
rm -fr aln/ce11_denovo_test_pbsim.sam

# Bias towards NGS; discarded
# conda activate yasim-detonate
# mkdir -p sim/ce11_denovo_test_rsem
# rsem-prepare-reference \
#     -p 40 \
#     --bowtie \
#     sim/ce11_denovo_test.fa \
#     sim/ce11_denovo_test_rsem/ref
# rsem-calculate-expression \
#     -p 40 \
#     --paired-end \
#     sim/ce11_denovo_test_art_1.fq \
#     sim/ce11_denovo_test_art_2.fq \
#     sim/ce11_denovo_test_rsem/ref \
#     sim/ce11_denovo_test_rsem/expr_ngs
# samtools sort \
#     -@ 40 \
#     sim/ce11_denovo_test_rsem/expr_ngs.transcript.bam \
#     -o sim/ce11_denovo_test_rsem/expr_ngs.transcript.sorted.bam
# samtools index sim/ce11_denovo_test_rsem/expr_ngs.transcript.sorted.bam

# RSEM do not support TGS
# rsem-calculate-expression \
#     -p 40 \
#     sim/ce11_denovo_test_pbsim.fq \
#     sim/ce11_denovo_test_rsem/ref \
#     sim/ce11_denovo_test_rsem/expr_tgs
# samtools sort \
#     -@ 40 \
#     sim/ce11_denovo_test_rsem/expr_tgs.transcript.bam \
#     -o sim/ce11_denovo_test_rsem/expr_tgs.transcript.sorted.bam
# samtools index sim/ce11_denovo_test_rsem/expr_tgs.transcript.sorted.bam

# ref-eval-estimate-true-assembly \
#     --reference sim/ce11_denovo_test_rsem/ref \
#     --expression sim/ce11_denovo_test_rsem/expr \
#     --assembly sim/ce11_denovo_test_rsem/estimated_best_assembly \
#     --alignment-policy best

# rsem-prepare-reference \
#     -p 40 \
#     --bowtie \
#     sim/ce11_denovo_test_rsem/estimated_best_assembly_0.fa \
#     sim/ce11_denovo_test_rsem/estimated_best_assembly_ref
# rsem-calculate-expression \
#     -p 40 \
#     --paired-end \
#     sim/ce11_denovo_test_art_1.fq \
#     sim/ce11_denovo_test_art_2.fq \
#     sim/ce11_denovo_test_rsem/estimated_best_assembly_ref \
#     sim/ce11_denovo_test_rsem/estimated_best_assembly_expr
# samtools sort \
#     -@ 40 \
#     sim/ce11_denovo_test_rsem/estimated_best_assembly_expr.transcript.bam \
#     -o sim/ce11_denovo_test_rsem/estimated_best_assembly_expr.transcript.sorted.bam
# samtools index sim/ce11_denovo_test_rsem/estimated_best_assembly_expr.transcript.sorted.bam

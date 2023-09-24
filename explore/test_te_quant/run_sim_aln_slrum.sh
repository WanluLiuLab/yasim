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
    --genomeDir ref/ce11_STAR_idx \
    --readFilesIn sim/ce11_denovo_test_art_1.fq sim/ce11_denovo_test_art_2.fq \
    --outSAMtype BAM Unsorted \
    --outSAMattributes All \
    --outFileNamePrefix aln/ce11_denovo_test_art_STAR/ \
    --twopassMode Basic
samtools sort \
    -@ 40 \
    -o aln/ce11_denovo_test_art.bam \
    aln/ce11_denovo_test_art_STAR/Aligned.out.bam
samtools index aln/ce11_denovo_test_art.bam

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

# minimap2 \
#     -a ref/ce11.te_loci.fa sim/ce11_denovo_test.fa \
#     > aln/ce11_denovo_test.minimap2.te_loci.sam
# samtools sort -@40 \
#     aln/ce11_denovo_test.minimap2.te_loci.sam \
#     -o aln/ce11_denovo_test.minimap2.te_loci.bam
# samtools index aln/ce11_denovo_test.minimap2.te_loci.bam

# pblat \
#     -t=dna \
#     -q=dna \
#     -threads=40 \
#     -out=blast9 \
#     ref/ce11.te_loci.2bit \
#     sim/ce11_denovo_test.fa \
#     aln/ce11_denovo_test.blat.te_loci.blast9

# minimap2 \
#     -a ref/teidx.pkl.xz.fa sim/ce11_denovo_test.fa \
#     > aln/ce11_denovo_test.minimap2.te_consensus.sam
# samtools sort -@40 \
#     aln/ce11_denovo_test.minimap2.te_consensus.sam \
#     -o aln/ce11_denovo_test.minimap2.te_consensus.bam
# samtools index aln/ce11_denovo_test.minimap2.te_consensus.bam

# pblat \
#     -t=dna \
#     -q=dna \
#     -threads=40 \
#     -out=blast9 \
#     ref/teidx.2bit \
#     sim/ce11_denovo_test.fa \
#     aln/ce11_denovo_test.blat.te_consensus.blast9

# perl \
#     ../../deps/RepeatMasker/RepeatMasker \
#     -species "Caenorhabditis elegans" \
#     -engine rmblast \
#     -parallel 40 \
#     -dir aln/ce11_denovo_test_rmsk_prep.rmsk_rmblast.d \
#     -gff -small -s \
#     sim/ce11_denovo_test_rmsk_prep.fa
perl \
    ../../deps/RepeatMasker/RepeatMasker \
    -species "Caenorhabditis elegans" \
    -engine hmmer \
    -parallel 40 \
    -dir aln/ce11_denovo_test_rmsk_prep.rmsk_hmmer.d \
    -gff -small -s \
    sim/ce11_denovo_test_rmsk_prep.fa

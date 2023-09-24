#!/usr/bin/env bash
#SBATCH --job-name=run_sim_rmsk
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=60
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=36:00:00
set -e

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"

conda activate yasim_dev

# for name in ce11.rmsk_loci ce11.rmsk_consensus; do
#     minimap2 \
#         -t 40 \
#         -a \
#         ref/"${name}".fa sim/ce11_denovo_test.fa \
#         > aln/"${name}".minimap2.sam
#     samtools sort \
#         -@40 \
#         -o aln/"${name}".minimap2.bam \
#         aln/"${name}".minimap2.sam
#     samtools index aln/"${name}".minimap2.bam

#     pblat \
#         -t=dna \
#         -q=dna \
#         -threads=40 \
#         -out=blast9 \
#         ref/"${name}".2bit \
#         sim/ce11_denovo_test.fa \
#         aln/"${name}".blat.blast9
# done

# minimap2 \
#     -t 40 \
#     -a \
#     ref/ce11.fa sim/ce11_denovo_test.fa \
#     > aln/ce11.minimap2.sam
# samtools sort \
#     -@40 \
#     -o aln/ce11.minimap2.bam \
#     aln/ce11.minimap2.sam
# samtools index aln/ce11.minimap2.bam
# mkdir -p aln/ce11.minimap2.fc.tsv.d
# featureCounts -L -O -M --primary --ignoreDup \
#     -a ref/ce11.rmsk_filtered.gtf \
#     -g transcript_id \
#     -t dispersed_repeat \
#     -o aln/ce11.minimap2.fc.tsv \
#     -R CORE \
#     --Rpath aln/ce11.minimap2.fc.tsv.d \
#     aln/ce11.minimap2.bam

# perl \
#     ../../deps/RepeatMasker/RepeatMasker \
#     -species "Caenorhabditis elegans" \
#     -engine rmblast \
#     -parallel 40 \
#     -dir aln/rmsk_rmblast.d \
#     -gff -small -s \
#     sim/ce11_denovo_test.fa
# perl \
#     ../../deps/RepeatMasker/RepeatMasker \
#     -species "Caenorhabditis elegans" \
#     -engine hmmer \
#     -parallel 40 \
#     -dir aln/rmsk_hmmer.d \
#     -gff -small -s \
#     sim/ce11_denovo_test.fa

# mkdir sim/ce11_denovo_test_split
# split --lines=16384 sim/ce11_denovo_test.fa sim/ce11_denovo_test_split/

for fn in sim/ce11_denovo_test_split/*; do
    echo "NHMMER ${fn}"
    nhmmer \
        --tblout "${fn}.nhmmer_original.out" \
        -o /dev/null \
        --notextw \
        --noali \
        --dna \
        --cpu 60 \
        ref/ce11.rmsk_consensus.hmm \
        "${fn}"
done

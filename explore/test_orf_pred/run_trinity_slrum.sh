#!/usr/bin/env bash
#SBATCH --job-name=run_trinity
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

singularity run ../../singularity/trinity.sif Trinity \
    --max_memory 100G \
    --verbose \
    --CPU 40 \
    --left sim/WBcel235_art_1.fq \
    --right sim/WBcel235_art_2.fq \
    --seqType fq \
    --output assmb/trinity_ngs_de_novo &>> \
    assmb/trinity_ngs_de_novo.log

singularity run ../../singularity/trinity.sif Trinity \
    --max_memory 100G \
    --verbose \
    --CPU 40 \
    --left sim/WBcel235_art_1.fq \
    --right sim/WBcel235_art_2.fq \
    --long_reads sim/WBcel235_pbsim.fa \
    --seqType fq \
    --output assmb/trinity_ngs_tgs_de_novo &>> \
    assmb/trinity_ngs_tgs_de_novo.log

singularity run ../../singularity/trinity.sif Trinity \
    --max_memory 100G \
    --verbose \
    --CPU 40 \
    --genome_guided_bam aln/WBcel235_art.bam \
    --genome_guided_max_intron 10000 \
    --seqType fq \
    --output assmb/trinity_ngs_referenced &>> \
    assmb/trinity_ngs_referenced.log

singularity run ../../singularity/trinity.sif Trinity \
    --max_memory 100G \
    --verbose \
    --CPU 40 \
    --genome_guided_bam aln/WBcel235_art.bam \
    --long_reads_bam aln/WBcel235_pbsim.bam \
    --genome_guided_max_intron 10000 \
    --seqType fq \
    --output assmb/trinity_ngs_tgs_referenced &>> \
    assmb/trinity_ngs_tgs_referenced.log

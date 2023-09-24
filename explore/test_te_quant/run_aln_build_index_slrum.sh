#!/usr/bin/env bash
#SBATCH --job-name=run_aln_build_index
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=4GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"

conda activate yasim_dev
# genomeSAindexNbases as-is adviced by STAR
STAR --runThreadN 40 \
    --runMode genomeGenerate \
    --genomeDir ref/ce11_STAR_idx \
    --genomeFastaFiles ref/ce11.fa \
    --sjdbGTFfile ref/ce11.ncbiRefSeq.gtf \
    --genomeSAindexNbases 12

mkdir -p ref/teidx_bwa_idx
faToTwoBit ref/teidx.pkl.xz.fa ref/teidx.2bit
bwa index -p ref/teidx_bwa_idx/teidx ref/teidx.pkl.xz.fa
mkdir -p ref/te_loci_bwa_idx
faToTwoBit ref/ce11.te_loci.fa ref/ce11.te_loci.2bit
bwa index -p ref/te_loci_bwa_idx/ce11.te_loci ref/ce11.te_loci.fa

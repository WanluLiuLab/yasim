#!/usr/bin/env bash
# shellcheck disable=SC2002
set -ue
mkdir -p ref aln
# Get Reference
cd ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
cat ce11.ncbiRefSeq.gtf.gz | pigz -d >ce11.ncbiRefSeq.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
cat ce11.fa.gz | pigz -d >ce11.fa
samtools faidx ce11.fa

wget https://dfam.org/releases/current/families/Dfam_curatedonly.h5.gz
cat Dfam_curatedonly.h5.gz | pigz -d >Dfam_curatedonly.h5

sbatch <slrum/run_rmsk_on_reference.sh

python -m yasim_scripts rmsk_filter \
    -i ref/ce11.rmsk_rmblast.d/ce11.fa.out.gff \
    -o ref/ce11.rmsk_filtered.gtf
python -m yasim_scripts rmsk_transcribe \
    -i ref/ce11.rmsk_filtered.gtf \
    -f ref/ce11.fa \
    -o ref/ce11.rmsk_loci.fa
python -m yasim generate_te_index \
    -i ref/Dfam_curatedonly.h5 \
    -o ref/teidx.pkl.xz \
    --dst_consensus_fa_path ref/ce11.rmsk_consensus.fa \
    --dst_hmm_path ref/ce11.rmsk_consensus.hmm \
    --txid 6239
samtools faidx ref/ce11.rmsk_loci.fa
samtools faidx ref/ce11.rmsk_consensus.fa

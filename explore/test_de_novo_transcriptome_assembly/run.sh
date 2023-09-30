#!/usr/bin/env bash
set -ue
mkdir -p logs ref sim docs assmb_final
singularity run ../../singularity/trinity.sif Trinity --help >docs/trinity.help.txt
singularity run ../../singularity/trinity.sif Trinity --show_full_usage_info >docs/trinity.help_full.txt

conda activate yasim-transabyss
transabyss --help >docs/transabyss.help.txt
conda deactivate

conda activate yasim-spades
rnaspades.py --help >docs/rnaspades.help.txt
conda deactivate

conda activate yasim-transrate
transrate --help >docs/transrate.help.txt
conda deactivate

# Get Reference
cd ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
cat ce11.ncbiRefSeq.gtf.gz | pigz -d >ce11.ncbiRefSeq.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
cat ce11.fa.gz | pigz -d >ce11.fa
samtools faidx ce11.fa

wget https://dfam.org/releases/current/families/Dfam_curatedonly.h5.gz
python -m yasim generate_te_index \
    -i Dfam_curatedonly.h5.gz \
    -o teidx.pkl.xz \
    -t 6239
cd ..
sbatch <run_aln_build_index_slrum.sh

# Simulating data

cd sim
python -m yasim generate_gene_te_fusion \
    -g ../ref/ce11.ncbiRefSeq.gtf \
    -f ../ref/ce11.fa \
    --tedb ../ref/teidx.pkl.xz \
    -o ce11_denovo_test \
    -d 150 \
    -n 10000
samtools faidx ce11_denovo_test.fa
python -m labw_utils.bioutils split_fasta ce11_denovo_test.fa
cd ..
python plot_n_frags.py

sbatch <run_art_slrum.sh
sbatch <run_pbsim_slrum.sh

seqkit fq2fa sim/ce11_denovo_test_pbsim.fq >sim/ce11_denovo_test_pbsim.fa

# Running Aligner
rm -fr aln
mkdir -p aln
sbatch <run_sim_aln_slrum.sh

# Running Assembler
rm -fr assmb
mkdir -p assmb
sbatch <run_trinity_slrum.sh
sbatch <run_transabyss_slrum.sh
sbatch <run_rnaspades_slrum.sh

rm -fr assmb_final
mkdir -p assmb_final
for fn in \
    assmb/rnaspades_*/transcripts.fasta \
    assmb/trinity_*.Trinity.fasta \
    assmb/*/Trinity-GG.fasta \
    assmb/transabyss_*/transabyss-final.fa; do
    fn_dst="$(echo "${fn}" | sed 's;^assmb/;;' | sed -E 's;^([a-z_]+).*$;\1;')"
    mv "${fn}" assmb_final/"${fn_dst}".fa
done

sbatch <run_transrate_slrum.sh
# sbatch <run_detonate_slrum.sh
python plot_transrate.py

# Run RepeatMasker
python prepare_for_rmsk.py sim/ce11_denovo_test_pbsim
python prepare_for_rmsk.py sim/ce11_denovo_test
sbatch <run_sim_rmsk_slrum.sh
python rmsk_cmp.py
python plot_rmsk.py
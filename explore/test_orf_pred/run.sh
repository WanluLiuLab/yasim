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
python -m labw_utils.bioutils transcribe \
    -f ce11.fa \
    -g ce11.ncbiRefSeq.gtf \
    -o ce11.trans.fa
samtools faidx ce11.trans.fa
cd ..
sbatch <run_aln_build_index_slrum.sh

# Simulating data

cd sim
python -m yasim generate_gene_depth \
    -g ../ref/ce11.ncbiRefSeq.gtf \
    -o ce11_denovo_test.gene_depth.tsv \
    -d 150
python -m yasim generate_isoform_depth \
    -g ../ref/ce11.ncbiRefSeq.gtf \
    -o ce11_denovo_test.isoform_depth.tsv \
    -d ce11_denovo_test.gene_depth.tsv
cd ..

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

# Copy generated files
rm -fr assmb_transcripts_final
mkdir -p assmb_transcripts_final
for fn in \
    assmb/rnaspades_*/transcripts.fasta \
    assmb/trinity_*.Trinity.fasta \
    assmb/*/Trinity-GG.fasta \
    assmb/transabyss_*/transabyss-final.fa; do
    fn_dst="$(echo "${fn}" | sed 's;^assmb/;;' | sed -E 's;^([a-z_]+).*$;\1;')"
    mv "${fn}" assmb_transcripts_final/"${fn_dst}".fa
done

# Running CDS Predictor
...

# Exit

sbatch <run_transrate_slrum.sh
python plot_transrate.py

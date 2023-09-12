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

wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz -O WBcel235.pep.faa.gz
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz -O WBcel235.genome.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.57.gtf.gz -O WBcel235.gtf.gz
pigz -cdf - < WBcel235.pep.faa.gz > WBcel235.pep.faa
pigz -cdf - < WBcel235.genome.fa.gz > WBcel235.genome.fa
pigz -cdf - < WBcel235.gtf.gz > WBcel235.gtf

python -m labw_utils.bioutils transcribe \
    -f WBcel235.genome.fa \
    -g WBcel235.gtf \
    -o WBcel235.trans.fa
samtools faidx WBcel235.genome.fa
samtools faidx WBcel235.pep.faa
samtools faidx WBcel235.trans.fa
cd ..
sbatch <run_aln_build_index_slrum.sh

# Simulating data

cd sim
python -m yasim generate_gene_depth \
    -g ../ref/WBcel235.gtf \
    -o WBcel235.gene_depth.tsv \
    -d 150
python -m yasim generate_isoform_depth \
    -g ../ref/WBcel235.gtf \
    -o WBcel235.isoform_depth.tsv \
    -d WBcel235.gene_depth.tsv
cd ..

sbatch <run_art_slrum.sh
sbatch <run_pbsim_slrum.sh

seqkit fq2fa sim/WBcel235_pbsim.fq >sim/WBcel235_pbsim.fa

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

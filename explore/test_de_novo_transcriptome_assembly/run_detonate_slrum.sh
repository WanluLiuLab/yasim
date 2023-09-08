#!/usr/bin/env bash
#SBATCH --job-name=run_detonate
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=40
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=20:00:00

echo "No longer used."
exit 1

eval "$("${HOME}"/conda/condabin/conda shell.bash hook)"
conda activate yasim-detonate
set -e

mkdir -p assmb_final/detonate/
NREADS="$(wc -l sim/ce11_denovo_test_art_1.fq | cut -f 1 -d ' ')"
EBA_PATH="sim/ce11_denovo_test_rsem/estimated_best_assembly_0.fa"

printf 'SAMPLE\tKC\n' >assmb_final/detonate/kc.txt
printf 'SAMPLE\tCONTIG\n' >assmb_final/detonate/contig.txt
printf 'SAMPLE\tNUCL\n' >assmb_final/detonate/nucl.txt
for fn in assmb_final/*.fa; do
    ref-eval \
        --scores kc \
        --A-seqs "${fn}" \
        --B-seqs "${EBA_PATH}" \
        --B-expr sim/ce11_denovo_test_rsem/estimated_best_assembly_expr.isoforms.results \
        --kmerlen 150 \
        --readlen 150 \
        --num-reads "${NREADS}" |
        awk '$1 == "kmer_compression_score"' |
        sed 's;kmer_compression_score;'"${fn}"';' \
            >>assmb_final/detonate/kc.txt

    pblat \
        -threads=40 \
        -minIdentity=80 \
        "${EBA_PATH}" \
        "${fn}" \
        assmb_final/detonate/eba_to_fn.psl
    pblat \
        -threads=40 \
        -minIdentity=80 \
        "${fn}" \
        "${EBA_PATH}" \
        assmb_final/detonate/fn_to_eba.psl
    ref-eval \
        --scores contig,nucl \
        --weighted no \
        --A-seqs "${fn}" \
        --B-seqs "${EBA_PATH}" \
        --A-to-B assmb_final/detonate/eba_to_fn.psl \
        --B-to-A assmb_final/detonate/fn_to_eba.psl \
        --min-frac-identity 0.90 \
        >assmb_final/detonate/tmp_contig_nucl.txt
    cat assmb_final/detonate/tmp_contig_nucl.txt |
        awk '$1 == "unweighted_contig_F1"' |
        sed 's;unweighted_contig_F1;'"${fn}"';' \
            >>assmb_final/detonate/contig.txt
    cat assmb_final/detonate/tmp_contig_nucl.txt |
        awk '$1 == "unweighted_nucl_F1"' |
        sed 's;unweighted_nucl_F1;'"${fn}"';' \
            >>assmb_final/detonate/nucl.txt
    rm -f \
        assmb_final/detonate/eba_to_fn.psl \
        assmb_final/detonate/fn_to_eba.psl \
        assmb_final/detonate/tmp_contig_nucl.txt
done

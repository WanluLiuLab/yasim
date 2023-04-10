#!/usr/bin/env bash
eval "$(conda 'shell.bash' 'hook' 2>/dev/null)"
conda activate yasim_c_elegans_as_depth_analysis

set -ue
export PYTHONPATH="../src"
# Install spladder
bash common/get_spladder.sh
bash common/get_index_c_elegans.sh
bash batch_c_elegans/get_data.sh

for fastq in ...; do
    bash common/procedure_ngs_paired_end.sh --FASTQ_BASE_NAME:"${fastq}" --GENE_REFERENCE:ce11.fa --TRANSCRIPT_REFERENCE:ce11.trans.fa --THREAD_NUM:36 --GENE_GTF:ce11.ncbiRefSeq.gtf
done

bash common/procedure_ngs_single_end.sh
bash common/procedure_tgs_single_end.sh

Rscript common/R/summarize_as.R --libfile common/R/lib.R --filename batch_c_elegans/all_as.tsv

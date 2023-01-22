#shellcheck shell=bash
SPLADDR_DIR="${FASTQ_BASE_NAME}"_SPLADDER

[ -d "${SPLADDR_DIR}" ] || \
spladder build \
-b "${FASTQ_BASE_NAME}".GENE.bam \
-o "${FASTQ_BASE_NAME}"_SPLADDER \
-a "${GENE_GTF}" \
--parallel "${THREAD_NUM}"

[ -f "${SPLADDR_DIR}/all_as.tsv" ] || \
    python3 "${SHDIR}"/src/extract_as_from_spladder.py \
    --spladder_dir "${SPLADDR_DIR}"

[ -f "${SPLADDR_DIR}/all_as.tsv.aggregated.tsv" ] || \
    Rscript "${SHDIR}"/R/summarize_as.R \
    --libfile "${SHDIR}"/R/lib.R \
    --filename "${SPLADDR_DIR}/all_as.tsv"

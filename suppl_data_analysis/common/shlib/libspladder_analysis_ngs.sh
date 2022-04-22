#shellcheck shell=bash
[ ! -d "${FASTQ_BASE_NAME}"_SPLADDER ] && \
spladder build \
-b "${FASTQ_BASE_NAME}".GENE.bam \
-o "${FASTQ_BASE_NAME}"_SPLADDER \
-a "${GENE_GTF}" \
--parallel "${THREAD_NUM}"

true # FIXME: Why have to add this?

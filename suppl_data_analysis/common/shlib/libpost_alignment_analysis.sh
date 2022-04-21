#shellcheck shell=bash
[ ! -f "${FASTQ_BASE_NAME}".GENE.bam.bai ] && samtools index "${FASTQ_BASE_NAME}".GENE.bam
[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam.bai ] && samtools index "${FASTQ_BASE_NAME}".TRANS.bam
[ ! -f "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv ] && samtools depth "${FASTQ_BASE_NAME}".TRANS.bam > "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv

[ ! -f "${FASTQ_BASE_NAME}".transcript.depth.tsv ] &&  \
    Rscript "${SHDIR}"/R/transform_depth_results.R \
    --libfile "${SHDIR}"/R/lib.R \
    --input "${FASTQ_BASE_NAME}".TRANS.bam.depth.tsv \
    --fa_stats "${TRANSCRIPT_REFERENCE}.stats" \
    --output "${FASTQ_BASE_NAME}".transcript.depth.tsv

[ ! -f "${FASTQ_BASE_NAME}".stringtie_unguided.gtf ] && \
stringtie "${FASTQ_BASE_NAME}".GENE.bam \
-p "${THREAD_NUM}" \
-o "${FASTQ_BASE_NAME}".stringtie_unguided.gtf

[ ! -f "${FASTQ_BASE_NAME}".stringtie_guided.gtf ] && \
stringtie "${FASTQ_BASE_NAME}".GENE.bam \
-p "${THREAD_NUM}" \
-G "${GENE_GTF}" \
-o "${FASTQ_BASE_NAME}".stringtie_guided.gtf

true

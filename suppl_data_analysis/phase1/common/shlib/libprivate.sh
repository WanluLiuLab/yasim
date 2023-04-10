# shellcheck shell=bash
# shellcheck disable=SC2034
VERSION=0.1

declare -a STDS

for opt in "${@}"; do
    if isopt "${opt}"; then
        case "${opt}" in
        "-h" | "--help")
            echo "
Usage: bash ${0} \
    --FASTQ_BASE_NAME:FASTQ_BASE_NAME \
    --GENE_REFERENCE:GENE_REFERENCE \
    --TRANSCRIPT_REFERENCE:TRANSCRIPT_REFERENCE \
    --THREAD_NUM:THREAD_NUM \
    --GENE_GTF:GENE_GTF

Where:
    - FASTQ_BASE_NAME is basename of a fastq file. Can be compressed by GNU GZip.
    If is pair_end, please ignore '_1' or '_2' suffix.
    E.g., SRR8717004_1.fastq.gz -> SRR8717004

    - GENE_REFERENCE is genome reference in FASTA format. Cannot be compressed. Should be indexed.

    - TRANSCRIPT_REFERENCE is transcriptome reference in FASTA format. Cannot be compressed. Should be indexed.

    - THREAD_NUM is number of thread to use.

    - GENE_GTF is GTF file of genes.

    All arguments are required.
"
            builtin exit 0
            ;;
        "-v" | "--version")
            infoh "Version ${VERSION}, compatiable with libdo Version 2 & 3"
            infoh "minimap2:"
            minimap2 --version
            infoh "samtools:"
            samtools version
            infoh "stringtie:"
            stringtie --version
            builtin exit 0
            ;;
        --FASTQ_BASE_NAME\:*)
            FASTQ_BASE_NAME=${opt:18}
            ;;
        --GENE_REFERENCE\:*)
            GENE_REFERENCE=${opt:17}
            ;;
        --TRANSCRIPT_REFERENCE\:*)
            TRANSCRIPT_REFERENCE=${opt:23}
            ;;
        --THREAD_NUM\:*)
            THREAD_NUM=${opt:13}
            ;;
        --GENE_GTF\:*)
            GENE_GTF=${opt:11}
            ;;
        "-x")
            builtin set -x
            ;;
        "-V" | "--verbose")
            builtin set -v
            ;;
        *)
            warnh "Option '${opt}' invalid. Ignored"
            ;;
        esac
    else
        STDS=("${STDS[@]}" "${opt}")
    fi
done

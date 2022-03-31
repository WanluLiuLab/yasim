# shellcheck shell=bash
# shellcheck disable=SC2034
VERSION=0.1

declare -a STDS

for opt in "${@}"; do
    if isopt "${opt}"; then
        case "${opt}" in
        "-h" | "--help")
            echo "
Usage: ${0} --FASTQ_NAME:FASTQ_NAME --GENE_REFERENCE:GENE_REFERENCE --TRANSCRIPT_REFERENCE:TRANSCRIPT_REFERENCE --THREAD_NUM:THREAD_NUM --GENE_GTF:GENE_GTF
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
        --FASTQ_NAME\:*)
            FASTQ_NAME=${opt:13}
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
            GENE_GTF=${opt:10}
            ;;
        "-x")
            builtin set -x
            ;;
        *)
            warnh "Option '${opt}' invalid. Ignored"
            ;;
        esac
    else
        STDS=("${STDS[@]}" "${opt}")
    fi
done

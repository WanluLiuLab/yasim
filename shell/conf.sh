#shell=bash
# shellcheck disable=SC2034
set -ue
ROOTDIR="$(dirname "${SHDIR}")"
cd "${ROOTDIR}" || exit 1
export PYTHONPATH="${ROOTDIR}/src:${PYTHONPATH:-}"

REFERENCE_FASTA=/home/yuzj/Desktop/BioRef/chr1.fa
REFERENCE_GTF=/home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq.gtf
SALMON_INDEX=/home/yuzj/Desktop/BioRef/chr1_salmon_index
THREADS=40
SEL_GTF=Homo_sapiens.GRCh38.105_chr1_sel.gtf
SEL_CDNA_FASTA=hg38.chr1.transcripts_sel.fa
DEPTH_TSV=hg38.chr1.gene.depth.tsv
FASTQ_BASENAME=hg38.chr1.transcripts
STAR_INDEX=/home/yuzj/Desktop/BioRef/chr1_STAR_index

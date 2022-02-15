
export ROOTDIR := $(abspath $(CURDIR)/../)
export PYTHONPATH := "$(ROOTDIR)/src:$(PYTHONPATH:-)"
export YASIM := PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python -m yasim
export REFERENCE_FASTA := /home/yuzj/Desktop/BioRef/chr1.fa
export REFERENCE_GTF := /home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq.gtf
export REFERENCE_CDNA := $(ROOTDIR)/hg38.chr1.reference_transcripts.fa
export THREADS := 40
export SEL_GTF := $(ROOTDIR)/Homo_sapiens.GRCh38.105_chr1_sel.gtf
export SEL_CDNA_FASTA := $(ROOTDIR)/hg38.chr1.transcripts_sel.fa
export DEPTH_TSV := $(ROOTDIR)/hg38.chr1.gene.depth.tsv
export FASTQ_BASENAME := $(ROOTDIR)/hg38.chr1.transcripts

# export R_ENVIRON_USER := $(ROOTDIR)/R/.Rprofile

.SECONDARY:

all: $(ROOTDIR)/yasim_to_salmon_quant_dwgsim.png \
$(ROOTDIR)/yasim_to_featureCounts_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_htseq_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_e_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_e_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_featureCounts_quant_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_htseq_quant_hisat2_dwgsim.png

include $(CURDIR)/general.mk
include $(CURDIR)/hisat2.mk
include $(CURDIR)/STAR.mk
include $(CURDIR)/salmon.mk
include $(CURDIR)/featureCounts.mk
include $(CURDIR)/htseq.mk
include $(CURDIR)/stringtie.mk

$(SEL_GTF):
	$(YASIM) alternative_splicing -f "$(REFERENCE_FASTA)" -g "$(REFERENCE_GTF)" -o $@

$(SEL_CDNA_FASTA): $(SEL_GTF)
	$(YASIM) transcribe -f "$(REFERENCE_FASTA)" -g "$(SEL_GTF)" -o $@

$(DEPTH_TSV): $(SEL_GTF)
	$(YASIM) dge -g "$(SEL_GTF)" -o "$(DEPTH_TSV)" -d 100 -l 20

$(REFERENCE_CDNA):
	$(YASIM) transcribe -g "$(REFERENCE_GTF)" -o $@ -f "$(REFERENCE_FASTA)"

$(FASTQ_BASENAME)_dwgsim_1.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) dwgsim -F "$(SEL_CDNA_FASTA)" -o "$(FASTQ_BASENAME)"_dwgsim -d "$(DEPTH_TSV)"


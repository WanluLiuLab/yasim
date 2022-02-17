
export ROOTDIR := $(abspath $(CURDIR)/../)
export DATADIR := $(ROOTDIR)/iter_2
export PYTHONPATH := "$(ROOTDIR)/src:$(PYTHONPATH:-)"
export YASIM := PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python -m yasim
export REFERENCE_FASTA := /home/yuzj/Desktop/BioRef/chr1.fa
export REFERENCE_GTF := /home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq.gtf
export REFERENCE_CDNA := $(DATADIR)/hg38.chr1.reference_transcripts.fa
export THREADS := 40
export SEL_GTF := $(DATADIR)/Homo_sapiens.GRCh38.105_chr1_sel.gtf
export SEL_CDNA_FASTA := $(DATADIR)/hg38.chr1.transcripts_sel.fa
export SEL_CDNA_FASTA_D := $(SEL_CDNA_FASTA).d
export DEPTH_TSV := $(DATADIR)/hg38.chr1.gene.depth.tsv
export FASTQ_BASENAME := $(DATADIR)/hg38.chr1.transcripts

ifeq "$(wildcard $(DATADIR) )" ""
$(shell mkdir -p $(DATADIR))
endif

# export R_ENVIRON_USER := $(ROOTDIR)/R/.Rprofile

.SECONDARY:
.DEFAULT: all

include $(CURDIR)/general.mk
include $(CURDIR)/simulators.mk
include $(CURDIR)/hisat2.mk
include $(CURDIR)/STAR.mk
include $(CURDIR)/minimap2.mk
include $(CURDIR)/salmon.mk
include $(CURDIR)/featureCounts.mk
include $(CURDIR)/htseq.mk
include $(CURDIR)/stringtie.mk
include $(CURDIR)/final_targets.mk

.PHONY: dist
dist:
	$(MAKE) -C $(ROOTDIR)/dist

$(SEL_GTF):
	$(YASIM) alternative_splicing -f "$(REFERENCE_FASTA)" -g "$(REFERENCE_GTF)" -o $@

$(DEPTH_TSV): $(SEL_GTF)
	$(YASIM) dge -g "$(SEL_GTF)" -o "$(DEPTH_TSV)" -d 100 -l 20

$(SEL_CDNA_FASTA_D): $(SEL_GTF) $(DEPTH_TSV)
	$(YASIM) transcribe -f "$(REFERENCE_FASTA)" -g "$(SEL_GTF)" -d $(DEPTH_TSV) -o $(SEL_CDNA_FASTA)

$(REFERENCE_CDNA):
	$(YASIM) transcribe -g "$(REFERENCE_GTF)" -o $@ -f "$(REFERENCE_FASTA)"

.PHONY: nothing
nothing:
	$(message Nothing made)

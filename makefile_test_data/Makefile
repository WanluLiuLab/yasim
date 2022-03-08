
export ROOTDIR := $(abspath $(CURDIR)/../)
export DATADIR := $(ROOTDIR)/iter_5
export PYTHONPATH := "$(ROOTDIR)/src:$(PYTHONPATH:-)"
export YASIM := PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python -m yasim
export YASIM_SCRIPTS := PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python -m yasim_scripts
export REFERENCE_FASTA := $(DATADIR)/ce11.fa
export REFERENCE_GTF := $(DATADIR)/ce11.ncbiRefSeq.gtf
export REFERENCE_CDNA := $(DATADIR)/ce11.reference_transcripts.fa
export THREADS := 40
export SEL_GTF := $(DATADIR)/ce11_sel.gtf
export SEL_CDNA_FASTA := $(DATADIR)/ce11_sel.fa
export SEL_CDNA_FASTA_D := $(SEL_CDNA_FASTA).d
export DEPTH_TSV := $(DATADIR)/ce11.gene.depth.tsv
export FASTQ_BASENAME := $(DATADIR)/ce11.transcripts

$(DATADIR)/ce11.fa.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz -o $@

$(DATADIR)/ce11.ncbiRefSeq.gtf.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz -o $@

ifeq "$(wildcard $(DATADIR) )" ""
$(shell mkdir -p $(DATADIR))
endif

# export R_ENVIRON_USER := $(ROOTDIR)/R/.Rprofile

.SECONDARY:
.DEFAULT_GOAL: all

include $(CURDIR)/final_targets.mk

include $(CURDIR)/general.mk
include $(CURDIR)/simulators.mk
include $(CURDIR)/hisat2.mk
# include $(CURDIR)/STAR.mk
include $(CURDIR)/minimap2.mk
include $(CURDIR)/salmon.mk
include $(CURDIR)/featureCounts.mk
# include $(CURDIR)/htseq.mk
include $(CURDIR)/stringtie.mk

$(SEL_GTF):$(REFERENCE_FASTA) $(REFERENCE_GTF)
	$(YASIM) alternative_splicing -f "$(REFERENCE_FASTA)" -g "$(REFERENCE_GTF)" -o $@

$(DEPTH_TSV): $(SEL_GTF)
	$(YASIM) dge -g "$(SEL_GTF)" -o "$(DEPTH_TSV)" -d 100

$(SEL_CDNA_FASTA): $(SEL_GTF) $(DEPTH_TSV) $(REFERENCE_FASTA)
	$(YASIM) transcribe -f "$(REFERENCE_FASTA)" -g "$(SEL_GTF)" -d $(DEPTH_TSV) -o $(SEL_CDNA_FASTA)

$(SEL_CDNA_FASTA).stats:$(SEL_CDNA_FASTA)

$(REFERENCE_CDNA): $(REFERENCE_FASTA) $(REFERENCE_GTF)
	$(YASIM) transcribe -g "$(REFERENCE_GTF)" -o $@ -f "$(REFERENCE_FASTA)"

.PHONY: nothing
nothing:
	$(message Nothing made)



$(FASTQ_BASENAME)_%_ground_truth.tsv:$(FASTQ_BASENAME)_dwgsim.fq.stats $(SEL_CDNA_FASTA).stats
	Rscript $(ROOTDIR)/R/generate_ground_truth.R --libfile "$(ROOTDIR)/R/lib.R" --fq_stats $(FASTQ_BASENAME)_dwgsim.fq.stats --fa_stats $(SEL_CDNA_FASTA).stats --output  $@



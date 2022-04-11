

export ALL_BAM_TGS_FILENAME := $(DATADIR)/minimap2_pbsim2_R94.bam \
$(DATADIR)/minimap2_badread_pacbio2016.bam \
$(DATADIR)/minimap2_badread_verybad.bam \
$(DATADIR)/minimap2_badread_verynice.bam \
$(DATADIR)/minimap2_pbsim2_R103.bam \
$(DATADIR)/minimap2_badread_nanopore2020.bam \
$(DATADIR)/minimap2_pbsim2_P5C3.bam \
$(DATADIR)/minimap2_badread_nanopore2018.bam \
$(DATADIR)/minimap2_pbsim_clr.bam \
$(DATADIR)/minimap2_pbsim_ccs.bam \
$(DATADIR)/minimap2_pbsim2_R95.bam \
$(DATADIR)/minimap2_pbsim2_P6C4.bam \
$(DATADIR)/minimap2_pbsim2_P4C2.bam \
# $(DATADIR)/minimap2_simlord.bam \

export ALL_BAM_NGS_FILENAME := $(DATADIR)/hisat2_dwgsim.bam \
$(DATADIR)/STAR_dwgsim.bam \

ALL_BAM_TGS: $(ALL_BAM_TGS_FILENAME)
ALL_BAM_NGS: $(ALL_BAM_NGS_FILENAME)

export BAM_FIGS_TGS_FIG_FILENAME := $(addsuffix .LEN.png,$(ALL_BAM_TGS_FILENAME))
export BAM_FIGS_NGS_FIG_FILENAME := $(addsuffix .LEN.png,$(ALL_BAM_NGS_FILENAME))
export BAM_FIGS_TGS_BW_FILENAME := $(addsuffix .bw,$(ALL_BAM_TGS_FILENAME))
export BAM_FIGS_NGS_BW_FILENAME := $(addsuffix .bw,$(ALL_BAM_NGS_FILENAME))


.PHONY:all

all:

BAM_FIGS: BAM_FIGS_TGS BAM_FIGS_NGS
BAM_BW: BAM_BW_TGS BAM_BW_NGS

$(DATADIR)/%.bam.tsv:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai
	$(YASIM_SCRIPTS) get_sam_satistics --sam $< --out $@

$(DATADIR)/%.bam.tsv:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai
	$(YASIM_SCRIPTS) get_sam_satistics --sam $< --out $@

$(DATADIR)/%.bam.LEN.png: $(DATADIR)/%.bam.tsv
	Rscript $(ROOTDIR)/R/test_bam_satistics.R --libfile "$(ROOTDIR)/R/lib.R" --ss_tsv $< --output  $(basename $(basename $@))

$(DATADIR)/%.bam.bw:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai
	bamCoverage -b $< -o $@

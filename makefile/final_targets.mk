

export ALL_BAM_TGS_FILENAME := $(DATADIR)/minimap2_pbsim2_R94.bam \
$(DATADIR)/minimap2_badread_pacbio2016.bam \
$(DATADIR)/minimap2_pbsim2_R103.bam \
$(DATADIR)/minimap2_badread_nanopore2020.bam \
$(DATADIR)/minimap2_pbsim2_P5C3.bam \
$(DATADIR)/minimap2_badread_nanopore2018.bam \
$(DATADIR)/minimap2_pbsim_clr.bam \
$(DATADIR)/minimap2_pbsim_ccs.bam \
$(DATADIR)/minimap2_pbsim2_R95.bam \
$(DATADIR)/minimap2_pbsim2_P6C4.bam \
$(DATADIR)/minimap2_simlord.bam \
$(DATADIR)/minimap2_pbsim2_P4C2.bam \

export ALL_BAM_NGS_FILENAME := $(DATADIR)/hisat2_dwgsim.bam \
$(DATADIR)/STAR_dwgsim.bam \

export BAM_FIGS_TGS_FIG_FILENAME := $(addsuffix .LEN.png,$(ALL_BAM_TGS_FILENAME))
export BAM_FIGS_NGS_FIG_FILENAME := $(addsuffix .LEN.png,$(ALL_BAM_NGS_FILENAME))
export BAM_FIGS_TGS_BW_FILENAME := $(addsuffix .bw,$(ALL_BAM_TGS_FILENAME))
export BAM_FIGS_NGS_BW_FILENAME := $(addsuffix .bw,$(ALL_BAM_NGS_FILENAME))


.PHONY:all DWGSIM PBSIM PBSIM2 BADREAD NGS TGS

all:NGS TGS

NGS: DWGSIM BAM_FIGS_NGS BAM_BW_NGS
TGS: PBSIM PBSIM2 BADREAD SIMLORD BAM_FIGS_TGS BAM_BW_TGS
BAM_FIGS: BAM_FIGS_NGS BAM_FIGS_TGS
BAM_BW: BAM_BW_NGS BAM_BW_TGS

DWGSIM: \
$(DATADIR)/yasim_to_salmon_quant_dwgsim.png \
$(DATADIR)/yasim_to_featureCounts_quant_STAR_dwgsim.png \
$(DATADIR)/yasim_to_htseq_quant_STAR_dwgsim.png \
$(DATADIR)/yasim_to_stringtie_quant_STAR_dwgsim.png \
$(DATADIR)/yasim_to_stringtie_quant_e_STAR_dwgsim.png \
$(DATADIR)/yasim_to_stringtie_quant_hisat2_dwgsim.png \
$(DATADIR)/yasim_to_stringtie_quant_e_hisat2_dwgsim.png \
$(DATADIR)/yasim_to_featureCounts_quant_hisat2_dwgsim.png \
$(DATADIR)/yasim_to_htseq_quant_hisat2_dwgsim.png

BAM_FIGS_NGS: $(BAM_FIGS_NGS_FIG_FILENAME)
BAM_FIGS_TGS: $(BAM_FIGS_TGS_FIG_FILENAME)
BAM_BW_NGS: $(BAM_FIGS_NGS_BW_FILENAME)
BAM_BW_TGS: $(BAM_FIGS_TGS_BW_FILENAME)

# -e is not that good in TGS, so not used.

PBSIM: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim_clr.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim_ccs.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim_clr.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim_ccs.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_clr.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_ccs.png

PBSIM2: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R103.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R95.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R94.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P6C4.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P5C3.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P4C2.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_R103.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_R95.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_R94.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_P6C4.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_P5C3.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_pbsim2_P4C2.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R103.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R95.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R94.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P6C4.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P5C3.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P4C2.png

BADREAD: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_badread_nanopore2018.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_badread_nanopore2020.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_badread_pacbio2016.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_badread_nanopore2018.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_badread_nanopore2020.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_badread_pacbio2016.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2018.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2020.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_pacbio2016.png \

SIMLORD: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_simlord.png \
$(DATADIR)/yasim_to_featureCounts_quant_minimap2_simlord.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_simlord.png \

$(DATADIR)/%.bam.tsv:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai
	$(YASIM) get_sam_satistics --sam $< --out $@

$(DATADIR)/%.bam.LEN.png: $(DATADIR)/%.bam.tsv
	Rscript $(ROOTDIR)/R/test_bam_satistics.R --libfile "$(ROOTDIR)/R/lib.R" --ss_tsv $< --output  $(basename $(basename $@))

$(DATADIR)/%.bam.bw:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai
	bamCoverage -b $< -o $@

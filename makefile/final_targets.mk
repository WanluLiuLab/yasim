.PHONY:all DWGSIM PBSIM PBSIM2 BADREAD NGS TGS

all:NGS TGS

NGS: DWGSIM
TGS: PBSIM PBSIM2 BADREAD

DWGSIM: \
$(ROOTDIR)/yasim_to_salmon_quant_dwgsim.png \
$(ROOTDIR)/yasim_to_featureCounts_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_htseq_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_e_STAR_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_stringtie_quant_e_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_featureCounts_quant_hisat2_dwgsim.png \
$(ROOTDIR)/yasim_to_htseq_quant_hisat2_dwgsim.png

# -e is not that good in TGS, so not used.

PBSIM: \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim_clr.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim_ccs.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_clr.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_ccs.png

PBSIM2: \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R103.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R95.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R94.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P6C4.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P5C3.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P4C2.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R103.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R95.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_R94.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P6C4.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P5C3.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_pbsim2_P4C2.png

BADREAD: \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_badread_nanopore2018.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_badread_nanopore2020.png \
$(ROOTDIR)/yasim_to_stringtie_quant_minimap2_badread_pacbio2016.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2018.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2020.png \
# $(ROOTDIR)/yasim_to_stringtie_quant_e_minimap2_badread_pacbio2016.png \



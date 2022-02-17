.PHONY:all DWGSIM PBSIM PBSIM2 BADREAD NGS TGS

all:NGS TGS

NGS: DWGSIM
TGS: PBSIM PBSIM2 BADREAD

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

# -e is not that good in TGS, so not used.

PBSIM: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim_clr.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim_ccs.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_clr.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_pbsim_ccs.png

PBSIM2: \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R103.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R95.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_R94.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P6C4.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P5C3.png \
$(DATADIR)/yasim_to_stringtie_quant_minimap2_pbsim2_P4C2.png \
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
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2018.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_nanopore2020.png \
# $(DATADIR)/yasim_to_stringtie_quant_e_minimap2_badread_pacbio2016.png \



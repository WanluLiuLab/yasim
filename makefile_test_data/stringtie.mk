
# $(DATADIR)/stringtie_e_STAR_%.gtf: $(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
# 	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

# $(DATADIR)/stringtie_STAR_%.gtf: $(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
# 	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

# $(DATADIR)/stringtie_e_hisat2_%.gtf: $(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
# 	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

# $(DATADIR)/stringtie_hisat2_%.gtf: $(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
# 	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

.PHONY: STRINGTIE_QUANT_ALL
STRINGTIE_QUANT_ALL:$(subst .bam,.tsv,$(subst $(DATADIR)/,$(DATADIR)/stringtie_quant_,$(wildcard $(DATADIR)/*.bam)))

$(DATADIR)/stringtie_quant_%.tsv:$(DATADIR)/stringtie_%.gtf
	$(YASIM_SCRIPTS) parse_stringtie_into_tsv -g $< -o $@

# $(DATADIR)/stringtie_e_minimap2_%.gtf: $(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
# 	stringtie -e -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_minimap2_%.gtf: $(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	stringtie -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

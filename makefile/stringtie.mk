
$(DATADIR)/stringtie_e_STAR_%.gtf: $(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_STAR_%.gtf: $(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_e_hisat2_%.gtf: $(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_hisat2_%.gtf: $(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_quant_%.tsv:$(DATADIR)/stringtie_%.gtf
	$(YASIM_SCRIPTS) parse_stringtie_into_tsv -g $< -o $@

$(DATADIR)/stringtie_e_minimap2_%.gtf: $(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	stringtie -e -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/stringtie_minimap2_%.gtf: $(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	stringtie -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(DATADIR)/yasim_to_stringtie_quant_e_%.png: $(DATADIR)/stringtie_quant_e_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_stringtie.R --libfile "$(ROOTDIR)/R/lib.R" -e --stringtie_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(DATADIR)/yasim_to_stringtie_quant_%.png: $(DATADIR)/stringtie_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_stringtie.R --libfile "$(ROOTDIR)/R/lib.R"  --stringtie_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

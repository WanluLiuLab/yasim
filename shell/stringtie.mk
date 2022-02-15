
$(ROOTDIR)/stringtie_e_STAR_%.gtf: $(ROOTDIR)/STAR_%.bam
	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/stringtie_STAR_%.gtf: $(ROOTDIR)/STAR_%.bam
	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/stringtie_e_hisat2_%.gtf: $(ROOTDIR)/hisat2_%.bam
	stringtie -e -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/stringtie_hisat2_%.gtf: $(ROOTDIR)/hisat2_%.bam
	stringtie -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/stringtie_quant_%.tsv:$(ROOTDIR)/stringtie_%.gtf
	PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python $(ROOTDIR)/src/scripts/parse_stringtie_into_tsv.py -g $< -o $@

$(ROOTDIR)/stringtie_e_minimap2_%.gtf: $(ROOTDIR)/minimap2_%.bam
	stringtie -e -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/stringtie_minimap2_%.gtf: $(ROOTDIR)/minimap2_%.bam
	stringtie -L -G "$(REFERENCE_GTF)" -o $@ -p $(THREADS) $<

$(ROOTDIR)/yasim_to_stringtie_quant_e_%.png: $(ROOTDIR)/stringtie_quant_e_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_stringtie.R --libfile "$(ROOTDIR)/R/lib.R" -e --stringtie_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(ROOTDIR)/yasim_to_stringtie_quant_%.png: $(ROOTDIR)/stringtie_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_stringtie.R --libfile "$(ROOTDIR)/R/lib.R"  --stringtie_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

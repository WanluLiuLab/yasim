# $(DATADIR)/featureCounts_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
# 	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

$(DATADIR)/yasim_to_featureCounts_quant_%.png: $(DATADIR)/featureCounts_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_featureCounts.R --libfile "$(ROOTDIR)/R/lib.R"  --featureCounts_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

# $(DATADIR)/featureCounts_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
# 	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

$(DATADIR)/featureCounts_quant_minimap2_%.tsv:$(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	featureCounts -t transcript -g transcript_id -a "$(REFERENCE_GTF)" -L -o $@ $<

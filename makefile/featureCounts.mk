$(DATADIR)/featureCounts_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

$(DATADIR)/yasim_to_featureCounts_quant_%.png: $(DATADIR)/featureCounts_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(DATADIR)/R/test_featureCounts.R --libfile "$(DATADIR)/R/lib.R"  --featureCounts_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(DATADIR)/featureCounts_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

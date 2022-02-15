$(ROOTDIR)/featureCounts_quant_STAR_%.tsv:$(ROOTDIR)/STAR_%.bam
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

$(ROOTDIR)/yasim_to_featureCounts_quant_%.png: $(ROOTDIR)/featureCounts_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_featureCounts.R --libfile "$(ROOTDIR)/R/lib.R"  --featureCounts_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(ROOTDIR)/featureCounts_quant_hisat2_%.tsv:$(ROOTDIR)/hisat2_%.bam
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

$(DATADIR)/featureCounts_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<
$(DATADIR)/featureCounts_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
	featureCounts -t transcript -g transcript_id -p -a "$(REFERENCE_GTF)" -T $(THREADS) -o $@ $<

.PHONY: FEATURECOUNTS_QUANT_ALL
FEATURECOUNTS_QUANT_ALL:$(subst .bam,.tsv,$(subst $(DATADIR)/,$(DATADIR)/featureCounts_quant_,$(wildcard $(DATADIR)/*.bam)))

$(DATADIR)/featureCounts_quant_minimap2_%.tsv:$(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	featureCounts -t transcript -g transcript_id -a "$(REFERENCE_GTF)" -L -o $@ $<

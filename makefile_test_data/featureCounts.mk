FEATURECOUNTS := featureCounts -O -M --primary --ignoreDup -a $(REFERENCE_GTF) -T $(THREADS)

$(DATADIR)/featureCounts_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g transcript_id -p -o $@ $<
$(DATADIR)/featureCounts_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g transcript_id -p -o $@ $<

$(DATADIR)/featureCounts_genelevel_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g gene_id -p -o $@ $<
$(DATADIR)/featureCounts_genelevel_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g gene_id -p -o $@ $<

.PHONY: FEATURECOUNTS_QUANT_ALL
FEATURECOUNTS_QUANT_ALL:$(subst .bam,.tsv,$(subst $(DATADIR)/,$(DATADIR)/featureCounts_quant_,$(wildcard $(DATADIR)/*.bam)))

$(DATADIR)/featureCounts_quant_minimap2_%.tsv:$(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g transcript_id -L -o $@ $<

$(DATADIR)/featureCounts_genelevel_quant_minimap2_%.tsv:$(DATADIR)/minimap2_%.bam $(REFERENCE_GTF)
	$(FEATURECOUNTS) -g gene_id -L -o $@ $<

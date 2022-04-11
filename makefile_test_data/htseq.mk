$(DATADIR)/htseq_quant_%.tsv:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai $(REFERENCE_GTF)
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"
$(DATADIR)/htseq_genelevel_quant_%.tsv:$(DATADIR)/%.bam $(DATADIR)/%.bam.bai $(REFERENCE_GTF)
	htseq-count --idattr gene_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"

$(DATADIR)/htseq_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(DATADIR)/STAR_%.bam.bai $(REFERENCE_GTF)
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"

$(DATADIR)/htseq_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(DATADIR)/hisat2_%.bam.bai $(REFERENCE_GTF)
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"

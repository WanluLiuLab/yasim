$(DATADIR)/htseq_quant_STAR_%.tsv:$(DATADIR)/STAR_%.bam $(DATADIR)/STAR_%.bam.bai
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"

$(DATADIR)/yasim_to_htseq_quant_%.png: $(DATADIR)/htseq_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(DATADIR)/R/test_htseq.R --libfile "$(DATADIR)/R/lib.R"  --htseq_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(DATADIR)/htseq_quant_hisat2_%.tsv:$(DATADIR)/hisat2_%.bam $(DATADIR)/hisat2_%.bam.bai
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"


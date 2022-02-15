$(ROOTDIR)/htseq_quant_STAR_%.tsv:$(ROOTDIR)/STAR_%.bam $(ROOTDIR)/STAR_%.bam.bai
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"

$(ROOTDIR)/yasim_to_htseq_quant_%.png: $(ROOTDIR)/htseq_quant_%.tsv $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_htseq.R --libfile "$(ROOTDIR)/R/lib.R"  --htseq_quant_tsv $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(ROOTDIR)/htseq_quant_hisat2_%.tsv:$(ROOTDIR)/hisat2_%.bam $(ROOTDIR)/hisat2_%.bam.bai
	htseq-count --idattr transcript_id -n $(THREADS) -c $@ $< "$(REFERENCE_GTF)"


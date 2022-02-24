$(DATADIR)/minimap2_%.bam:$(FASTQ_BASENAME)_%.fq $(REFERENCE_FASTA)
	minimap2 -x splice -a -t $(THREADS)  $(REFERENCE_FASTA) $< | samtools sort - -@ 40 -o $@

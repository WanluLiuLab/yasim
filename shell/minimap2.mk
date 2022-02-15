$(ROOTDIR)/minimap2_%.bam:$(FASTQ_BASENAME)_%.fq
	minimap2 -x splice -a -t $(THREADS)  $(REFERENCE_FASTA) $< | samtools sort - -@ 40 -o $@

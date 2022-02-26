export MINIMAP2_INDEX := $(REFERENCE_FASTA).mmi

.PHONY: MINIMAP2_INDEX
MINIMAP2_INDEX: $(MINIMAP2_INDEX)

$(MINIMAP2_INDEX): $(REFERENCE_FASTA)
	minimap2 -d $@ $<

$(DATADIR)/minimap2_%.bam:$(FASTQ_BASENAME)_%.fq $(MINIMAP2_INDEX)
	minimap2 -x splice -a -t $(THREADS)  $(MINIMAP2_INDEX) $< | samtools sort - -@ 40 -o $@

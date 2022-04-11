export SALMON_INDEX := $(DATADIR)/ce11_salmon_index

.PHONY: SALMON_INDEX
SALMON_INDEX: $(SALMON_INDEX)

$(SALMON_INDEX): $(REFERENCE_CDNA)
	salmon index --transcripts "$(REFERENCE_CDNA)" --index $@ -p $(THREADS)

$(DATADIR)/salmon_quant_%/quant.sf:$(FASTQ_BASENAME)_%_1.fq $(SALMON_INDEX)
	salmon quant -i "$(SALMON_INDEX)" --threads $(THREADS) -l IU -1 "$<" -2 $(subst _1.fq,_2.fq,$<) -o $(dir $@)

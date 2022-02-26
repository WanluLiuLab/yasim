export SALMON_INDEX := $(DATADIR)/ce11_salmon_index

.PHONY: SALMON_INDEX
SALMON_INDEX: $(SALMON_INDEX)

$(SALMON_INDEX): $(REFERENCE_CDNA)
	salmon index --transcripts "$(REFERENCE_CDNA)" --index $@ -p $(THREADS)

$(DATADIR)/salmon_quant_%/quant.sf:$(FASTQ_BASENAME)_%_1.fq $(SALMON_INDEX)
	salmon quant -i "$(SALMON_INDEX)" --threads $(THREADS) -l IU -1 "$<" -2 $(subst _1.fq,_2.fq,$<) -o $(dir $@)

$(DATADIR)/yasim_to_salmon_quant_%.png: $(DATADIR)/salmon_quant_%/quant.sf $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_salmon.R --libfile "$(ROOTDIR)/R/lib.R"  --salmon_quant_sf $< --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

export SALMON_INDEX := /home/yuzj/Desktop/BioRef/chr1_salmon_index


$(SALMON_INDEX): $(REFERENCE_CDNA)
	salmon index --transcripts "$(REFERENCE_CDNA)" --index $@ -p $(THREADS)

$(ROOTDIR)/salmon_quant_%:$(FASTQ_BASENAME)_%_1.fq $(SALMON_INDEX)
	salmon quant -i "$(SALMON_INDEX)" --threads $(THREADS) -l IU -1 "$<" -2 $(subst _1.fq,_2.fq,$<) -o $@

$(ROOTDIR)/yasim_to_salmon_quant_%.png: $(ROOTDIR)/salmon_quant_% $(DEPTH_TSV)
	 Rscript $(ROOTDIR)/R/test_salmon.R --libfile "$(ROOTDIR)/R/lib.R"  --salmon_quant_sf $</quant.sf --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

CPPTETGS := $(HOME)/Documents/cpptetgs_experimental/cmake-build-release-llvm/cpptetgs

$(DATADIR)/cpptetgs_quant_%.tsv:$(DATADIR)/%.bam $(REFERENCE_SORTED_GTF)
	$(CPPTETGS) -G "$(REFERENCE_SORTED_GTF)" --nt-io 10 --nt-wk 2 -o $@ -B $<

export HISAT2_INDEX_BASENAME := /home/yuzj/Desktop/BioRef/chr1_HISAT2_index/chr1
export HISAT2_INDEX_DIRNAME := $(dir $(HISAT2_INDEX_BASENAME))
export HISAT2_INDEX_EXON := $(HISAT2_INDEX_DIRNAME)exon
export HISAT2_INDEX_SS := $(HISAT2_INDEX_DIRNAME)ss
export HISAT2_INDEX_FILENAME := /home/yuzj/Desktop/BioRef/chr1_HISAT2_index/chr1.1.ht2

.PHONY: HISAT2_INDEX
HISAT2_INDEX:$(HISAT2_INDEX_FILENAME)

#$(HISAT2_INDEX_EXON): $(REFERENCE_GTF)
#	mkdir -p $(HISAT2_INDEX_DIRNAME)
#	hisat2_extract_exons.py $< > $@
#
#$(HISAT2_INDEX_SS): $(REFERENCE_GTF) $(HISAT2_INDEX_DIRNAME)
#	mkdir -p $(HISAT2_INDEX_DIRNAME)
#	hisat2_extract_splice_sites.py $< > $@
#
#$(HISAT2_INDEX_FILENAME): $(REFERENCE_FASTA) $(HISAT2_INDEX_EXON) $(HISAT2_INDEX_SS) # FIXME: This step executes all the time!
#	hisat2-build -p $(THREADS) --ss $(HISAT2_INDEX_SS) --exon $(HISAT2_INDEX_EXON) $< $(HISAT2_INDEX_BASENAME)

$(DATADIR)/hisat2_%.bam:$(FASTQ_BASENAME)_%_1.fq $(HISAT2_INDEX_FILENAME)
	hisat2 -p $(THREADS) -x $(HISAT2_INDEX_BASENAME) -1 "$<" -2 $(subst _1.fq,_2.fq,$<) |\
	samtools sort -@ $(THREADS) -o $@

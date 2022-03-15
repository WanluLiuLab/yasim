export STAR_INDEX := $(DATADIR)/ce11_STAR_index
.PHONY: STAR_INDEX

STAR_INDEX:$(STAR_INDEX)


$(STAR_INDEX): $(REFERENCE_FASTA) $(REFERENCE_GTF)
	STAR --runThreadN $(THREADS) \
	--runMode genomeGenerate \
	--genomeDir "$@" \
	--genomeFastaFiles "$(REFERENCE_FASTA)" \
	--sjdbGTFfile "$(REFERENCE_GTF)"

$(DATADIR)/STAR_%.bam:$(FASTQ_BASENAME)_%_1.fq $(STAR_INDEX)
	STAR --genomeDir "$(STAR_INDEX)" \
	--runThreadN $(THREADS) \
	--readFilesIn "$<" $(subst _1.fq,_2.fq,$<)  \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix $(basename $@)/ \
	--quantMode GeneCounts
	ln $(basename $@)/Aligned.sortedByCoord.out.bam $@

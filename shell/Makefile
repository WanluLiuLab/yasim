
export ROOTDIR := $(abspath $(CURDIR)/../)
export PYTHONPATH := "$(ROOTDIR)/src:$(PYTHONPATH:-)"
export YASIM := PYTHONPATH="$(ROOTDIR)/src:$(PYTHONPATH:-)" python -m yasim
export REFERENCE_FASTA := /home/yuzj/Desktop/BioRef/chr1.fa
export REFERENCE_GTF := /home/yuzj/Desktop/BioRef/chr1.ncbiRefSeq.gtf
export REFERENCE_CDNA := $(ROOTDIR)/hg38.chr1.reference_transcripts.fa
export SALMON_INDEX := /home/yuzj/Desktop/BioRef/chr1_salmon_index
export THREADS := 40
export SEL_GTF := $(ROOTDIR)/Homo_sapiens.GRCh38.105_chr1_sel.gtf
export SEL_CDNA_FASTA := $(ROOTDIR)/hg38.chr1.transcripts_sel.fa
export DEPTH_TSV := $(ROOTDIR)/hg38.chr1.gene.depth.tsv
export FASTQ_BASENAME := $(ROOTDIR)/hg38.chr1.transcripts
export STAR_INDEX := /home/yuzj/Desktop/BioRef/chr1_STAR_index
# export R_ENVIRON_USER := $(ROOTDIR)/R/.Rprofile

$(SEL_GTF):
	$(YASIM) alternative_splicing -f "${REFERENCE_FASTA}" -g "${REFERENCE_GTF}" -o $@

$(SEL_CDNA_FASTA): $(SEL_GTF)
	$(YASIM) transcribe -f "${REFERENCE_FASTA}" -g "${SEL_GTF}" -o $@

$(DEPTH_TSV): $(SEL_GTF)
	$(YASIM) dge -g "$(SEL_GTF)" -o "$(DEPTH_TSV)" -d 100 -l 20

$(REFERENCE_CDNA):
	$(YASIM) transcribe -g "$(REFERENCE_GTF)" -o $@ -f "$(REFERENCE_FASTA)"

${SALMON_INDEX}: $(REFERENCE_CDNA)
	salmon index --transcripts "$(REFERENCE_CDNA)" --index $@ -p "$(THREADS)"

$(STAR_INDEX):
	STAR --runThreadN "$(THREADS)" \
	--runMode genomeGenerate \
	--genomeDir "$@" \
	--genomeFastaFiles "$(REFERENCE_FASTA)" \
	--sjdbGTFfile "$(REFERENCE_GTF)"

$(FASTQ_BASENAME)_dwgsim_1.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) dwgsim -F "$(SEL_CDNA_FASTA)" -o "$(FASTQ_BASENAME)"_dwgsim -d "$(DEPTH_TSV)"

$(ROOTDIR)/salmon_quant_%:$(FASTQ_BASENAME)_%_1.fq $(SALMON_INDEX)
	salmon quant -i "$(SALMON_INDEX)" --threads "${THREADS}" -l IU -1 "$<" -2 $(subst _1.fq,_2.fq,$<) -o $@

$(ROOTDIR)/yasim_to_salmon_quant_%.png: $(ROOTDIR)/salmon_quant_%
	 Rscript $(ROOTDIR)/R/test_salmon.R --libfile "$(ROOTDIR)/R/lib.R"  --salmon_quant_sf $</quant.sf --yasim_tsv "$(DEPTH_TSV)" --output $(basename $@)

$(ROOTDIR)/STAR_%.bam:$(FASTQ_BASENAME)_%_1.fq $(STAR_INDEX)
	STAR --genomeDir "$(STAR_INDEX)" \
	--runThreadN "${THREADS}" \
	--readFilesIn "$<" $(subst _1.fq,_2.fq,$<)  \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix $(basename $@)/
	ln -s $(basename $@)/Aligned.sortedByCoord.out.bam $@

$(ROOTDIR)/featureCounts_quant_%.tsv:$(ROOTDIR)/STAR_%.bam $(SALMON_INDEX)
	featureCounts -t transcript -p -a "$(REFERENCE_GTF)" -T "${THREADS}" -o $@ $<

# --libfile "$(ROOTDIR)/R/lib.R"

all: $(ROOTDIR)/yasim_to_salmon_quant_dwgsim.png $(ROOTDIR)/featureCounts_quant_dwgsim.tsv



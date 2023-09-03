.PHONY: all
all: deps llrg

.PHONY: llrg
llrg: chrM.STAR.ref.d chrM.trans.fa chrI.trans.fa chrI.STAR.ref.d
	$(MAKE) -C llrg

.PHONY: deps
deps: ce11.ncbiRefSeq.chrM.gtf chrM.fa ce11.ncbiRefSeq.chrM.as.gtf ce11.ncbiRefSeq.chrI.as.gtf

.SUFFIXES:
.SECONDARY: chromFa.tar.gz \
	ce11.ncbiRefSeq.gtf \
	ce11.ncbiRefSeq.gtf.gz

chromFa.tar.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz

chrM.fa: chromFa.tar.gz
	tar xzvf chromFa.tar.gz
	touch chrI.fa chrII.fa chrIII.fa chrIV.fa chrM.fa
	gatk CreateSequenceDictionary -R chrI.fa
	gatk CreateSequenceDictionary -R chrM.fa

ce11.ncbiRefSeq.gtf.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
	touch ce11.ncbiRefSeq.gtf.gz

ce11.ncbiRefSeq.gtf: ce11.ncbiRefSeq.gtf.gz
	gunzip -f ce11.ncbiRefSeq.gtf.gz
	touch ce11.ncbiRefSeq.gtf

ce11.ncbiRefSeq.chrM.gtf: ce11.ncbiRefSeq.gtf
	grep '^chrM\s' ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chrM.gtf

ce11.ncbiRefSeq.chrI.gtf: ce11.ncbiRefSeq.gtf
	grep '^chrI\s' ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chrI.gtf

ce11.ncbiRefSeq.chrM.as.gtf:ce11.ncbiRefSeq.chrM.gtf chrM.fa
	python -m yasim generate_as_events \
		-g ce11.ncbiRefSeq.chrM.gtf \
		-f chrM.fa \
		-c 5 \
		-o ce11.ncbiRefSeq.chrM.as.gtf

ce11.ncbiRefSeq.chrI.as.gtf:ce11.ncbiRefSeq.chrI.gtf chrI.fa
	python -m yasim generate_as_events \
		-g ce11.ncbiRefSeq.chrI.gtf \
		-f chrI.fa \
		-c 5 \
		-o ce11.ncbiRefSeq.chrI.as.gtf

.PHONY: clean
clean:
	$(MAKE) -C llrg clean
	rm -f *.gvpkl.xz


chrM.STAR.ref.d: chrM.fa ce11.ncbiRefSeq.chrM.gtf
	STAR \
		--runThreadN 40 \
		--runMode genomeGenerate \
		--genomeDir chrM.STAR.ref.d \
		--genomeFastaFiles chrM.fa \
		--sjdbGTFfile ce11.ncbiRefSeq.chrM.gtf

chrI.STAR.ref.d: chrI.fa ce11.ncbiRefSeq.chrI.gtf
	STAR \
		--runThreadN 40 \
		--runMode genomeGenerate \
		--genomeDir chrI.STAR.ref.d \
		--genomeFastaFiles chrI.fa \
		--sjdbGTFfile ce11.ncbiRefSeq.chrI.gtf

chrM.trans.fa: chrM.fa ce11.ncbiRefSeq.chrM.gtf
	python -m labw_utils.bioutils transcribe \
		-g ce11.ncbiRefSeq.chrM.gtf \
		-f chrM.fa \
		-o chrM.trans.fa

chrI.trans.fa: chrI.fa ce11.ncbiRefSeq.chrI.gtf
	python -m labw_utils.bioutils transcribe \
		-g ce11.ncbiRefSeq.chrI.gtf \
		-f chrI.fa \
		-o chrI.trans.fa

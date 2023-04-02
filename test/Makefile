.PHONY: all
all: deps llrg

.PHONY: llrg
llrg:
	$(MAKE) -C llrg

.PHONY: deps
deps: ce11.ncbiRefSeq.chrM.gtf chrM.fa

.SUFFIXES:
.SECONDARY: chromFa.tar.gz \
	chrI.fa \
	chrII.fa \
	chrIII.fa \
	chrIV.fa \
	chrM.fa \
	ce11.ncbiRefSeq.chrI.gtf \
	ce11.ncbiRefSeq.chrII.gtf \
	ce11.ncbiRefSeq.chrIII.gtf \
	ce11.ncbiRefSeq.chrIV.gtf \
	ce11.ncbiRefSeq.chrM.gtf \
	ce11.ncbiRefSeq.chrV.gtf \
	ce11.ncbiRefSeq.chrX.gtf \
	ce11.ncbiRefSeq.gtf

chromFa.tar.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz

chrM.fa: chromFa.tar.gz
	tar xzvf chromFa.tar.gz
	touch chrI.fa chrII.fa chrIII.fa chrIV.fa chrM.fa

ce11.ncbiRefSeq.gtf.gz:
	axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz

ce11.ncbiRefSeq.gtf: ce11.ncbiRefSeq.gtf.gz
	gunzip ce11.ncbiRefSeq.gtf.gz

ce11.ncbiRefSeq.chrM.gtf: ce11.ncbiRefSeq.gtf
	for chrname in I II III IV V X M; do \
		grep '^chr'"${chrname}"'\s' ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr"${chrname}".gtf; \
		done

.PHONY: clean
clean:
	$(MAKE) -C llrg clean

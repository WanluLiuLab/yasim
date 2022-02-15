
all: $(CURDIR)/pacbio_mkindex.pl $(CURDIR)/PaSS_exe

$(CURDIR)/pacbio_mkindex.pl: $(CURDIR)/PaSS
	cp $</pacbio_mkindex.pl ./

$(CURDIR)/PaSS_exe:$(CURDIR)/PaSS
	$(CC) -lm -lpthread $</PaSS.c -o $@

$(CURDIR)/PaSS.tar.gz:
	axel https://cgm.sjtu.edu.cn/PaSS/src/PaSS.tar.gz

$(CURDIR)/PaSS:$(CURDIR)/PaSS.tar.gz
	gzip -dcf $< | tar xvf -
	rm -rf PaSS/C.elegan/ PaSS/E.coli/ PaSS/Arabidopsis/

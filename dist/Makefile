
all: $(CURDIR)/pacbio_mkindex.pl $(CURDIR)/PaSS_exe

$(CURDIR)/pacbio_mkindex.pl: $(CURDIR)/PaSS
	cp $</pacbio_mkindex.pl ./

$(CURDIR)/PaSS_exe:$(CURDIR)/PaSS
	gcc -lm -lpthread $</PaSS.c -o $@

$(CURDIR)/PaSS.tar.gz:
	axel https://cgm.sjtu.edu.cn/PaSS/src/PaSS.tar.gz

$(CURDIR)/PaSS:PaSS.tar.gz
	gzip -dcf $< | tar xvf -
	rm -rf PaSS/C.elegan/ PaSS/E.coli/ PaSS/Arabidopsis/

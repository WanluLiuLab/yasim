# Makefile for reference preparation
all: TE_gene_transcript_class.tsv TE_family_count.tsv # hg38.fai hg38.ncbiRefSeq.gtf

#rmsk.txt.gz:
#	axel https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

%: %.gz
	gzip -dfk $<

GRCh38_GENCODE_rmsk_TE.gtf.gz:
	axel http://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz

TE_gene_transcript_class.tsv: GRCh38_GENCODE_rmsk_TE.gtf
	printf "gene_id\ttranscript_id\tclass_id\n" > $@
	cat $< | cut -f 9 | cut -f 2,4,8 -d \" | tr '"' '\t' >> $@

TE_family_count.tsv: GRCh38_GENCODE_rmsk_TE.gtf
	printf "gene_id\ttotal_count\n" > $@
	cat $< | cut -f 9 | cut -f 2 -d \" | awk -f count_te_families.awk >> $@
#
#hg38.fa.gz:
#	axel https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#
#hg38.fai: hg38.fa
#	samtools faidx $<
#
#hg38.ncbiRefSeq.gtf.gz:
#	axel https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
#
#hg38_STAR_ref: hg38.fa hg38.ncbiRefSeq.gtf
#	STAR --runThreadN 36 --runMode genomeGenerate --genomeDir $PWD/hg38_STAR_ref --genomeFastaFiles hg38.fa --sjdbGTFfile hg38.ncbiRefSeq.gtf

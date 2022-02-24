
%.bam.bai: %.bam
	samtools index $<

%:%.gz
	gzip -dkf $<

%.fa.fai:%.fa
	samtools faidx $<

%.fasta.fai:%.fasta
	samtools faidx $<

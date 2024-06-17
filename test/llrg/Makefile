export PYTHONPATH := ../../src:../../deps/labw_utils/src

.PHONY: all
all: pbsim3_genome_mode \
	pbsim3_trans_mode \
	pbsim3_ccs_mode \
	pbsim3_low_accuracy_mode \
	pbsim3_errhmm_mode \
	art_mode \
	art_pe_mode \
	dtgs_mode

.PHONY: distclean
distclean: clean

.PHONY: clean
clean:
	rm -rf *mode* chrM.depth.tsv chrM.depth_gene.tsv *.log

chrM.depth_gene.tsv: ../ce11.ncbiRefSeq.chrM.gtf
	python -m yasim generate_gene_depth \
		-g ../ce11.ncbiRefSeq.chrM.gtf \
		-o chrM.depth_gene.tsv \
		--mu 50

chrM.depth.tsv: ../ce11.ncbiRefSeq.chrM.gtf chrM.depth_gene.tsv
	python -m yasim generate_isoform_depth \
		-g ../ce11.ncbiRefSeq.chrM.gtf \
		-o chrM.depth.tsv \
		-d chrM.depth_gene.tsv

pbsim3_genome_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-e pbsim \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_genome_mode \
		-j 40 \
		--strategy wgs
	minimap2 -a -x splice ../chrM.fa pbsim3_genome_mode.fq | \
		samtools sort -o pbsim3_genome_mode.bam
	samtools index pbsim3_genome_mode.bam
	samtools flagstat pbsim3_genome_mode.bam > pbsim3_genome_mode.bam.stat
	touch pbsim3_genome_mode

pbsim3_genome_ccs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-e pbsim \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_genome_ccs_mode \
		-j 40 \
		--ccs_pass 10 \
		--strategy wgs
	minimap2 -a -x splice ../chrM.fa pbsim3_genome_ccs_mode.fq | \
		samtools sort -o pbsim3_genome_ccs_mode.bam
	samtools index pbsim3_genome_ccs_mode.bam
	samtools flagstat pbsim3_genome_ccs_mode.bam > pbsim3_genome_ccs_mode.bam.stat
	touch pbsim3_genome_ccs_mode

pbsim3_trans_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-e pbsim \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_trans_mode \
		-j 40 \
		--strategy trans
	minimap2 -a -x splice ../chrM.fa pbsim3_trans_mode.fq | \
		samtools sort -o pbsim3_trans_mode.bam
	samtools index pbsim3_trans_mode.bam
	samtools flagstat pbsim3_trans_mode.bam > pbsim3_trans_mode.bam.stat
	touch pbsim3_trans_mode

pbsim3_ccs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-e pbsim \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_ccs_mode \
		-j 40 \
		--ccs_pass 10
	minimap2 -a -x splice ../chrM.fa pbsim3_ccs_mode.fq | \
		samtools sort -o pbsim3_ccs_mode.bam
	samtools index pbsim3_ccs_mode.bam
	samtools flagstat pbsim3_ccs_mode.bam > pbsim3_ccs_mode.bam.stat
	touch pbsim3_ccs_mode

pbsim3_low_accuracy_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-e pbsim \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_low_accuracy_mode \
		-j 40 \
		--accuracy-mean 0.4
	minimap2 -a -x splice ../chrM.fa pbsim3_low_accuracy_mode.fq | \
		samtools sort -o pbsim3_low_accuracy_mode.bam
	samtools index pbsim3_low_accuracy_mode.bam
	samtools flagstat pbsim3_low_accuracy_mode.bam > pbsim3_low_accuracy_mode.bam.stat
	touch pbsim3_low_accuracy_mode

pbsim3_errhmm_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m SEQUEL \
		-e pbsim \
		-M errhmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_errhmm_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim3_errhmm_mode.fq | \
		samtools sort -o pbsim3_errhmm_mode.bam
	samtools index pbsim3_errhmm_mode.bam
	samtools flagstat pbsim3_errhmm_mode.bam > pbsim3_errhmm_mode.bam.stat
	touch pbsim3_errhmm_mode

art_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim art \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o art_mode \
		-j 40
	touch art_mode

art_pe_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim art \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o art_pe_mode \
		--is_pair_end \
		--pair_end_fragment_length_mean 200 \
		--pair_end_fragment_length_std 5 \
		-j 40
	STAR \
		--twopassMode Basic \
		--quantMode TranscriptomeSAM GeneCounts \
		--runThreadN 40 \
		--genomeDir ../chrM.STAR.ref.d \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix art_pe_mode_2_STAR/ \
		--readFilesIn art_pe_mode_1.fq art_pe_mode_2.fq
	samtools sort art_pe_mode_2_STAR/Aligned.out.bam -o art_pe_mode_2_STAR/Aligned.sortedByCoord.out.bam
	samtools index art_pe_mode_2_STAR/Aligned.sortedByCoord.out.bam
	samtools flagstat art_pe_mode_2_STAR/Aligned.sortedByCoord.out.bam > art_pe_mode_2_STAR.stat
	touch art_pe_mode

dtgs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim dtgs \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o dtgs_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa dtgs_mode.fq | \
		samtools sort -o dtgs_mode.bam
	samtools index dtgs_mode.bam
	samtools flagstat dtgs_mode.bam > dtgs_mode.bam.stat
	touch dtgs_mode

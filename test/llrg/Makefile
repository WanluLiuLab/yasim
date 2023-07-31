export PYTHONPATH := ../../src:../../deps/labw_utils/src

.PHONY: all
all: pbsim3_genome_mode \
	pbsim3_trans_mode \
	pbsim3_ccs_mode \
	pbsim3_low_accuracy_mode \
	pbsim3_errhmm_mode \
	dwgsim_mode \
	art_mode \
	art_pe_mode \
	pbsim_ccs_mode \
	pbsim_clr_mode \
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

chrI.depth_gene.tsv: ../ce11.ncbiRefSeq.chrI.gtf
	python -m yasim generate_gene_depth \
		-g ../ce11.ncbiRefSeq.chrI.gtf \
		-o chrI.depth_gene.tsv \
		--mu 5

chrI.depth.tsv: ../ce11.ncbiRefSeq.chrI.gtf chrI.depth_gene.tsv
	python -m yasim generate_isoform_depth \
		-g ../ce11.ncbiRefSeq.chrI.gtf \
		-o chrI.depth.tsv \
		-d chrI.depth_gene.tsv

pbsim3_genome_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_genome_mode \
		-j 40 \
		--strategy wgs
	minimap2 -a -x splice ../chrM.fa pbsim3_genome_mode.fq | \
		samtools sort -o pbsim3_genome_mode.bam
	samtools index pbsim3_genome_mode.bam
	touch pbsim3_genome_mode

pbsim3_genome_ccs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
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
	touch pbsim3_genome_ccs_mode

pbsim3_trans_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_trans_mode \
		-j 40 \
		--strategy trans
	minimap2 -a -x splice ../chrM.fa pbsim3_trans_mode.fq | \
		samtools sort -o pbsim3_trans_mode.bam
	samtools index pbsim3_trans_mode.bam
	touch pbsim3_trans_mode

pbsim3_ccs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_ccs_mode \
		-j 40 \
		--ccs_pass 10
	minimap2 -a -x splice ../chrM.fa pbsim3_ccs_mode.fq | \
		samtools sort -o pbsim3_ccs_mode.bam
	samtools index pbsim3_ccs_mode.bam
	touch pbsim3_ccs_mode

pbsim3_low_accuracy_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_ccs_mode \
		-j 40 \
		--accuracy-mean 0.4
	minimap2 -a -x splice ../chrM.fa pbsim3_low_accuracy_mode.fq | \
		samtools sort -o pbsim3_low_accuracy_mode.bam
	samtools index pbsim3_low_accuracy_mode.bam
	touch pbsim3_low_accuracy_mode

pbsim3_errhmm_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim3 \
		-m SEQUEL \
		-M errhmm \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim3_errhmm_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim3_errhmm_mode.fq | \
		samtools sort -o pbsim3_errhmm_mode.bam
	samtools index pbsim3_errhmm_mode.bam
	touch pbsim3_errhmm_mode

dwgsim_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim dwgsim \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o dwgsim_mode \
		-j 40
	touch dwgsim_mode

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
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix art_pe_mode_2_STAR/ \
		--readFilesIn art_pe_mode_1.fq art_pe_mode_2.fq
	samtools index art_pe_mode_2_STAR/Aligned.sortedByCoord.out.bam
	touch art_pe_mode

pbsim2_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim2 \
		-m P4C2 \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim2_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim2_mode.fq | \
		samtools sort -o pbsim2_mode.bam
	samtools index pbsim2_mode.bam
	touch pbsim2_mode

pbsim_clr_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim_clr_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim_clr_mode.fq | \
		samtools sort -o pbsim_clr_mode.bam
	samtools index pbsim_clr_mode.bam
	touch pbsim_clr_mode

pbsim_ccs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim_ccs_mode \
		--ccs \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim_ccs_mode.fq | \
		samtools sort -o pbsim_ccs_mode.bam
	samtools index pbsim_ccs_mode.bam
	touch pbsim_ccs_mode

pbsim_ccs_mode_pi: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim pbsim \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o pbsim_ccs_mode \
		--ccs \
		-j 40 \
		--preserve_intermediate_files
	minimap2 -a -x splice ../chrM.fa pbsim_ccs_mode.fq | \
		samtools sort -o pbsim_ccs_mode.bam
	samtools index pbsim_ccs_mode.bam
	touch pbsim_ccs_mode_pi

dtgs_mode: ../chrM.trans.fa chrM.depth.tsv
	python -m yasim dtgs \
		-F ../chrM.trans.fa.d \
		-d chrM.depth.tsv \
		-o dtgs_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa dtgs_mode.fq | \
		samtools sort -o dtgs_mode.bam
	samtools index dtgs_mode.bam
	touch dtgs_mode

chrI.dtgs_mode: ../chrI.trans.fa chrI.depth.tsv
	python -m yasim dtgs \
		-F ../chrI.trans.fa.d \
		-d chrI.depth.tsv \
		-o chrI.dtgs_mode \
		-j 40
	minimap2 -t 40 -a -x splice ../chrI.fa chrI.dtgs_mode.fq | \
		samtools sort -o chrI.dtgs_mode.bam
	samtools index chrI.dtgs_mode.bam
	touch chrI.dtgs_mode

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
	rm -rf chrM.trans.fa chrM.trans.fa.d chrM.trans.fa.stats depth.tsv depth_gene.tsv

.PHONY: clean
clean:
	rm -rf *mode*


chrM.trans.fa: ../chrM.fa ../ce11.ncbiRefSeq.chrM.gtf
	python -m yasim transcribe \
		-g ../ce11.ncbiRefSeq.chrM.gtf \
		-f ../chrM.fa \
		-o chrM.trans.fa

depth_gene.tsv: ../ce11.ncbiRefSeq.chrM.gtf
	python -m yasim generate_gene_depth \
		-g ../ce11.ncbiRefSeq.chrM.gtf \
		-o depth_gene.tsv \
		--mu 50

depth.tsv: ../ce11.ncbiRefSeq.chrM.gtf depth_gene.tsv
	python -m yasim generate_isoform_depth \
		-g ../ce11.ncbiRefSeq.chrM.gtf \
		-o depth.tsv \
		-d depth_gene.tsv

pbsim3_genome_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_genome_mode \
		-j 40 \
		--strategy wgs
	minimap2 -a -x splice ../chrM.fa pbsim3_genome_mode.fq | \
		samtools sort -o pbsim3_genome_mode.bam
	samtools index pbsim3_genome_mode.bam
	touch pbsim3_genome_mode

pbsim3_genome_ccs_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_genome_ccs_mode \
		-j 40 \
		--ccs_pass 10 \
		--strategy wgs
	minimap2 -a -x splice ../chrM.fa pbsim3_genome_ccs_mode.fq | \
		samtools sort -o pbsim3_genome_ccs_mode.bam
	samtools index pbsim3_genome_ccs_mode.bam
	touch pbsim3_genome_ccs_mode

pbsim3_trans_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_trans_mode \
		-j 40 \
		--strategy trans
	minimap2 -a -x splice ../chrM.fa pbsim3_trans_mode.fq | \
		samtools sort -o pbsim3_trans_mode.bam
	samtools index pbsim3_trans_mode.bam
	touch pbsim3_trans_mode

pbsim3_ccs_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_ccs_mode \
		-j 40 \
		--ccs_pass 10
	minimap2 -a -x splice ../chrM.fa pbsim3_ccs_mode.fq | \
		samtools sort -o pbsim3_ccs_mode.bam
	samtools index pbsim3_ccs_mode.bam
	touch pbsim3_ccs_mode

pbsim3_low_accuracy_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m RSII \
		-M qshmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_ccs_mode \
		-j 40 \
		--accuracy-mean 0.4
	minimap2 -a -x splice ../chrM.fa pbsim3_low_accuracy_mode.fq | \
		samtools sort -o pbsim3_low_accuracy_mode.bam
	samtools index pbsim3_low_accuracy_mode.bam
	touch pbsim3_low_accuracy_mode

pbsim3_errhmm_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim3 \
		-m SEQUEL \
		-M errhmm \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim3_errhmm_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim3_errhmm_mode.fq | \
		samtools sort -o pbsim3_errhmm_mode.bam
	samtools index pbsim3_errhmm_mode.bam
	touch pbsim3_errhmm_mode

dwgsim_mode: chrM.trans.fa depth.tsv
	python -m yasim dwgsim \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o dwgsim_mode \
		-j 40
	touch dwgsim_mode

art_mode: chrM.trans.fa depth.tsv
	python -m yasim art \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o art_mode \
		-j 40
	touch art_mode

art_pe_mode: chrM.trans.fa depth.tsv
	python -m yasim art \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o art_pe_mode \
		--is_pair_end \
		--pair_end_fragment_length_mean 200 \
		--pair_end_fragment_length_std 5 \
		-j 40
	touch art_pe_mode

pbsim2_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim2 \
		-m P4C2 \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim2_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim2_mode.fq | \
		samtools sort -o pbsim2_mode.bam
	samtools index pbsim2_mode.bam
	touch pbsim2_mode


pbsim_clr_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim_clr_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim_clr_mode.fq | \
		samtools sort -o pbsim_clr_mode.bam
	samtools index pbsim_clr_mode.bam
	touch pbsim_clr_mode

pbsim_ccs_mode: chrM.trans.fa depth.tsv
	python -m yasim pbsim \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o pbsim_ccs_mode \
		--ccs \
		-j 40
	minimap2 -a -x splice ../chrM.fa pbsim_ccs_mode.fq | \
		samtools sort -o pbsim_ccs_mode.bam
	samtools index pbsim_ccs_mode.bam
	touch pbsim_ccs_mode

dtgs_mode: chrM.trans.fa depth.tsv
	python -m yasim dtgs \
		-F chrM.trans.fa.d \
		-d depth.tsv \
		-o dtgs_mode \
		-j 40
	minimap2 -a -x splice ../chrM.fa dtgs_mode.fq | \
		samtools sort -o dtgs_mode.bam
	samtools index dtgs_mode.bam
	touch dtgs_mode

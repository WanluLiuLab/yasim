$(FASTQ_BASENAME)_dwgsim_1.fq: $(SEL_CDNA_FASTA)
	$(YASIM) dwgsim -F "$(SEL_CDNA_FASTA)".d -o "$(FASTQ_BASENAME)"_dwgsim

$(FASTQ_BASENAME)_pbsim_clr.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$(SEL_CDNA_FASTA)".d -o "$(FASTQ_BASENAME)"_pbsim_clr
	rm -rf $(FASTQ_BASENAME)_pbsim_clr.d

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$(SEL_CDNA_FASTA)".d -o "$(FASTQ_BASENAME)"_pbsim_ccs --ccs
	rm -rf $(FASTQ_BASENAME)_pbsim_ccs.d

$(FASTQ_BASENAME)_pbsim2_%.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim2 -F "$(SEL_CDNA_FASTA)".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@)))
	rm -rf $(basename $@).d

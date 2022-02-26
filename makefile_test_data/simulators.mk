$(FASTQ_BASENAME)_dwgsim_1.fq: $(SEL_CDNA_FASTA_D)
	$(YASIM) dwgsim -F "$<" -o "$(FASTQ_BASENAME)"_dwgsim

$(FASTQ_BASENAME)_pbsim_clr.fq: $(SEL_CDNA_FASTA_D)
	$(YASIM) pbsim -F "$<" -o "$(FASTQ_BASENAME)"_pbsim_clr
	# rm -rf $(FASTQ_BASENAME)_pbsim_clr.d

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(SEL_CDNA_FASTA_D)
	$(YASIM) pbsim -F "$<" -o "$(FASTQ_BASENAME)"_pbsim_ccs --ccs
	# rm -rf $(FASTQ_BASENAME)_pbsim_ccs.d

$(FASTQ_BASENAME)_pbsim2_%.fq: $(SEL_CDNA_FASTA_D)
	$(YASIM) pbsim2 -F "$<" -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@)))
	# rm -rf $(basename $@).d

$(FASTQ_BASENAME)_badread_%.fq: $(SEL_CDNA_FASTA_D)
	$(YASIM) badread -F "$<" -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@)))
	# rm -rf $(basename $@).d

# $(FASTQ_BASENAME)_simlord.fq: $(SEL_CDNA_FASTA_D)
# 	$(YASIM) simlord -F "$<" -o $(basename $@)
# 	# rm -rf $(basename $@).d

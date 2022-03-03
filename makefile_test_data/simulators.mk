.PHONY:DWGSIM_FQ
DWGSIM_FQ: $(FASTQ_BASENAME)_dwgsim_1.fq

$(FASTQ_BASENAME)_dwgsim_1.fq: $(SEL_CDNA_FASTA)
	$(YASIM) dwgsim -F "$<".d -o "$(FASTQ_BASENAME)"_dwgsim

$(FASTQ_BASENAME)_dwgsim.fq.stats:$(FASTQ_BASENAME)_dwgsim_1.fq

$(FASTQ_BASENAME)_pbsim_clr.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_clr

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_ccs --ccs

$(FASTQ_BASENAME)_pbsim2_%.fq: $(SEL_CDNA_FASTA)
	$(YASIM) pbsim2 -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@)))

$(FASTQ_BASENAME)_badread_%.fq: $(SEL_CDNA_FASTA)
	$(YASIM) badread -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@)))

# $(FASTQ_BASENAME)_simlord.fq: $(SEL_CDNA_FASTA)
# 	$(YASIM) simlord -F "$<".d -o $(basename $@)
# 	# rm -rf $(basename $@).d

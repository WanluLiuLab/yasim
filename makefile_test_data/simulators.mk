.PHONY:DWGSIM_FQ
DWGSIM_FQ: $(FASTQ_BASENAME)_dwgsim_1.fq

$(FASTQ_BASENAME)_dwgsim_1.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) dwgsim -F "$<".d -o "$(FASTQ_BASENAME)"_dwgsim -d $(DEPTH_TSV)

$(FASTQ_BASENAME)_dwgsim.fq.stats:$(FASTQ_BASENAME)_dwgsim_1.fq

$(FASTQ_BASENAME)_pbsim_clr.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_clr -d $(DEPTH_TSV) -j 1200

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_ccs --ccs -d $(DEPTH_TSV) -j 1200

$(FASTQ_BASENAME)_pbsim2_%.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim2 -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@))) -d $(DEPTH_TSV) -j 1200

$(FASTQ_BASENAME)_badread_%.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) badread -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@))) -d $(DEPTH_TSV)

# $(FASTQ_BASENAME)_simlord.fq: $(SEL_CDNA_FASTA)
# 	$(YASIM) simlord -F "$<".d -o $(basename $@)
# 	# rm -rf $(basename $@).d

$(FASTQ_BASENAME)_dwgsim_1.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) dwgsim -F "$(SEL_CDNA_FASTA)" -o "$(FASTQ_BASENAME)"_dwgsim -d "$(DEPTH_TSV)"

$(FASTQ_BASENAME)_pbsim_clr.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$(SEL_CDNA_FASTA)" -o "$(FASTQ_BASENAME)"_pbsim_clr -d "$(DEPTH_TSV)"
	rm -rf $(FASTQ_BASENAME)_pbsim_clr.d

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) pbsim -F "$(SEL_CDNA_FASTA)" -o "$(FASTQ_BASENAME)"_pbsim_ccs -d "$(DEPTH_TSV)" --ccs
	rm -rf $(FASTQ_BASENAME)_pbsim_ccs.d

$(FASTQ_BASENAME)_pbsim2_%.fq: $(DEPTH_TSV) $(SEL_CDNA_FASTA)
	$(YASIM) pbsim2 -F "$(SEL_CDNA_FASTA)" -o $(basename $@) -d "$(DEPTH_TSV)" -m $(lastword $(subst _,\ ,$(basename $@)))
	rm -rf $(basename $@).d

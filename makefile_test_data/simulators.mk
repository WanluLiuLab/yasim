.PHONY: ALL_FQ
ALL_FQ: DWGSIM_FQ PBSIM_FQ PBSIM2_FQ BADREAD_FQ

.PHONY:DWGSIM_FQ
DWGSIM_FQ: $(FASTQ_BASENAME)_dwgsim_1.fq

$(FASTQ_BASENAME)_dwgsim_1.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) dwgsim -F "$<".d -o "$(FASTQ_BASENAME)"_dwgsim -d $(DEPTH_TSV)

$(FASTQ_BASENAME)_dwgsim.fq.stats:$(FASTQ_BASENAME)_dwgsim_1.fq

.PHONY:PBSIM_FQ
PBSIM_FQ: $(FASTQ_BASENAME)_pbsim_clr.fq $(FASTQ_BASENAME)_pbsim_ccs.fq

$(FASTQ_BASENAME)_pbsim_clr.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_clr -d $(DEPTH_TSV) -j 1200

$(FASTQ_BASENAME)_pbsim_ccs.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim -F "$<".d -o "$(FASTQ_BASENAME)"_pbsim_ccs --ccs -d $(DEPTH_TSV) -j 1200


.PHONY:PBSIM2_FQ
PBSIM2_FQ: \
$(FASTQ_BASENAME)_pbsim2_P5C3.fq \
$(FASTQ_BASENAME)_pbsim2_P4C2.fq \
 $(FASTQ_BASENAME)_pbsim2_R94.fq \
 $(FASTQ_BASENAME)_pbsim2_R95.fq \
 $(FASTQ_BASENAME)_pbsim2_R103.fq \
 $(FASTQ_BASENAME)_pbsim2_P6C4.fq

$(FASTQ_BASENAME)_pbsim2_%.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) pbsim2 -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@))) -d $(DEPTH_TSV) -j 1200


.PHONY:BADREAD_FQ
BADREAD_FQ: \
$(FASTQ_BASENAME)_badread_nanopore2018.fq \
$(FASTQ_BASENAME)_badread_nanopore2020.fq \
$(FASTQ_BASENAME)_badread_pacbio2016.fq \
$(FASTQ_BASENAME)_badread_verybad.fq \
$(FASTQ_BASENAME)_badread_verynice.fq \


$(FASTQ_BASENAME)_badread_%.fq: $(SEL_CDNA_FASTA) $(DEPTH_TSV)
	$(YASIM) badread -F "$<".d -o $(basename $@) -m $(lastword $(subst _,\ ,$(basename $@))) -d $(DEPTH_TSV)

# $(FASTQ_BASENAME)_simlord.fq: $(SEL_CDNA_FASTA)
# 	$(YASIM) simlord -F "$<".d -o $(basename $@)
# 	# rm -rf $(basename $@).d

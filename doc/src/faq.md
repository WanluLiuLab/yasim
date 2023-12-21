# FAQs

## Experimental Design

- **Q: Can YASIM be used to benchmark tools that identify Differentially Expressed Genes (DEGs)?**

    A: No. Generated depth of two runs are totally independent without biological assumptions held by tools like [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

---

- **Q: Can I use YASIM to generate scRNA-Seq data?**

    A: You might be interested in [scTCR-Seq variant](https://github.com/WanluLiuLab/yasim-sctcr.git) of YASIM, which contains a module that generates scRNA-Seq data.

---

- **Q: Can I use custom count matrices in YASIM?**

    A: Yes. However, please make sure that:

    - Your count matrix is isoform-level. If not, you may either pass it to `generate_isoform_depth` or re-quantify your sample in isoform-level.
    - Your count matrix is a Tab-Separated Value (TSV) which consists of 2 rows named `TRANSCRIPT_ID` and `DEPTH`.
    - Your count matrix uses same trancript ID nomenclature. If not, please convert it using tools like [BiomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).

## Guide to Parameters

- **Q: How do I specify read completeness?**

    A: Use `--truncate_ratio_3p` or `--truncate_ratio_5p` to truncate a portion of every read when invoking the LLRGs. Note, the 3' and 5' ends are in respect to the sequencer instead of cDNA or reference genome.

---

- **Q: How should I decide on the complexity of AS events?**

    A: Use `--complexity` while using `generate_as_events`. The larger for more complex genome annotation (i.e., larger mean number of transcripts inside a gene).

---

- **Q: How should I choose the sample-level sequencing depth?**

    A: To adjust the sequencing depth for the entire sample, use `--depth` while invoking `generate_gene_depth`. Note that the depth is not equal to metrics like read counts, RPKM, TPM, etc.

    If you want to customize sequencing depth for individual genes or individual isoform inside a gene, there's currently no way to do it. You may manually edit the depth TSV produced using an editor.

---

- **Q: How do I set low- and high-cutoffs for sequencing depth?**

    A: You are allowed to set low- and high-cutoffs in `generate_gene_depth` and `generate_isoform_depth`. The lower cutoff is designed for LLRGs that do not perform well if given a depth smaller than 1. Setting this value higher will result in a smoother execution (i.e., less number of errors) but will impair fold changes.

    The higher cutoff ratio is designed to reduce extremely high expression levels. For example, the GMM may generate large numbers like 100,000 (which would take LLRGs years to finish) even if your targeted sample-level mean depth is only 20. The real cutoff will be higher cutoff ratio times targeted sample-level mean depth. While using CCS-based LLRGs, you are recommended to lower this value since calculation of consensus sequence is computational-intensive.

    Both values affect the maximal fold-change of a sample. For example, you set low-cutoff to 0.1, mean depth to 50 and high cutoff ratio to 100. The generated expression levels should at a range of 0.1 to 5,0000 which allows a 50,000-fold change between most and least expressed genes.

## Compatibility

- **Q:What kind of reference genome would YASIM accept?**

    A: For FASTA \& GTF, this simulator does not support fancy genomic features like decoy sequences, HLA sequences, EBV sequences, fix or novel patches, unplaced or unlocalized scaffolds. You are free to use references with those features, but generated data may not be biologically meaningful.

    For GTF, you shold use UCSC references **WITHOUT** being processed (e.g., sorted by `bedtools`) for compatibility. Reference genome annotation from NCBI RefSeq official website (i.e., <https://www.ncbi.nlm.nih.gov/genome/?term=txid6239[orgn]>) is explicitly incompatible, and those from Ensembl may experience bugs.
---

- **Q: Since the `transcribe` step generates cDNA, how is it different from cDNA we downloaded from Ensembl?**

    A: Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.

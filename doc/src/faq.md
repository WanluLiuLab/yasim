# Frequently Asked Questions (FAQs)

- **Q:What kind of reference genome should be used?**

    A: For FASTA \& GTF, this simulator does not support fancy genomic features like decoy sequences, HLA sequences, EBV sequences, fix or novel patches, unplaced or unlocalized scaffolds. You are free to use references with those features but generated data may not be biologically meaningful.

    For GTF, you shold use UCSC references **WITHOUT** being processed (e.g., sorted by `bedtools`) for compatibility. Reference genome annotation from NCBI RefSeq official website (i.e., <https://www.ncbi.nlm.nih.gov/genome/?term=txid6239[orgn]>) is explicitly incompatible, and those from Ensembl may experience bugs.

---

- **Q: Can YASIM be used to benchmark tools that identifies Differentially Expressed Genes (DEGs)?**

    A: No. Generated depth of two runs are totally independent without biological assumptions held by tools like [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

---

- **Q: Since the `transcribe` step generates cDNA, how is it different from cDNA we downloaded from Ensembl?**

    A: Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.

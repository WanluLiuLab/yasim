# FAQs

## Experimental Design

- **Q: Can YASIM be used to benchmark tools that identifies Differentially Expressed Genes (DEGs)?**

    A: No. Generated depth of two runs are totally independent without biological assumptions held by tools like [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

---

- **Q: Can I use YASIM to generate scRNA-Seq data?**

    A: You might be interested in [scTCR-Seq variant](https://github.com/WanluLiuLab/yasim-sctcr.git) of YASIM, which contains a module that generates scRNA-Seq data.

## Guide to Parameters

- **Q: How do I specify read completeness?**

    A: Use `--truncate_ratio_3p` or `--truncate_ratio_5p` to truncate a portion of every read when invoking the LLRGs. Note, the 3' and 5' ends are in respect to the sequencer instead of cDNA or reference genome.

---

- **Q: How should I decide on the complexity of AS events?**

    A: Use `--complexity` while generating AS events. The larger for more complex genome annotation (i.e., more number of transcripts inside a gene).

---

- **Q: How should I choose the Sequencing Depth of Isoforms?**

    A: Use `--depth` while generating gene-level depth file instead of isoform-level depth file.

## Compatibility

- **Q:What kind of reference genome would YASIm accept?**

    A: For FASTA \& GTF, this simulator does not support fancy genomic features like decoy sequences, HLA sequences, EBV sequences, fix or novel patches, unplaced or unlocalized scaffolds. You are free to use references with those features but generated data may not be biologically meaningful.

    For GTF, you shold use UCSC references **WITHOUT** being processed (e.g., sorted by `bedtools`) for compatibility. Reference genome annotation from NCBI RefSeq official website (i.e., <https://www.ncbi.nlm.nih.gov/genome/?term=txid6239[orgn]>) is explicitly incompatible, and those from Ensembl may experience bugs.
---

- **Q: Since the `transcribe` step generates cDNA, how is it different from cDNA we downloaded from Ensembl?**

    A: Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.

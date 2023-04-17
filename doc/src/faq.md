# Frequently Asked Questions (FAQs)

## What kind of reference genome should be used?

For FASTA \& GTF, this simulator does not support fancy genomic features like decoy sequences, HLA sequences, EBV sequences, fix or novel patches, unplaced or unlocalized scaffolds. You are free to use references with those features but generated data may not be biologically meaningful.

For GTF, you shold use UCSC references **WITHOUT** being processed (e.g., sorted by `bedtools`) for compatibility. Reference genome annotation from NCBI RefSeq official website (i.e., <https://www.ncbi.nlm.nih.gov/genome/?term=txid6239[orgn]>) is explicitly incompatible, and those from Ensembl may experience bugs.

## How is YASIM V3API differes from V2API?

The YASIM V3 and V2 API differes in depth generation. In YASIM V2, we generate isoforms regardless of base expression level of genes. In YASIM V3, however, we first assign base level of each gene and then generate isoform-level abundance upon them. The YASIM V2 API also does not allow zeros to be present.

The YASIM V2 API was preserved. You can access it through `generate_depth_v2`.

## Can YASIM be used to benchmark tools that identifies Differentially Expressed Genes (DEGs)?

No. Generated depth of two runs are totally independent without biological assumptions held by tools like [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Since the `transcribe` step generates cDNA, how is it different from cDNA we downloaded from Ensembl?

Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.

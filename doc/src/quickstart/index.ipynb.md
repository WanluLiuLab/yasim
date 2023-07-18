---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Quickstart

This is a quickstart documentation for YASIM. In this documentation, you would generate Third-generation Sequencing (TGS) RNA-Seq reads from *C. Elegans* (worm) reference genome with Alternative Splicing (AS) events.

This tutorial assumes a basic understanding of RNA-Seq and Shell scripting.

+++

**How to read this documentation**:

Here would list version information of each component used in this tutorial for reproductive purposes.

| Software | Version   |
|----------|-----------|
| GNU Bash | 5.2.15(1) |
| GNU Grep | 3.8       |
| GNU Wget | 1.21.3    |
| pbsim    | 1.0.4     |
| yasim    | 3.1.5     |

## YASIM Data Flow Diagram

The following diagram lists the full YASIM workflow. You are not limited to this and can start at any position.

```{figure} yasim_toplevel.svg
:alt: YASIM data flow
:width: 100%
```

+++

## Step 0. Retrieve Reference Genome Sequence and Annotation

Inside the example, chrI of *C. Elegans* reference genome sequence (in FASTA format) and annotation (in GTF format) from UCSC will be used.

Get chrI of the CE11 reference genome sequence and annotation from UCSC.

```shell
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
gunzip ce11.ncbiRefSeq.gtf.gz
grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ce11.fa.gz
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa
```

## Step 1. Generation of AS Events

This step would generate alternative splicing events. It would take reference genome annotation as input and generate GTF with AS events as output.

The following example generates AS events from chromosome 1 of CE11 reference genome annotation with Transcriptome Complexity Index 5:

```shell
python -m yasim generate_as_events \
    -f ce11.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq.chr1.as.gtf \
    -c 5
```

**Generates:**

- `ce11.ncbiRefSeq.chr1.as.gtf`, the generated GTF with AS events.
- `ce11.ncbiRefSeq.chr1.gtf.0.4.gvpkl.xz`, if not exist. This is a cache file for the `labw_utils` GTF parser.

The generated GTF should be seen as ground truth for benchmarking AS detectors. The reference genome should have a Transcriptome Complexity Index between 1 and 2 with real organisms in about 2.

+++

## Step 2. Generate Sequencing Depth of Gene

This step would generate base gene expression levels in coverage (aka., depth) for each gene over some GTF.

```{warning}
The generated coverage is **NOT** number of reads generated! It cannot be used as ground truth to assess quantification software! The number of reads ground truth will be provided by LLRG UIs introduced below.
```

Example of coverage generation on genome annotation with AS events just generated with mean depth 5:

```shell
python -m yasim generate_gene_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_gene_depth.tsv \
    -d 5
```

**Generates:**

- `ce11_gene_depth.tsv`, a TSV file with the following columns:
  - `GENE_ID`, the `gene_id` field in GTF.
  - `DEPTH`, gene expression level in coverage.
- `ce11.ncbiRefSeq.chr1.as.gtf.0.4.gvpkl.xz`, if not present.

+++

## Step 3. Generate Sequencing Depth of Isoform

Here assigns expression level in coverage to isoforms. The mean expression level of each isoform in some genes should be equal to the base expression level assigned to that gene in the previous step.

Following is a typical example. The parameter `alpha` (default to 4) was used to adjust the evenness between isoform expression levels inside a gene.

```shell
python -m yasim generate_isoform_depth \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -d ce11_gene_depth.tsv \
    -o ce11_isoform_depth.tsv
```

**Generates:**

- `ce11_isoform_depth.tsv`, a TSV file with the following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `DEPTH`, isoform expression levels in coverage.

+++

## Step 4. Transcribe GTF to FASTA

This step is general-purpose. It can be used to transcribe (**stranded**) any GTF that contains some isoforms into a FASTA file of all cDNA and a directory with all cDNAs in separate FASTA files and skip those isoforms whose region is not defined in the genomic sequence (FASTA). It would not add post-transcriptional modifications.

For those who are familiar with [BedTools](https://bedtools.readthedocs.io), It should generate similar output with:

``` shell
bedtools getfasta -nameOnly -s -fi [FASTA] -bed [GTF] > [OUT]
```

Example:

```shell
python -m labw_utils.bioutils transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.as.gtf \
    -o ce11_trans_as.fa
```

**Generates:**

- `ce11_transcripts.fa`, the generated cDNA sequence FASTA.
- `ce11_transcripts.fa.d`, the directory where every cDNA is stored as separate FASTA files.
- `ce11_transcripts.fa.stats`, a TSV file with the following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `GENE_ID`, the `gene_id` field in GTF.
  - `SEQNAME`, chromosome \& scaffold \& contig name.
  - `START`, the `start` field in GTF, 1-based inclusive.
  - `END`, the `end` field in GTF, 1-based inclusive.
  - `STRAND`, the `strand` field in GTF.
  - `ABSOLUTE_LENGTH`, is `START` - `END` + 1.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `GC`, GC content of the cDNA in percentage.

+++

## Step 5. Invocation of LLRGs: Use PBSIM for Example

The Low-Level Read Generators (LLRGs) are programs that simulate DNA-Seq on some reference genome sequences by clipping reads in appropriate lengths and introducing sequencing errors. YASIM invokes LLRG on stranded cDNA sequences to generate RNA-Seq data. Here we would demonstrate their usage with the following examples:

```{warning}
The official build of PBSIM, PBSIM2 and PBSIM3 shares a common executable anme (`pbsim`) but with different argument layout. For convenience, I renamed executable of PBSIM2 to `pbsim2` and PBSIM3 to `pbsim3`. If you do not use this in your computer, please use the `-e` option.
```

PBSIM is a general-purpose TGS DNA-Seq simulator that supports the simulation of PacBio RS sequencers. It can generate both Continuous Long Reads (CLR) and Circular Consensus Sequence (CCS) data.

Compared to NGS simulators, TGS simulators have `truncate_ratio_3p` and `truncate_ratio_5p`. These two parameters are used to set hard limits at two sides that allow the simulation of incomplete reads due to reasons like 3' truncation.

Following is an example of a simulation of CCS data using the PacBio RS model:

```shell
python -m yasim pbsim \
    -F ce11_trans_as.fa.d \
    -j 40 \
    --ccs \
    -d ce11_isoform_depth.tsv \
    -o pbsim_mode
```

**Generates:**

- `pbsim_mode.fq`, simulated Single-End FASTQ.
- `pbsim_mode.d`, a temporary directory for diagnostic purposes that can be safely deleted.
- `pbsim_mode.fq.stats`, statistics of simulated FASTQ. a TSV containing the following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `INPUT_DEPTH`, isoform expression levels in coverage provided by the upstream source.
  - `SIMULATED_N_OF_READS`, simulated number of reads. **This value can be used in assessing read quantifiers**.
  - `SIMULATED_N_OF_BASES`, simulated number of bases.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `SIMULATED_DEPTH`, simulated depth.

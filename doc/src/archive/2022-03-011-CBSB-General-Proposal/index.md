---
orphan: true
---

# Detailed CBSB Project Proposal

2022-03-11, Zhejian YU, Yaqi SU, Ruihong YUAN

## Introduction

Alternative splicing (AS) is a molecular mechanism that modifies pre-mRNA constructs before translation. This process can produce a diversity of mRNAs from a single gene by arranging coding sequences (exons) from recently spliced RNA transcripts into different combinations. The mechanisms of AS help to explain how one gene can be encoded into numerous proteins with various functions. This complexity helps drive the cellular differentiation and diversity observed throughout biology.

RNA-Seq is commonly used for AS identification and numerous tools are available for this kind of analysis. To objectively assess their sensitivity and specificity together with scalability, gold standard datasets are required. Since the true number and types of AS events (such as exon skipping, intron retention, etc.) in an experimental dataset is not known, suitable datasets need to be simulated. Moreover, Third-Generation Sequencing (TGS), which was developed in recent years, allows sequencing of single transcripts without splitting into short reads. With this method, the transcribed exons can be better preserved, allowing us to perform AS analysis more accurately.

This project is designed to write a simulator that generates artificial novel alternative splicing (AS) events and real gene expression patterns in universal organisms for the purpose of profiling tools that is claimed to be able to detect AS events and expression levels.

## What is Done and What is Left

The Infrastructure Construction like GTF Parser, FASTA/FASTQ Parser, Gene-Transcript-Exon Three-Tier Abstraction ahd almost done.

What is left is mainly the CBSB [^CBSB] part, that is, to fit distributions.

[^CBSB]: Computational Biology \& Systems Biology is a course at University of Edinburgh.


# L4 C. Elegans RNA-seq Analysis


## The Data

### SRP095467: L4, day 1, day 7 and day 10 old total RNA-Seq from Caenorhabditis elegans

Using Illumina HiSeq 2500

SRR5123644 SRR5123648 SRR5123649

### ERP114391: Nanopore based direct RNA sequencing across development of C elegans.

Using Oxford Nanopore GridION

ERR3245471 ERR3245470

### SRP185743: Full-length mRNA sequencing reveals principles of poly(A) tail length control

Using PacBio SMRT Sequel

SRR8568877 SRR8568878

## Get Data from EBI:

```shell
cat << EOF | while read line; do [ ! -f $(basename "${line}") ] && axel -n 20 "${line}"; done
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/001/ERR3245471/ERR3245471.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/000/ERR3245470/ERR3245470.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/007/SRR8568877/SRR8568877_subreads.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/008/SRR8568878/SRR8568878_subreads.fastq.gz
EOF
```

## Align to Transcriptome

The reason why we align to transcriptome instead ofgenome is because using this method we may get coverage from aligned data without transformation.

Using minimap2 2.17-r941, bwa 0.7.17-r1188, samtools 1.11 (using htslib 1.11-4)

```shell
bwa index ce11.reference_transcripts.fa -p CE11_TRANSCRIPTOME_BWA_INDEX
minimap2 -d ce11.reference_transcripts.fa.mmi ce11.reference_transcripts.fa -t 16

minimap2 -x splice -a -t 16 ce11.reference_transcripts.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz | \
samtools sort -@ 16 -o nanopore_transcriptome_1.bam # 84.24% mapping rate
minimap2 -x splice -a -t 16 ce11.reference_transcripts.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz | \
samtools sort -@ 16 -o pacbio_transcriptome_1.bam # 96.65% mapping rate

bwa mem -t 16 CE11_TRANSCRIPTOME_BWA_INDEX <(zcat SRR5123644_1.fastq.gz SRR5123648_1.fastq.gz SRR5123649_1.fastq.gz) \
<(zcat SRR5123644_2.fastq.gz SRR5123648_2.fastq.gz SRR5123649_2.fastq.gz) |\
samtools sort -@ 16 -o ngs_transcript.bam

for fn in *.bam; do samtools index $fn; samtools depth $fn > $fn.depth.tsv; done
```


## Merge results

```shell
Rscript R/dge_merge.R \
--fa_stats suppl_data_analysis/c_elegans_transcriptome/ce11.reference_transcripts.fa.stats \
--featureCounts_tsv suppl_data_analysis/c_elegans_transcriptome/pacbio_1.bam.tsv suppl_data_analysis/c_elegans_transcriptome/nanopore_1.bam.tsv \
--salmon_quant_sf suppl_data_analysis/c_elegans_transcriptome/ngs_1_quant.sf \
-o suppl_data_analysis/c_elegans_transcriptome/all \
--libfile R/lib.R
```

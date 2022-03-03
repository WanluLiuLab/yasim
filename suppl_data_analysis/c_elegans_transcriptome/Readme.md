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

```shell
bwa index ce11_reference_transcript.fa -p CE11_TRANSCRIPTOME_BWA_INDEX
minimap2 -d ce11_reference_transcript.fa.mmi ce11_reference_transcript.fa -t 16

minimap2 -a -t 16 ce11_reference_transcript.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz | \
samtools sort -@ 16 -o nanopore_transcriptome_1.bam
minimap2 -a -t 16 ce11_reference_transcript.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz | \
samtools sort -@ 16 -o pacbio_transcriptome_1.bam

bwa mem -t 16 CE11_TRANSCRIPTOME_BWA_INDEX <(zcat SRR5123644_1.fastq.gz SRR5123648_1.fastq.gz SRR5123649_1.fastq.gz) \
<(zcat SRR5123644_2.fastq.gz SRR5123648_2.fastq.gz SRR5123649_2.fastq.gz) |\
samtools sort -@ 16 -o ngs_transcript.bam

samtools index *.bam

for fn in *.bam; do samtools depth $fn > $fn.depth.tsv; done
```



## Quantification using Salmon 1.4.0:

```shell
salmon quant --gcBias -i SALMON_CE11_INDEX --threads 16 -l IU -o salmon_quant \
-1 SRR5123644_1.fastq.gz SRR5123648_1.fastq.gz SRR5123649_1.fastq.gz \
-2 SRR5123644_2.fastq.gz SRR5123648_2.fastq.gz SRR5123649_2.fastq.gz
```

## Alignment using minimap2

```shell
minimap2 -d ce11.fa.mmi ce11.fa

minimap2 -a -x splice -t40 ce11.fa.mmi ERR3245471.fastq.gz ERR3245470.fastq.gz | samtools sort -@ 16 -o nanopore_1.bam
minimap2 -a -x splice -t40 ce11.fa.mmi SRR8568877_subreads.fastq.gz SRR8568878_subreads.fastq.gz | samtools sort -@ 16 -o pacbio_1.bam
samtools index *.bam
```

Output of `samtools flagstat nanopore_1.bam`:

```text
557297 + 0 in total (QC-passed reads + QC-failed reads)
16942 + 0 secondary
157 + 0 supplementary
0 + 0 duplicates
473439 + 0 mapped (84.95% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Output of `samtools flagstat pacbio_1.bam`

```text
883105 + 0 in total (QC-passed reads + QC-failed reads)
54767 + 0 secondary
251714 + 0 supplementary
0 + 0 duplicates
851954 + 0 mapped (96.47% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```


## Quantification using featureCounts

```shell
for bam_file in *.bam; do
  featureCounts -t transcript -g transcript_id -a ce11.ncbiRefSeq.gtf -L -o "${bam_file}".tsv "${bam_file}"
done
```


## Merge results

```shell
Rscript R/dge_merge.R \
--fa_stats suppl_data_analysis/c_elegans_transcriptome/ce11_reference_transcript.fa.stats \
--featureCounts_tsv suppl_data_analysis/c_elegans_transcriptome/pacbio_1.bam.tsv suppl_data_analysis/c_elegans_transcriptome/nanopore_1.bam.tsv \
--salmon_quant_sf suppl_data_analysis/c_elegans_transcriptome/ngs_1_quant.sf \
-o suppl_data_analysis/c_elegans_transcriptome/all \
--libfile R/lib.R
```

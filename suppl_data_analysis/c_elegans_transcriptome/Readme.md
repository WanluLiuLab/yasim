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

## Align to Transcriptome

The reason why we align to transcriptome instead ofgenome is because using this method we may get coverage from aligned data without transformation.

Using minimap2 2.17-r941, bwa 0.7.17-r1188, samtools 1.11 (using htslib 1.11-4)

The alignment details can be found at `get_data.sh`.

## Merge Results


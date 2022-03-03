# C. Elegans RNA-seq Analysis

## Get Data from SRA:

```shell
cat << EOF | while read line; do axel -n 20 "http://${line}"; done
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/002/SRR5123642/SRR5123642_1.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/002/SRR5123642/SRR5123642_2.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_1.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/009/SRR5123649/SRR5123649_2.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_1.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/004/SRR5123644/SRR5123644_2.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_1.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/008/SRR5123648/SRR5123648_2.fastq.gz
EOF
```

## Quantification using Salmon 1.4.0:

```shell
salmon quant -i SALMON_CE11_INDEX --threads 16 -l IU -1 *_1.fastq.gz -2 *_2.fastq.gz -o salmon_quant
```

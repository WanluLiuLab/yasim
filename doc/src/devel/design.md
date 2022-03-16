# Design of YASIM

## Aim

This project is designed to write a simulator that generates artificial novel alternative splicing (AS) events and real gene expression patterns in universal organisms for the purpose of profiling tools that is claimed to be able to detect AS events and expression levels.

## Abbreviations

AS
: Alternative Splicing

CLI
: Command Line Interface

GEP
: Gene Expression Profile

LLRG
: Low-Level Read Generators

nReads
: Number of Reads

refGTF
: Reference GTF

RPKM
: Read Per Kilo base per Million mapped reads

TPM
: Transcript per Million

## Definitions and Terminologies

Here lists definitions and terminologies used in API docs, etc.

LLRG
: A program that simulates the behavior of real sequencers. It can be abstracted into a step, which takes input sequence, sequencing depth, some error profiles, and generates simulated sequencing data with its statistical parameters (e.g., read length, truncation, GC bias, etc.) like data generated by real sequencers.

LLRG Python Adapter
: A Python class that wraps LLRG, or LLRG Shell Adapter, to integrate it into the main simulation process. It is called by LLRG UI (introduced below).

LLRG Shell Adapter
: A Shell script that handles installation and maintenance of corresponding LLRG.

LLRG User Interface (LLRG UI)
: A CLI provided by YASIM, allowing users to produce simulated reads using their preferred LLRGs.

Absolute Abundance ($Y$) of a Transcript/Gene
: RPKM/TPM/nReads/ Coverage of a gene or transcript. We denote abundance of transcript $t$ as $Y_{t}$.

GEP
: A profile with Transcript ID - Absolute Abundance of that Transcript ($t_{ID} \rightarrow (Y_{t})$) pair.

Ground Truth GEP TSV
: A TSV that stores ground truth of GEP in **nReads**.

refGTF
: Reference GTF can be downloaded at UCSC, NCBI, Ensembl or WormBase (for _C. Elegans_). Support of GFF3 format is purposed **but not implemented**. Following is a table where you can download refGTF:

| Name    | Genome FASTA                                                 | Genome GTF                                                   | cDNA FASTA                                                   |
| ------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| UCSC    | [hg38.p13](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | [ncbiRefSeq](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz) | NA                                                           |
| Ensembl | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz) | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) |
| NCBI    | [GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz) | [GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gtf.gz) | NA                                                           |

AS Ground Truth GTF
: GTF with AS events, which is generated from refGTF by YASIM.

```{caution}
Shouls we filter:

- Alternative Regopns `alt`, Patches `patch` and other non-main-chromosome contigs?
- None Protein-Coding RNAs?
```

# Steps of Simulation

## `generate_as_events`: Simulate AS Events based on a Statistical Model Given by Mining Raw Data and Works of Literature

- INPUT: refGTF.
- OUTPUT: AS Ground Truth GTF.

```{caution}
This part is waiting for satistical tests.
```

Firstly, we parse the GTF and get a list of AS-able transcripts. Simulated from pre-fitted distribution of corresponding event type, AS events are then added to these transcripts.

## `generate_depth`: Simulate GEP based on a Statistical Model Given by Mining Raw Data and Works of Literature

- INPUT: AS Ground Truth GTF.
- OUTPUT: A `TranscriptID -  Coverage` TSV.

```{caution}
This part is waiting for satistical tests.
```

From literature of MicroArray [^Furusawa2003] [^Ueda2004] [^Lu2005] [^Brown2007], we can infer that the distribution of GEP of general organisms obeys Zipf's Law [^Zipf1949] and can be fitted into a Log-Normal Distribution. This opinion was challenged by [^Nowrousian2013].

The reason why we use "coverage" instead of more formal parameters like RPKM or TPM is that most LLRGs (introduced below) are designed for DNA-Seq simulation and accept coverage only (a minority also accepts a number of reads).

```{warning}
The generated coverage is **NOT** number of reads generated! It cannot be used as ground truth to assess quantification software! The number of reads ground truth will be provided by LLRG UIs introduced below.
```

## `transcribe`: Transcribe GTF into AS cDNA

- INPUT: AS Ground Truth GTF
- OUTPUT: Transcribed AS cDNA

This step is a general-purposed transcriber that can be used to convert any GTF that contains some transcript into a FASTA of all cDNA and a directory with all cDNAs in separate FASTA.

The seqname of FASTA records is `transcript_id` of AS Ground Truth GTF.

```{note}
Although this software can be used to generate reference cDNAs for software like Salmon, there are differences between transcribed cDNA and Ensembl-provided cDNA. Ensembl-provided cDNA does not include small features like lncRNA, while `yasim` transcribed cDNA includes all transcripts and exons inside GTF.
```

## Use LLRGs to Generate Raw Reads for Next- or Third-Generation Sequencing

- INPUT: A `TranscriptID -  Coverage` TSV, reference FASTA, AS Ground Truth GTF
- OUTPUT: LLRG-Generated FASTQ, Ground Truth GEP TSV

This step calls LLRG UIs to "sequence" simulated cDNAs to simulated reads with depth. A ground truth including the actual number of reads generated will be produced.

The LLRGs are all from third-parties and what YASIM do is to write Shell \& Python adapters for them. The adapters are designed in a way that new LLRGs can be added only by inheriting from some base classes.

# Core Data Analysis

The data analysis in this part is mainly _C. Elegans_ Data. From previous studies, we have collected three L4 _C. Elegans_ RNA-Seq data. They are ONT dRNA (NANOPORE), PacBio cDNA (PACB), and Illumina HiSeq 2500 cDNA (3 replicates, ILLM1, ILLM2, ILLM3 respectively). Followings are how they will be analysed:

## GEP

From pior knowledge, we know that:

1. The GEP across transcripts in one sample can be fitted into a log-normal (LNorm) distribution.

    $$
    Y_{t, s} \sim LNorm(mean=\mu_{s}, variance=\sigma_{s})
    $$

    Where Y is absolute abundance of a transcript, $t$ is a transcript ID, $s$ is s sample ID.

2. The GEP across transcripts in one sample obeys Ziff's Law.

    $$
    R(Y_{t, s}) = Y_{t, s}^{-k}
    $$

    which means:

    $$
    \log(R(Y_{t, s})) \approx -k * \log(Y_{t, s})
    $$

    where $R$ is the rank of abundance $Y$.

    This is given in a new form in FluxSimulator [^Griebel2012]:

    $$
    Y_{t} = \max_{t}(Y_{t}) R(Y_{t})^{k}\exp(\frac{R(Y_{t})}{a}\left(\frac{R(Y_{t})}{b}\right)^2)
    $$

    Where $k$, $a$ and $b$ are parameters with $k \in (‐0.6, ­‐0.9 )$ and $a=b$ and $a \sim 10^{4}$.

    ```r
    dFluxSim <- function(R_Y, k, Y_max, a, b){
        return (Y_max * R_Y^k*exp(Y/a*(Y/b)^2))
    }
    ```

3. Given the relative abundance of a transcript, its abundance across samples is negative binomial (NB) distributed [^Anders2010] [^Robinson2010].

    $$
    Y_{t, s} \sim NB(mean = \mu_{t}, size = r_{t})
    $$

So, we hypothesis that:

$$
Y_{t, s} \sim SomeDistribution(GC_{t}, LEN_{t}, POS_{t}, mean = \mu, \ldots)
$$

## AS Events

From [ASimulatoR](https://github.com/biomedbigdata/ASimulatoR) [^Quirin2021] we get that:

$$
E_{c, t, s} \sim SomeKnownDistribution(some\_mean = \mu, some\_se = \sigma)
$$

Where $E_{c}$ is AS Event type $c$.

[^Furusawa2003]: Furusawa, C., & Kaneko, K. (2003). Zipf's law in gene expression. Physical review letters, 90(8), 088102. <https://doi.org/10.1103/PhysRevLett.90.088102>
[^Ueda2004]: Ueda, H. R., Hayashi, S., Matsuyama, S., Yomo, T., Hashimoto, S., Kay, S. A., Hogenesch, J. B., & Iino, M. (2004). Universality and flexibility in gene expression from bacteria to human. Proceedings of the National Academy of Sciences of the United States of America, 101(11), 3765–3769. <https://doi.org/10.1073/pnas.0306244101>
[^Lu2005]: Lu, T., Costello, C. M., Croucher, P. J., Häsler, R., Deuschl, G., & Schreiber, S. (2005). Can Zipf's law be adapted to normalize microarrays?. BMC bioinformatics, 6, 37. <https://doi.org/10.1186/1471-2105-6-37>
[^Nowrousian2013]: Nowrousian M. (2013). Fungal gene expression levels do not display a common mode of distribution. BMC research notes, 6, 559. <https://doi.org/10.1186/1756-0500-6-559>
[^Zipf1949]: Zipf, G.K. (1949). Human Behaviour and the Principle of Least Effort: an Introduction to Human Ecology. <https://doi.org/10.2307/2572028>
[^Brown2007]: Brown, R. J. C. (2007). The use of Zipf’s law in the screening of analytical data: A step beyond Benford. Analyst, 132(4), 344–349. <https://doi.org/10.1039/B618255K>
[^Anders2010]: Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. Genome biology, 11(10), R106. <https://doi.org/10.1186/gb-2010-11-10-r106>
[^Robinson2010]: Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics (Oxford, England), 26(1), 139–140. <https://doi.org/10.1093/bioinformatics/btp616>
[^Griebel2012]: Griebel, T., Zacher, B., Ribeca, P., Raineri, E., Lacroix, V., Guigó, R., & Sammeth, M. (2012). Modelling and simulating generic RNA-Seq experiments with the flux simulator. Nucleic acids research, 40(20), 10073–10083. <https://doi.org/10.1093/nar/gks666>
[^Quirin2021]: Quirin Manz, Olga Tsoy, Amit Fenn, Jan Baumbach, Uwe Völker, Markus List, Tim Kacprowski, ASimulatoR: splice-aware RNA-Seq data simulation, Bioinformatics, Volume 37, Issue 18, 15 September 2021, Pages 3008–3010, <https://doi.org/10.1093/bioinformatics/btab142>
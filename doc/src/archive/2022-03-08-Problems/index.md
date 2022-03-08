# Core Problems and Tasks before Project

2022-03-08, Zhejian YU, Yaqi SU, Ruihong YUAN

## Tasks needed to be Done

### Find Datasets

- Analyse _C. Elegans_ Data. From previous studies, we have collected three L4 _C. Elegans_ RNA-Seq data. They are:  ONT dRNA (NANOPORE), PacBio cDNA (PACB) and Illumina HiSeq 2500 cDNA (3 replicates, ILLM1, ILLM2, ILLM3 respectively).
- Analyse Human RNA-Seq Data. Data collection needed.

### Analyse DGE Data

We hypothesis that:

$$
Y_{t, s} \sim SomeDistribution(GC_{t}, LEN_{t}, POS_{t}, mean = \mu, \ldots)
$$

Where $Y$ is transcript abundance (Coverage) ,$t$ is Transcript ID, $s$ is Sample ID, $GC$ is GC content of given transcript, $LEN$ is length of given transcript and $POS$ is whether given transcript is at special positions (e.g., terminal region). $\mu$ is eman coverage accross this simulation and $\ldots$ measn other arbitrary arguments, like dispersion parameter, etc.

This distribution is designed to fit any arbitrary transcript or sample of specific organism. That is, for any case with known $GC$, $LEN$, $POS$ and setted $\mu$, the distribution should be transcript and sample irrelevant.

Here are known:

$$
Y_{t, s} \sim NB(mean = \mu_{t}, size = r_{t})
$$

That is, given relative abundance of a transcript, its abundance accross samples are negative bionomial distributed.

Here are the problems:

1. How to find such a distribution? A probable way is to fit all known distributions, and pick the one with highest Akaike Information Criterion (AIC).
2. We hypothesis that $GC_{t}$, $LEN_{t}$ and $POS_{t}$ is related to abundance. How to prove it? How to prove their independence (or, not independent)?

With minor problems like:

- How to get coverage. Two methods are available:
  - Align reads to genome using **spliced** aligner (STAR and minimap2 with `-x splice`), then use featureCounts to determine number of mapped reads. Coverage can be calculated by: $$Y_{t} = \bar{L} * N_{t} / LEN_{t}$$, where $\bar{L}$ is mean read length.
  - Align reads to transcriptome using **unspliced** aligner (BWA and minimap2), then use `samtools depth` to get per-base depth. Coverage can be calculated by: $$Y_{t} = \sum_{pos \in t} D_{pos} / LEN_{t}$$, where $pos$ is base position in transcript $t$ and $D_{pos}$ is depth created by `samtools depth`.

### Analyse AS Events Data

From ASimulatoR we get that:

$$
E_{c, t, s} \sim SomeKnownDistribution(some\_mean = \mu, some\_se = \sigma)
$$

Where $E_{c}$ is AS Event type $c$.

So we may hypothesis that:

$$
E_{c, t, s} \sim SomeDistribution(GC_{t}, LEN_{t}, POS_{t}, \ldots)
$$

Here are the problems:

1. How to get AS event for each transcript?
2. Other same satistical questions as is in DGE simulation.

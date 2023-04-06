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

# YASIM Tutorial

This is the tutorial documentation for YASIM. In this documentation, you would generate Next-Generation Sequencing (NGS) and Third-generation Sequencing (TGS) RNA-Seq reads from _C. Elegans_ (worm) reference genome with Alternative Splicing (AS) events, quantify them on an isoform level using Salmon and compare them with ground-truth data.

This tutorial assumes basic understandings on RNA-Seq, Shell scripting and Python programming with Pandas and Seaborn.

**How to read this documentation**: Code block without leading `!` are Python code blocks. For example:

```{code-cell}
import os; os.name
```

Code block with leading `!` are Shell code blocks. For example:

```{code-cell}
!ls -lFh | grep ipynb
```

```{figure} ../fig/yasim_toplevel.svg
:width: 100%
:align: left
:alt: Common YASIM Workflow
```

## Environment Specifications

Here would list version information of each component used in this tutorial for reproductive purposes.

+++

### Python Packages

```{code-cell}
import sys
import itertools

import numpy as np
import pandas as pd
try:
    import pyarrow as pa
    CSV_ENGINE = "pyarrow"
except ImportError:
    CSV_ENGINE = "c"
import seaborn as sns
import matplotlib.pyplot as plt

import yasim
import labw_utils

print("Python version: " + sys.version.replace("\n", ""))
print("Pandas version: " + pd.__version__)
print("Numpy version: " + np.__version__)
print("Seaborn version: " + sns.__version__)

print("YASIM version: " + yasim.__version__)
print("labw_utils version: " + labw_utils.__version__)
```

### Shell and Bioinformatics Utilities

Here we would assume that PBSIM3 is installed at `/home/yuzj/bin/pbsim3`.

```{code-cell}
!bash --version | head -n 1
!grep --version | head -n 1
!axel --version | head -n 1
!art_illumina 2>&1 | head -n 6
!ccs --version | head -n 1
!samtools --version | head -n 1
!salmon --version
```

## Preparations

Inside the example, chromosome 1 of _C. Elegans_ reference genome sequence (in FASTA format) and annotation (in GTF format) from UCSC is used.

```{warning}
This version of YASIM uses `labw_utils` 0.1.X GTF parser. This parser is NOT stable and is to be deprecated. You shold use **UNSORTED** UCSC references for compatibility. Reference genome annotation from NCBI RefSeq official website (i.e., <https://www.ncbi.nlm.nih.gov/genome/?term=txid6239[orgn]>) is explicitly incompatible, and those from Ensembl may experience bugs.
```

Get first chromosome of CE11 reference genome sequence and annotation from UCSC.

```{code-cell}
! if [ ! -f ce11.ncbiRefSeq.chr1.gtf ]; then \
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz; \
    gunzip ce11.ncbiRefSeq.gtf.gz; \
    grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf; \
else \
    echo "ce11.ncbiRefSeq.gtf.gz already exists."; \
fi
! if [ ! -f ce11.chr1.fa ]; then \
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz; \
    gunzip ce11.fa.gz; \
    head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa; \
else \
    echo "ce11.fa.gz already exists."; \
fi
```

## Generation of AS Events: `generate_as_events`

This step would generate alternative splicing events. It would take reference genome annotation as input and generate GTF with AS events as output.

The generated GTF should be seen as ground truth for benchmarking AS detectors. If you wish to benchmark quantifiers only, this step can be safely omitted.

The help message is as follows:

```{code-cell}
!python -m yasim generate_as_events --help
```

Example (Some warnings generated during GTF parsing was filtered out):

```{code-cell}
!if [ ! -f ce11.ncbiRefSeq.chr1.as.gtf ]; then \
    python -m yasim generate_as_events \
        -f ce11.fa \
        -g ce11.ncbiRefSeq.chr1.gtf \
        -o ce11.ncbiRefSeq.chr1.as.gtf \
        -c 5 \
    2>&1 | grep -v 'inferred from feature transcript'; \
else \
    echo "ce11.ncbiRefSeq.chr1.as.gtf already exists."; \
fi
```

Generates:

- `ce11.ncbiRefSeq.chr1.as.gtf`, the generated GTF with AS events.
- `ce11.ncbiRefSeq.chr1.gtf.0.4.gvpkl.xz`, if not exist. This is a cache file for the `labw_utils` GTF parser.

+++

Before proceeding, it would be desired to see difference between generated GTF and reference GTF. We would use `GeneView` class in `labw_utils` to read them into current Python environment:

```{code-cell}
from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory

as_gtf = GeneViewFactory.from_file("ce11.ncbiRefSeq.chr1.as.gtf")
ref_gtf = GeneViewFactory.from_file("ce11.ncbiRefSeq.chr1.gtf")
```

For a simple description:

```{code-cell}
print(f"AS GTF have {as_gtf.number_of_genes} Genes and {as_gtf.number_of_transcripts} Isoforms")
print(f"Reference GTF have {ref_gtf.number_of_genes} Genes and {ref_gtf.number_of_transcripts} Isoforms")

total_num_isoforms = len(set(
    itertools.chain(ref_gtf.iter_transcript_ids(), as_gtf.iter_transcript_ids())
))
num_new_isoforms = as_gtf.number_of_transcripts + ref_gtf.number_of_transcripts - total_num_isoforms
print(f"Generated {num_new_isoforms} new isoforms")
```

Now we would plot the distribution of Number of Isoforms per Gene (NIpG) for each GTF.

```{code-cell}
as_nipg = [gene.number_of_transcripts for gene in as_gtf.iter_genes()]
ref_nipg = [gene.number_of_transcripts for gene in ref_gtf.iter_genes()]

fig, axs = plt.subplots(nrows=2, sharex=True)
sns.violinplot(x=ref_nipg, ax=axs[0])
sns.violinplot(x=as_nipg, ax=axs[1])
axs[0].set_title("Reference GTF")
axs[1].set_title("Generated GTF")
```

It is clear that the generated GTF is more complex in terms of NIpG.

+++

## Generate Sequencing Depth of Gene: `generate_gene_depth`

This step would generate base gene expression level in coverage for each gene over some GTF.

```{note}
The YASIM V3 and V2 API differes in depth generation. In YASIM V2, we generate isoforms regardless of base expression level of genes. In YASIM V3, however, we first assign base level of each gene and then generate isoform-level abundance upon them. The YASIM V2 API also does not allow zeros to be present.

The YASIM V2 API was preserved. You can access it through `generate_depth_v2`.
```

```{warning}
The generated coverage is **NOT** number of reads generated! It cannot be used as ground truth to assess quantification software! The number of reads ground truth will be provided by LLRG UIs introduced below.
```

The help message is as follows:

```{code-cell}
!python -m yasim generate_gene_depth --help
```

Example:

```{code-cell}
!if [ ! -f ce11_gene_depth.tsv ]; then \
    python -m yasim generate_gene_depth \
        -g ce11.ncbiRefSeq.chr1.as.gtf \
        -o ce11_gene_depth.tsv \
        -d 5; \
else \
    echo "ce11_gene_depth.tsv already exists."; \
fi
```

Generates:

- `ce11_gene_depth.tsv`, a TSV file with following columns:
  - `GENE_ID`, the `gene_id` field in GTF.
  - `DEPTH`, gene expression level in coverage.
- `ce11.ncbiRefSeq.chr1.as.gtf.0.4.gvpkl.xz`, if not present.

Before proceed into next step, let's have a look at generated coverage file:

```{code-cell}
ce11_gene_depth = pd.read_table("ce11_gene_depth.tsv", engine=CSV_ENGINE)
```

```{code-cell}
ce11_gene_depth.head()
```

The depth distribution should satisfy a zero-inflated Gaussian Mixture Model (GMM) with its mean value similar to targeted mean.

```{code-cell}
sns.histplot(ce11_gene_depth, x="DEPTH")
```

```{code-cell}
ce11_gene_depth.DEPTH.mean() # Should be close to 5.
```

## Generate Sequencing Depth of Isoform: `generate_isoform_depth`

Here assigns expression level in coverage to isoforms. The mean expression level of each isoform in some gene should be equal to the base expression level assigned to that gene in previous step.

The help message is as follows:

```{code-cell}
!python -m yasim generate_isoform_depth --help
```

Example:

```{code-cell}
!if [ ! -f ce11_isoform_depth.tsv ]; then \
    python -m yasim generate_isoform_depth \
        -g ce11.ncbiRefSeq.chr1.as.gtf \
        -d ce11_gene_depth.tsv \
        -o ce11_isoform_depth.tsv \
        --alpha 4; \
else \
    echo "ce11_isoform_depth.tsv already exists."; \
fi
```

Generates:

- `ce11_isoform_depth.tsv`, a TSV file with following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `DEPTH`, isoform expression level in coverage.

+++

## Transcribe GTF to FASTA: `transcribe`

```{warning}
The `yasim transcribe` module was be deprecated. Use `labw_utils.bioutils transcribe` instead.
```

This step would transcribe the input genome GTF and genome FASTA into **stranded** transcriptome FASTA. It is designed to be general-purposed, i.e., can be applied on any matching GTF and FASTA. It should generate similar output with `bedtools getfasta -nameOnly -s -fi [FASTA] -bed [GTF] > [OUT]`

```{note}
Although this software can be used to generate reference cDNAs for software like Salmon, there are differences between transcribed cDNA and Ensembl-provided cDNA. Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.
```

The help message is as follows:

```{code-cell}
!python -m labw_utils.bioutils transcribe --help
```

Example:

```{code-cell}
!if [ ! -f ce11_trans_as.fa ]; then \
    python -m labw_utils.bioutils transcribe \
        -f ce11.chr1.fa \
        -g ce11.ncbiRefSeq.chr1.as.gtf \
        -o ce11_trans_as.fa; \
else \
    echo "ce11_trans_as.fa already exists."; \
fi
```

Generates:

- `ce11_transcripts.fa`, the generated cDNA sequence FASTA.
- `ce11_transcripts.fa.d`, the directory where every cDNA is stored as separate FASTA.
- `ce11_transcripts.fa.stats`, a TSV file with following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `GENE_ID`, the `gene_id` field in GTF.
  - `SEQNAME`, chromosome or contig name.
  - `START`, the `start` field in GTF, 1-based inclusive.
  - `END`, the `end` field in GTF, 1-based inclusive.
  - `STRAND`, the `strand` field in GTF.
  - `ABSOLUTE_LENGTH`, is `START` - `END` + 1.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `GC`, GC content of the cDNA in percentage.

Following is an example of `ce11_transcripts.fa.stats`:

```{code-cell}
cdna_stats = pd.read_table("ce11_trans_as.fa.stats", engine=CSV_ENGINE)
cdna_stats.head()
```

This can be used in multiple quality control scenarios, like histogram of GC content:

```{code-cell}
sns.histplot(cdna_stats, x="GC")
```

## Invocation of LLRGs

The Low-Level Read Generators (LLRGs) are programs that simulates DNA-Seq on some reference genome sequence by clipping reads in appropriate length and introducing sequencing errors. YASIM invokes LLRG on stranded cDNA sequences to generate RNA-Seq data. Here we would demonstrate their usage with following examples:

+++

### Invocation of NGS LLRGs: Use `art` for Example

`art_illumnina` is a general-purposed DNA-Seq simulator for Illumina NGS reads. It can be used to simulate error profiles for multiple Illumina sequencers and can perform pair-end (PE) sequencing.

The help message is as follows:

```{code-cell}
!python -m yasim art --help
```

Below is an example of using `art_illumina` to generate PE reads for Illumina HiSeq 2500 with default read length.

```{code-cell}
!if [ ! -f art_mode_1.fq ]; then \
    python -m yasim art \
        -F ce11_trans_as.fa.d \
        -j 40 \
        -d ce11_isoform_depth.tsv \
        -o art_mode \
        --sequencer_name HS25 \
        --pair_end_fragment_length_mean 300 \
        --pair_end_fragment_length_std 5 \
        --is_pair_end \
        2>&1 | grep -v WARNING; \
else \
    echo "art_mode_1.fq art_mode_2.fq already exists."; \
fi
```

Generates:

- `art_mode_1.fq` \& `art_mode_2.fq`, simulated PE FASTQ.
- `art_mode.d`, temporary directory that can be safely deleted.
- `art_mode.fq.stats`, statistics of simulated FASTQ. a TSV containing following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `INPUT_DEPTH`, isoform expression level in coverage provided by the upstream source.
  - `SIMULATED_N_OF_READS`, simulated number of reads. **This value can be used in assessing read quantifiers**.
  - `SIMULATED_N_OF_BASES`, simulated number of bases.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `SIMULATED_DEPTH`, simulated depth.

```{warning}
If the LLRG failed on some isoform, it would **NOT** appear in output FASTQ or its sattistics.
```

Example of `art_mode.fq.stats`:

```{code-cell}
art_mode_fq_stats = pd.read_table("art_mode.fq.stats", engine=CSV_ENGINE)
```

```{code-cell}
art_mode_fq_stats.head()
```

### Invocation of TGS LLRGs: Use `pbsim3` for Example

```{warning}
The official build of PBSIM, PBSIM2 and PBSIM3 shares a common executable anme (`pbsim`) but with different argument layout. For convenience, I renamed executable of PBSIM2 to `pbsim2` and PBSIM3 to `pbsim3`. If you do not use this in your computer, please use the `-e` option.
```

`pbsim3` is a general-purposed TGS DNA- and RNA-Seq simulator that supports multiple PacBio and Oxford Nanopore sequencers. It can generate Circular Consensus Sequence (CCS)/HiFi data.

Compared to NGS simulators, TGS simulators have `truncate_ratio_3p` and `truncate_ratio_5p`. These two parameters are used to set hard limits at two sides that allows simulation of incomplete reads due to reasons like 3' truncation.

The help message is as follows:

```{code-cell}
!python -m yasim pbsim3 --help
```

Following is an example of simulation of CCS data using PacBio RS II model:

```{code-cell}
!if [ ! -f pbsim3_mode.fq ]; then \
    python -m yasim pbsim3 \
        -F ce11_trans_as.fa.d \
        -j 40 \
        -e /home/yuzj/bin/pbsim3 \
        --ccs_pass 20 \
        -d ce11_isoform_depth.tsv \
        -m RSII \
        -M qshmm \
        --strategy trans \
        -o pbsim3_mode \
        2>&1 | grep -v WARNING; \
else \
    echo "pbsim3_mode.fq already exists."; \
fi
```

Generated files similar to NGS, so omitted.

+++

### Quality Control with Ease

Using utilities bundled in `labw_utils`, it is easy to perform quality control on reads you just generated without the need of installation of additional software! Look at following example:

```{code-cell}
!python -m labw_utils.bioutils describe_fastq pbsim3_mode.fq art_mode_1.fq
```

This generates:

- `pbsim3_mode.fq.stats.d` \& `art_mode_1.fq.stats.d`, the directory where quality control files are located. They are:
  - `all.tsv`, the per-read quality control file, with following columns:
    - `SEQID`, the FASTQ read ID.
    - `GC`, per-read GC content in absolute value.
    - `LEN`, actual read length.
    - `MEANQUAL`, mean sequencing quality using Phread33 score.
  - `extension_stat.tsv`, mean per-base quality of bases on each read from 5' to 3', mainly for NGS.
    - `POS`, position of each base on every transcript from 5' to 3'. If not present, would be omitted.
    - `QUAL`, mean sequencing quality using Phread33 score.

Example using `all.tsv`:

```{code-cell}
pbsim3_mode_all_qc = pd.read_table(os.path.join("pbsim3_mode.fq.stats.d", "all.tsv"), engine=CSV_ENGINE)
```

```{code-cell}
pbsim3_mode_all_qc.head()
```

For example, the read length distribution:

```{code-cell}
sns.histplot(pbsim3_mode_all_qc, x="LEN")
```

For example, the GC distribution:

```{code-cell}
sns.histplot(pbsim3_mode_all_qc, x="GC")
```

`extension_stat.tsv` may indicate whether clipping of terminal low-quality regions in NGS reads using [CutAdapt](https://cutadapt.readthedocs.io/en/stable) or [Trimmomatic](www.usadellab.org/cms/?page=trimmomatic) are required. For example:

```{code-cell}
art_mode_extension_qc = pd.read_table(os.path.join("art_mode_1.fq.stats.d", "extension_stat.tsv"), engine=CSV_ENGINE)
```

```{code-cell}
art_mode_extension_qc.head()
```

Following is a plot of mean per-base quality of all reads from 5' end to 3' end:

```{code-cell}
sns.lineplot(art_mode_extension_qc, x="POS", y="QUAL")
```

Which can be seen as suggestion of trimming leading and trailing bases. This is not done in this example due to time limit.

+++

## Salmon Quasi Alignment and Quantification

Salmon is used to precisely align and quantify reads mapped to transcriptome on an isoform-specific manner. Here we would make Salmon align to cDNA of ground truth GTF and see whether the quantification is accurate.

```{code-cell}
!if [ ! -d SALMON_IDX ]; then \
    salmon index \
        -i SALMON_IDX \
        -p 40 \
        -t ce11_trans_as.fa \
        &> salmon_index.log; \
else \
    echo "SALMON_IDX already exists."; \
fi
!salmon quant \
    -i SALMON_IDX \
    -l IU \
    -1 art_mode_1.fq -2 art_mode_2.fq \
    --validateMappings \
    -o art_mode_salmon \
    &> art_mode_salmon.log
!salmon quant \
    -i SALMON_IDX \
    -l U \
    -r pbsim3_mode.fq \
    --validateMappings \
    -o pbsim3_mode_salmon \
    &> pbsim3_mode_salmon.log
```

Compare them to ground truth.

```{code-cell}
def calculate_tpm(n_reads: pd.Series, transcribed_length: pd.Series) -> pd.Series:
    rpk = 1E3 * n_reads / transcribed_length
    return 1E6 * rpk / rpk.sum()

def read_salmon(file_path: str) -> pd.DataFrame:
    return (
    pd.read_table(file_path, comment="#")
    [["Name", "NumReads"]].
    rename(columns={"Name": "TRANSCRIPT_ID", "NumReads": "REAL_N_OF_READS"}).
    set_index('TRANSCRIPT_ID')
)

ngs_salmon_data = read_salmon(os.path.join("art_mode_salmon", "quant.sf"))
tgs_salmon_data = read_salmon(os.path.join("pbsim3_mode_salmon", "quant.sf"))

ngs_df = (
    pd.read_table("art_mode.fq.stats", engine=CSV_ENGINE).
    join(ngs_salmon_data, on="TRANSCRIPT_ID").
    fillna(0)
    [["TRANSCRIPT_ID", "SIMULATED_N_OF_READS", "REAL_N_OF_READS", "TRANSCRIBED_LENGTH"]]
)
tgs_df = (
    pd.read_table("pbsim3_mode.fq.stats", engine=CSV_ENGINE).
    join(tgs_salmon_data, on="TRANSCRIPT_ID").
    fillna(0)
    [["TRANSCRIPT_ID", "SIMULATED_N_OF_READS", "REAL_N_OF_READS", "TRANSCRIBED_LENGTH"]]
)
ngs_df["GENERATION"] = "NGS"
tgs_df["GENERATION"] = "TGS"

# Calculate TPM for NGS and TGS data.
for df in (ngs_df, tgs_df):
    df["TPM_SIM"] = calculate_tpm(
        df["SIMULATED_N_OF_READS"],
        df["TRANSCRIBED_LENGTH"]
    )
    df["TPM_REAL"] = calculate_tpm(
        df["REAL_N_OF_READS"],
        df["TRANSCRIBED_LENGTH"]
    )
# Merge and filter low-expression isoforms.
ngs_tgs_merged_df = pd.concat((ngs_df, tgs_df)).query(
    "TPM_SIM > 10 & TPM_REAL > 10 & TPM_SIM < 3000 & TPM_REAL < 3000"
)
# Calculate Log 2 Fold Change (L2FC).
ngs_tgs_merged_df["TPM_L2FC"] = np.log2(
    ngs_tgs_merged_df["TPM_SIM"] / ngs_tgs_merged_df["TPM_REAL"]
)
```

Plotting simulated (`TPM_SIM`) vs. actual (`TPM_REAL`) data.

```{code-cell}
# Getting axis limits.
xylim = max(ngs_tgs_merged_df["TPM_SIM"].max(), ngs_tgs_merged_df["TPM_REAL"].max())
g = sns.FacetGrid(
    ngs_tgs_merged_df,
    col="GENERATION",
    xlim=(0, xylim),
    ylim=(0, xylim),
    height=5,
    aspect=1
)
g.map(sns.scatterplot, "TPM_SIM", "TPM_REAL", alpha=0.4)
```

Plotting log 2 Fold Change (`TPM_L2FC`) vs. base mean (`TPM_SIM`) data.

```{code-cell}
g = sns.FacetGrid(
    ngs_tgs_merged_df,
    col="GENERATION",
    xlim=(-4, 4),
    ylim=(0, xylim),
    height=5,
    aspect=1
)
g.map(sns.scatterplot, "TPM_L2FC",  "TPM_SIM" , alpha=0.4)
```

From above plot, it is evident that in more complex genomes, TGS data could outperform NGS ones.

+++

List of files:

```{code-cell}
!ls -lFh
```

```{code-cell}

```

---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Additional Tutorial on LLRGs

Although most LLRGs are straightforward, there are still some LLRGs that takes a long time to understand. This is a tutorial to correctly use different LLRGs to simulate different sequencers. We assume that you had already read _YASIM Tutorial_.

Here we would demonstrate features of commonly-used LLRGs using mitochondria of _C. Elegans_.

+++

## Introduction

Low-Level Read Generators (LLRGs) are programs used to generate machine errors. Normally, a LLRG can be represented in 2 forms: As a standalone third-party executable (e.g., PBSIM can benevoked by `pbsim` command) or as a Python module that can be imported (e.g., dTGS). These two forms are unified using LLRG adapter, a builtin middleware inside YASIM that performs execution and error handling of LLRGs. The LLRG adapters are evoked by LLRG frontend (i.e., the `python -m yasim [LLRG]` command) in bulk or single-cell RNA-Seq experiments. Following is a workflow of LLRG step:

```{figure} ../fig/llrg_step.svg
:width: 100%
:align: left
:alt: LLRG Step Flow Chart

LLRG Step Flow Chart

This figure demonstrates the basic workflow of the LLRG step. It is, in details, as follows:

1. The input cDNA sequences and depths are separated by isoform's transcript ID.
2. For each transcript, its cDNA sequence and depth are passed into LLRG adapter, which performs DNA-Seq on those cDNAs and generates reads in FASTQ format.
3. The LLRG adapters are executed in parallel.
4. Generated reads would be merged into one (SE) or two (PE) files.
```

So after knowing these knowledge, it would be easy to decode LLRG parameter specifications. Among all possible parameters, the most important are `--simulator_name` and `-e`. The former specifies what would appear in FASTQ SEQID, and the latter specifies path to LLRG executable. For example, if you installed PBSIM3 in `/home/yuzj/bin/pbsim3` instead of normal locations, you should specify `-e /home/yuzj/bin/pbsim3`.

Compared to NGS simulators, TGS simulators supports two additional parameter: `--truncate_ratio_3p` and `--truncate_ratio_5p`. These two parameters specifies 3' and 5' truncation where 3' and 5' are **IN RESPECT TO SEQUENCER**. For example, suppose that a isoform is defined on forward strand. PBSIM1 would take its cDNA (on forward strand) and generate reads in **both** forward strand and reverse strand. At this time, if we specified a 3' truncation, in respect to reference genome, the forward cDNA-Seq read would be clipped on 3' end while reverse cDNA-Seq strand would be clipped on 5' end.

The argument parser of LLRG frontend supports pass-through. That is, arguments or options that cannot be recognized would be applied to all LLRG. For example, if we wish to adjust the accuracy of PBSIM, we can add `--accuracy-mean [VALUE]` parameter at the end of `python -m yasim pbsim ...` command. Since `--accuracy-mean` is not a recognizable parameter of YASIM, it would be passed to all `pbsim` process.

+++

## Preparations

Here generates all steps before proceeding into LLRG. It would download CE11 reference genome sequence and annotation from UCSC and generate sequencing depth using V3API. Note that generation of AS events is optional so not used.

```{code-cell}
:tags: [skip-execution]

%%bash
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz
tar xzvf chromFa.tar.gz
rm -f chrI.fa chrII.fa chrIII.fa chrIV.fa

axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
gzip -cfd ce11.ncbiRefSeq.gtf.gz | grep '^chrM\s' >chrM.ncbiRefSeq.gtf
python -m yasim generate_gene_depth \
    -g chrM.ncbiRefSeq.gtf \
    -o gene_depth.tsv \
    -d 60
python -m yasim generate_isoform_depth \
    -g chrM.ncbiRefSeq.gtf \
    -d gene_depth.tsv \
    -o isoform_depth.tsv \
    --alpha 4
python -m yasim generate_gene_depth \
    -g chrM.ncbiRefSeq.gtf \
    -o gene_low_depth.tsv \
    -d 5
python -m yasim generate_isoform_depth \
    -g chrM.ncbiRefSeq.gtf \
    -d gene_low_depth.tsv \
    -o isoform_low_depth.tsv \
    --alpha 4
```

To interactively seen LLRG simulation statistics, following Python modules are imported:

```{code-cell}
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
```

## ART

ART is a general-purposed NGS DNA-Seq simulator. YASIM uses ART's Illumina mode (Executable `art_illumina`) with version information as follows:

```{code-cell}
%%bash
art_illumina | head -n 6
```

### Specifying Sequencer Model and Read Length

The output of `python -m yasim art --help` may be hard to read, and we would provide an explanation here. Although ART simulator can simulate Illumina sequencers with different model and read length, the choises of read length is limited by the chose of models. For example, if we choose Illumina GenomeAnalyzer I (`GA1` as is in parameter), the valid read length would be 36 or 44. If the read length is not valid or not specified, the default read length (the first valid read length) would be used.

Following is a table of supported ART sequencers and read length:

```{code-cell}
:tags: [hide-input]

from yasim.llrg_adapter.art import AVAILABLE_ILLUMINA_ART_SEQUENCER
art_compatibility_matrix = pd.DataFrame(AVAILABLE_ILLUMINA_ART_SEQUENCER).transpose()
art_compatibility_matrix.columns = ("Real Name", "Allowed Read lengths")
art_compatibility_matrix
```

So, for example, to simulate single-end (SE) HiSeq 2500 with 150 read length, you need to:

```shell
python -m yasim art \
    -F [FASTAS DIR] \
    -d [ISOFORM DEPTH] \
    -o [OUTPUT_NAME] \
    --sequencer_name HS25 \
    --read_length 150
```

+++

### Supporting Pair-End (PE) Mode

ART also supports pair-end (PE) simulation. Under that circumstance, one additional flag, `--is_pair_end`, needs to be set and two additional parameters, `--pair_end_fragment_length_mean` and `--pair_end_fragment_length_std`, needs to be specified.

### ART Errors

Art may raise following errors:

```text
OpenBLAS blas thread_init: pthread_create failed for thread 14 of 128: Resource temporarily unavailable
OpenBLAS blas thread init: RLIMIT_NPROC 4096 current, 2061441 max
```

Solve this error by setting `-j` parameter to a smaller value or wait until the server becomes less busy.

+++

## DWGSIM

DWGSIM is a NGS simulator. YASIM uses its illumina PE mode.

+++

## Badread

Badread is a TGS simulator implemented in Python. Compared to other simulator, it is relatively slow.

+++

## PBSIM

PBSIM version 1 can simulate data generated by PacBio RS, a very old sequencer that is not commonly seen nowadays.

PBSIM has a CCS mode. You can use this mode by specifying `--ccs` flag.

```{warning}
PBSIM generates CCS with a mechanism similiar to CLR. It does **NOT** require invocation of PacBio `ccs` utility.
```

+++

## PBSIM2

PBSIM2 can simulate early PacBio RS II and ONT models. It does **NOT** support CCS generation.

Using PBSIM2 for TGS simulation is recommended. It is fast, reliable and accurate.

+++

## PBSIM3

PBSIM3 is one of the most complicated LLRGs used in this software. Comparing to PBSIM3, it has following additional parameters:

### Strategy

PBSIM3 have 2 strategy: `wgs` and `trans`. Their difference is as follows:

- The `wgs` strategy is same as what is used in PBSIM2, PBSIM, and other LLRGs. i.e., PBSIM3 was used as a DNA-Seq simulator that reads cDNA and outputs sequences.
- The `trans` strategy is new in PBSIM3. It would generate reads based on a new PBSIM3-specific isoform record format.

An example is as follows:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy wgs \
    -F chrm_trans.fa.d \
    -d isoform_depth.tsv \
    -o chrm_pbsim3_wgs \
    -j 40
python -m yasim pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy trans \
    -F chrm_trans.fa.d \
    -d isoform_depth.tsv \
    -o chrm_pbsim3_trans \
    -j 40
python -m labw_utils.bioutils describe_fastq chrm_pbsim3_trans.fq chrm_pbsim3_wgs.fq
python -m labw_utils.bioutils describe_gtf chrM.ncbiRefSeq.gtf
```

```{code-cell}
pbsim3_trans_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_trans.fq.stats.d", "all.tsv"),
    quotechar="'"
)
pbsim3_wgs_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_wgs.fq.stats.d", "all.tsv"),
    quotechar="'"
)
pbsim3_trans_fq_stats["TRANSCRIPT_ID"] = pbsim3_trans_fq_stats["SEQID"].apply(lambda s: s.split(":")[0])
pbsim3_wgs_fq_stats["TRANSCRIPT_ID"] = pbsim3_wgs_fq_stats["SEQID"].apply(lambda s: s.split(":")[0])

ref_genome_stats = pd.read_table("chrM.ncbiRefSeq.gtf.transcripts.tsv").set_index("TRANSCRIPT_ID")

pbsim3_trans_fq_stats = pbsim3_trans_fq_stats.join(ref_genome_stats, on="TRANSCRIPT_ID")
pbsim3_wgs_fq_stats = pbsim3_wgs_fq_stats.join(ref_genome_stats, on="TRANSCRIPT_ID")
```

Following is a plot of read length on FASTQ vs. transcribed length on GTF.

```{code-cell}
fig, axs = plt.subplots(nrows=2, sharex=True)
sns.boxplot(pbsim3_trans_fq_stats, y="LEN", x="TRANSCRIBED_LENGTH", ax=axs[0])
sns.boxplot(pbsim3_wgs_fq_stats, y="LEN", x="TRANSCRIBED_LENGTH", ax=axs[1])
axs[0].set_title("trans mode")
axs[1].set_title("wgs mode")
```

```{warning}
The `trans` mode in YASIM can only simulate cDNA of corresponding strand. It would **NOT** generate reads on the reverse strand.
```

+++

### HMM Method

PBSIM3 have 2 HMM method: Error HMM and Quality Score HMM. Please refer to the paper for more details. A apparent difference is that the former would not generate meaningful quality scores in output FASTQ (all zero). See following example:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim pbsim3 \
    -m RSII \
    -M errhmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy trans \
    -F chrm_trans.fa.d \
    -d isoform_low_depth.tsv \
    -o chrm_pbsim3_errhmm \
    -j 40
python -m yasim pbsim3 \
    -m RSII \
    -M qshmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy trans \
    -F chrm_trans.fa.d \
    -d isoform_low_depth.tsv \
    -o chrm_pbsim3_qshmm \
    -j 40
python -m labw_utils.bioutils describe_fastq chrm_pbsim3_errhmm.fq chrm_pbsim3_qshmm.fq
```

```{code-cell}
pbsim3_errhmm_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_errhmm.fq.stats.d", "all.tsv"),
    quotechar="'"
)
pbsim3_qshmm_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_qshmm.fq.stats.d", "all.tsv"),
    quotechar="'"
)
```

Plotting of mean quality score per read.

```{code-cell}
fig, axs = plt.subplots(nrows=2, sharex=True)
sns.histplot(pbsim3_errhmm_fq_stats, x="MEANQUAL", ax=axs[0])
sns.histplot(pbsim3_qshmm_fq_stats, x="MEANQUAL", ax=axs[1], binwidth=1)
axs[0].set_title("errhmm mode")
axs[1].set_title("qshmm mode")
```

The QSHMM does not support PacBio Sequel model since they do not have sequence quality information.

+++

### CCS

Circular Consensus Sequence (CCS)/HiFi Reads are commonly used in PacBio sequencing as it provides user with higher accuracy that is comparable to NGS reads. See [official site](https://ccs.how/) for more details.

YASIM can generate CCS FASTQ. The `--ccs_pass` parameter determines number of passes to perform when simulating CCS reads. To simulate CLR reads, do not set `--ccs_pass` or set `--ccs_pass` to 1. If this parameter is set to other value, simulated data would be in CCS.

On CCS generation, YASIM would firstly invoke PBSIM3 to generate PacBio CLR reads, and then call CCS using `ccs` utility (which is slow) from PacBio. The MAF generated in CCS mode is paired with generated CLR reads and cannot reflect error status of generated CCS reads.

CCS generation requires `samtools` and `ccs` to be present. You may set their path in corresponding parameters. For IsoSeq-based pipeline that requires CCS BAM, please refer to the appendix. 3' and 5' truncation set in YASIM parameters are applicable for CCS FASTQ but not applicable to CCS BAM.

See following example for details:

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim pbsim3 \
    -m RSII \
    -M qshmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy trans \
    -F chrm_trans.fa.d \
    -d isoform_low_depth.tsv \
    -o chrm_pbsim3_clr \
    --ccs_pass 1 \
    -j 40
fi
python -m yasim pbsim3 \
    -m RSII \
    -M qshmm \
    -e /home/yuzj/bin/pbsim3 \
    --strategy trans \
    -F chrm_trans.fa.d \
    -d isoform_low_depth.tsv \
    -o chrm_pbsim3_ccs \
    --ccs_pass 10 \
    -j 40
fi
python -m labw_utils.bioutils describe_fastq chrm_pbsim3_clr.fq chrm_pbsim3_ccs.fq
```

```{code-cell}
pbsim3_ccs_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_ccs.fq.stats.d", "all.tsv"),
    quotechar="'"
)
pbsim3_clr_fq_stats = pd.read_table(
    os.path.join("chrm_pbsim3_clr.fq.stats.d", "all.tsv"),
    quotechar="'"
)
```

Plotting of mean quality score per read.

```{code-cell}
fig, axs = plt.subplots(nrows=2, sharex=True)
sns.histplot(pbsim3_ccs_fq_stats, x="MEANQUAL", ax=axs[0], binwidth=1)
sns.histplot(pbsim3_clr_fq_stats, x="MEANQUAL", ax=axs[1], binwidth=1)
axs[0].set_title("CCS mode")
axs[1].set_title("CLR mode")
```

## dTGS

dTGS (Dumb TGS) is a TGS simulator that outputs all given contig with given depth times without any error (Phread33 score `K`, the highest one). It is ultra-fast and accurate but without biological meanings.

+++

## Appendices

+++

### Using PacBio IsoSeq Pipelines

```{warning}
In this mode, `truncate_ratio_5p` and `truncate_ratio_3p` cannot be effective.
```

CCS reads generated by `pbsim3` can be used in officially supported PacBio [IsoSeq](https://isoseq.how) pipelines. To finish this tutorial, you need to install PacBio [SMRTLink](https://www.pacb.com/support/software-downloads/) or its [community version](https://github.com/PacificBiosciences/pbbioconda) (recommended). Version of Dependencies:

```{code-cell}
%%bash
pbmerge --version | head -n 1
pbindex --version | head -n 1
samtools --version | head -n 1
ccs --version | head -n 1
pbmm2 --version | head -n 1
isoseq3 --version | head -n 1
```

Generation of CCS reads. We would use PacBio Sequel for example.

```{code-cell}
:tags: [skip-execution]

%%bash
python -m yasim pbsim3 \
    -m SEQUEL \
    -M errhmm \
    -e /home/yuzj/bin/pbsim3 \
    -F chrm_trans.fa.d \
    -d isoform_low_depth.tsv \
    -o chrm_ccs_isoseq \
    -j 40 \
    --ccs_pass 10
```

Merge all small CCS BAMs into single CCS BAM. The file `merge.py` is provided as follows:

```{literalinclude} merge.py
:language: python
```

```{code-cell}
:tags: [skip-execution]

%%bash
python merge.py chrm_ccs_isoseq.ccs.bam chrm_ccs_isoseq.d/*/tmp.ccs.bam
pbindex chrm_ccs_isoseq.ccs.bam
samtools index chrm_ccs_isoseq.ccs.bam
```

Then you can use standard PacBio IsoSeq pipeline. For example:

```{code-cell}
:tags: [skip-execution]

%%bash
isoseq3 cluster \
    chrm_ccs_isoseq.ccs.bam \
    chrm_ccs_isoseq.transcripts.xml \
    --log-level INFO \
    --num-threads 40
pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --log-level INFO \
    chrm_ccs_isoseq.transcripts.xml.hq.bam \
    chrM.fa \
    chrm_ccs_isoseq.aln.bam
isoseq3 collapse \
    --do-not-collapse-extra-5exons \
    --log-level INFO \
    chrm_ccs_isoseq.aln.bam \
    chrm_ccs_isoseq.ccs.bam \
    chrm_ccs_isoseq.collapse.gff
```

The generated annotation file would be available at `chrm_ccs_isoseq.collapse.gff`. You are free to use `gffcompare` or `pigeon` for further analysis.

+++

### Interpretation of LLRG Exceptions

After each simulation, the LLRG adapter would print a line like this: `2023-04-03 15:16:41,070	[INFO] Status of errors: {'NORMAL': 5, 'LLRGFail': 6}`. This line indicates LLRG exception status. Below the definition of commonly-seen exceptions:

- `NORMAL`: If no exception occurs
- `EmptyOutFile`: If LLRG exited normally but with empty output file.
- `NoOutputFile`: If LLRG exited normally but with no output file.
- `LLRGFail`: If LLRG exited abnormally.
- `InitFail`: if pre-execution of initialization hook failed.
- `UNKNOWN`: Other errors.

```{figure} ../fig/llrg_adapter.svg
:width: 100%
:align: left
:alt: LLRG Adapter
```

+++

### Hint on Management of LLRGs

You may use wrapper scripts for LLRGs that requires complex prerequisites.

Following is a wrapper for Badread. This script would:

1. Search for `badread` executable. If succeeded, would execute that executable.
2. Search for `badread` Conda environment. If succeeded, would activate that environment and use `badread` executable inside.
3. Setup `badread` Conda environment and use `badread` executable inside.

```shell
#!/usr/bin/env bash
set -e
if which badread &>> /dev/null; then
    exec badread "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!" >&2
    exit 127
fi

if ! conda env list | grep ^badread &>> /dev/null; then
    conda create -y -n badread -c bioconda badread=0.2.0 python-edlib
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate badread
exec badread "${@}"
```

# `yasim` -- Yet Another SIMulator

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

YASIM is a read simulator for Next- and Third-Generation bulk RNA Sequencing with Alternative Splicing and realistic Gene Expression Profile. It can be used to benchmark various of tools that claimed to be able to detect or quantify isoforms. Here we would provide detailed guidance on usage of YASIM.

This documentation is a small instruction for users of YASIM.

## Installation

### Using pre-built Library from PYPA

You need Python interpreter (CPython implementation) >= 3.7 and latest [`pip`](https://pip.pypa.io/) to install this software from PYPA. Command:

```shell
pip install yasim==1.0.0
```

You are recommended to use this application inside a virtual environment like [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io) or [`poetry`](https://python-poetry.org).

### Build from Source

You need Python interpreter (CPython implementation) >= 3.7, latest PYPA [`build`](https://pypa-build.readthedocs.io) and latest [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build the simulator using:

```shell
python3 -m build
pip install dist/yasim-1.0.0.tar.gz
```

## YASIM Quickstart

### Third-Party Files

This program relies on **matched** reference **genomic** GTF and FASTA as input. You may get them from authentic authorities. Following is a table where you can download Reference GTF and FASTA:

| Name    | Genome FASTA                                                                                                                                                                  | Genome GTF                                                                                                                                                                    | cDNA FASTA                                                                                                      |
|---------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| UCSC    | [hg38.p13](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)                                                                                                | [ncbiRefSeq](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)                                                                            | NA                                                                                                              |
| Ensembl | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)                                                    | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz)                                                                          | [GRCh38.105](http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) |
| NCBI    | [GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz) | [GRCh38.p14](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gtf.gz) | NA                                                                                                              |

```{warning}
This simulator does not support fancy genomic features like decoy sequences, HLA sequences, EBV sequences, pathes or alternative loci. You are free to use references with those features but generated data may not be biologically meaningful.
```

### List All Available Subcommands

Coexistence of different functionality in YASIM works through subcommands. Through `python3 -m yasim lscmd`, you would see some output like:

```text
2023-01-19 04:20:27,739 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 0.2.13
2023-01-19 04:20:27,739 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py lscmd
2023-01-19 04:20:27,739 [INFO] Listing modules...
badread
dwgsim
generate_as_events
generate_depth
pbsim
pbsim2
transcribe
```

Inside the output, you would find 3 lines of log and 7 subcommands sorted in alphabetical order. We would introduce all subcommands in logical order below.

### The Example

Inside the example, _C. Elegans_ reference data from UCSC is used. See:

- <https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz>
- <https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz>

To reduce time consumption, only chromosome 1 was used.

```{warning}
This version of YASIM uses `labw_utils` 0.1.X GTF parser. This parser is NOT stable and is to be removed. You shold use **UNSORTED** UCSC references for compatibility.
```

### Generate Alternative Splicing Events: `generate_as_events`

This step would generate alternative splicing events. It would take reference GTF (as those downloaded from UCSC, etc.) as input and generate GTF with AS events as output.

The generated GTF should be seen as ground truth for benchmarking AS detectors. If you wish to benchmark quantifiers only, this step can be safely omitted.

Synopsis: `python3 -m yasim generate_as_events [-h] -f [FASTA] -g [GTF] -o [OUT]`

Required arguments:

- `-f [FASTA]`, `--fasta [FASTA]` Reference genome, in FASTA format
- `-g [GTF]`, `--gtf [GTF]` Reference genome, in GTF format
- `-o [OUT]`, `--out [OUT]` Output GTF

Optional arguments:

- `-h`, `--help` show this help message and exit

Example:

```shell
python -m yasim generate_as_events -f ce11.fa -g ce11.ncbiRefSeq_sorted.gtf -o ce11.ncbiRefSeq_as.gtf
```

```text
2023-01-19 04:16:07,405 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 0.2.13
2023-01-19 04:16:07,405 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py generate_as_events -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq_as.gtf
Reading ce11.ncbiRefSeq.gtf:   0%|                                             | 0/74066 [00:00<?, ?L/s]
2023-01-19 04:16:07,737 [WARNING] Gene Y74C9A.6 is inferred from feature transcript
2023-01-19 04:16:07,738 [WARNING] Gene homt-1 is inferred from feature transcript
2
[...]
Reading ce11.ncbiRefSeq.gtf: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 74066/74066 [00:01<00:00, 53647.41L/s]
2023-01-19 04:16:09,118 [INFO] Pickling to gvpkl...
2023-01-19 04:16:10,827 [INFO] Loaded 3961 genes with 5300 transcript
Generating isoforms...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 3961/3961 [00:10<00:00, 380.15it/s]
2023-01-19 04:16:21,274 [INFO] Will remove 1135 genes out of 3961
2023-01-19 04:16:21,341 [INFO] Will remove genes FIN
Standardizing transcripts: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 3841/3841 [00:00<00:00, 89509.25it/s]
Standardizing genes: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 3922/3922 [00:00<00:00, 85025.10it/s]
```

Generates:

- `ce11.ncbiRefSeq_as.gtf`.
- `ce11.ncbiRefSeq.gtf.0.4.gvpkl.xz`, if not exist. This is a cache file for the `labw_utils` GTF parser.

### Generate Sequencing Depth of Gene: `generate_depth`

This step would generate Gene Expression Profile (GEP) over one GTF. The input can be reference GTF or GTF from previous step.

Synopsis: `python3 -m yasim generate_depth [-h] -f [FASTA] -g [GTF] -o [OUT]`

Required arguments:

- `-g [GTF]`, `--gtf [GTF]` Reference genome, in GTF format
- `-o [OUT]`, `--out [OUT]` Output TSV
- `-d [MU]`, `--mu [MU]` Average depth

Optional arguments:

- `-h`, `--help` show this help message and exit

```{warning}
The generated coverage is **NOT** number of reads generated! It cannot be used as ground truth to assess quantification software! The number of reads ground truth will be provided by LLRG UIs introduced below.
```

Example:

```shell
python -m yasim generate_depth -g ce11.ncbiRefSeq_as.gtf -o ce11_depth.tsv -d 5
```

```text
2023-01-19 04:40:16,317 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 0.2.13
2023-01-19 04:40:16,317 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py generate_depth -g ce11.ncbiRefSeq_as.gtf -o ce11_depth.tsv -d 5
Reading ce11.ncbiRefSeq_as.gtf.0.4.gvpkl.xz: 100%|XXXXXXXXXXXXXXXXX| 5.97M/5.97M [00:00<00:00, 24.6MB/s]
Simulating...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 3841/3841 [00:00<00:00, 411355.37it/s]
```

Generates:

- `ce11_depth.tsv`, a TSV file with following columns:
  - `TRANSCRIPT_ID`
  - `DEPTH`
- `ce11.ncbiRefSeq_as.gtf.0.4.gvpkl.xz`, if not present.

### Transcribe GTF to FASTA: `transcribe`

This step would transcribe the input genome GTF and genome FASTA into stranded transcriptome FASTA.

This step is designed to be general-purposed. It can be applied on any matching GTF and FASTA.

This step should generate similiar output with `bedtools getfasta -nameOnly -s -fi [FASTA] -bed [GTF] > [OUT]`

Synopsis: `python3 -m yasim transcribe -f [FASTA] -g [GTF] -o [OUT]`

Optional arguments:

- `-h`, `--help` show this help message and exit

Required arguments:

- `-f [FASTA]`, `--fasta [FASTA]` Reference genome, in FASTA format
- `-g [GTF]`, `--gtf [GTF]` GTF to be transcribed
- `-o [OUT]`, `--out [OUT]` Output FASTA. Except from `[OUT]`, this step would also generate `[OUT].d`, which is a directory of FASTAs with only one transcript inside.

```{note}
Although this software can be used to generate reference cDNAs for software like Salmon, there are differences between transcribed cDNA and Ensembl-provided cDNA. Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.
```

Example:

```shell
python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq_as.gtf -o ce11_transcripts.fa
```

```text
2023-01-19 04:26:07,014 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 0.2.13
2023-01-19 04:26:07,014 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py transcribe -f ce11.fa -g ce11.ncbiRefSeq_as.gtf -o ce11_transcripts.fa
Reading ce11.ncbiRefSeq_as.gtf.0.4.gvpkl.xz: 100%|XXXXXXXXXXXXXXXXX| 5.97M/5.97M [00:00<00:00, 39.9MB/s]
Reading ce11.fa: 301451L [00:00, 881722.09L/s]
Transcribing GTF...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 3841/3841 [00:01<00:00, 2877.34it/s]
```

Generates:

- `ce11_transcripts.fa`.
- `ce11_transcripts.fa.d`, the directory where every transcipt is stored as separate FASTA.
- `ce11_transcripts.fa.stats`, a TSV file with following columns:
  - `TRANSCRIPT_ID`
  - `GENE_ID`
  - `SEQNAME`
  - `START`
  - `END`
  - `STRAND`
  - `ABSOLUTE_LENGTH`, is `START` - `END`
  - `TRANSCRIBED_LENGTH`, length of the CDNA
  - `GC`

### Simulate: Use Badread for Example

This step would add machine noise and abundance information to the transcribed FASTA.

The machine noises are add by something called Low-Level Read Generators (LLRGs). They are independent third-party programs that experts in simulating sequencer noises but does not provide pre-transcriptional modification of abundance simulation functionalities. You should have them installed in advance.

Currently supported LLRGs are [Badread](https://github.com/rrwick/Badread), [PBSIM](https://github.com/pfaucon/PBSIM-PacBio-Simulator), [PBSIM2](https://github.com/yukiteruono/pbsim2) and [DWGSIM](https://github.com/nh13/DWGSIM).

Synopsis: `python3 -m yasim badread -F [FASTAS] -g [GTF] -o [OUT] -d [DEPTH] -e [EXENAME] -j [JOBS]`

Required arguments:

- `-F [FASTAS]`, `--fastas [FASTAS]` Directory of transcribed DGE FASTAs from `transcribe` step
- `-g [GTF]`, `--gtf [GTF]` GTF to be transcribed
- `-o [OUT]`, `--out [OUT]` Output transcript prefix
- `-d [DEPTH]`, `--depth [DEPTH]` Depth generated by `generate_depth` step
- `-m [MODEL_NAME]`, `--model_name [MODEL_NAME]`. Used model. Should be one of  -m `nanopore2018`,`nanopore2020`,`pacbio2016`,`verybad`,`verynice`

Optional arguments:

- `-h`, `--help` show this help message and exit
- `-e [EXENAME]`, `--exename [EXENAME]` Executable name or absolute path.
- `-j [JOBS]`, `--jobs [JOBS]` Number of threads used.

```{warning}
The official build of PBSIM and PBSIM2 shares a common executable anme (`pbsim`) but with different argument layout. For convenience, I renamed executable of PBSIM2 to `pbsim2`. If you do not use this in your computer, please use the `-e` option.
```

````{hint}
You may use wrapper scripts for LLRGs that requires complex prerequisites.

Following is a wrapper for Badread. This script would:

1. Search for `badread` executable. If succeed, would execute that executable.
2. Search for `badread` Conda environment. If succeed, would activate that environment and use `badread` executable inside.
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

Similiar script for PBSIM:

```shell
#!/usr/bin/env bash
set -e
if which pbsim &>> /dev/null; then
    exec pbsim "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!" >&2
    exit 127
fi

if ! conda env list | grep ^pbsim &>> /dev/null; then
    conda create -y -n pbsim -c bioconda pbsim=1.0.3
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pbsim
exec pbsim "${@}"
```
````

Example:

```shell
python3 -m yasim pbsim2 \
    -F ce11.fa.d \
    -o hg38_pbsim \
    -d ce11_depth.tsv \
    --hmm_model P4C2 \
    --exename /home/yuzj/bin/pbsim2 \
    -j 40
```

```text
2023-01-19 05:04:29,970 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 0.2.13
2023-01-19 05:04:29,971 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py pbsim2 -F ce11.fa.d -o hg38_pbsim -d ce11_depth.tsv --hmm_model P4C2 --exename /home/yuzj/bin/pbsim2 -j 40
Reading depth file...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 716/716 [00:00<00:00, 1861823.72it/s]
Submitting jobs...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 716/716 [00:00<00:00, 88949.76it/s]
Simulating jobs: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 716/716 [00:16<00:00, 43.38it/s]
Merging...: 100%|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 716/716 [00:00<00:00, 81860.16it/s]
```

## NEWS

- 2022-11-09: The `yasim.commonutils` and `yasim.bioutils` package would be separated into a new package (still under internal development) and `yasim` of later versions would depend on it.

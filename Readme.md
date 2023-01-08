# `yasim` -- Yet Another SIMulator

Read simulator for Next- and Third-Generation Sequencing with Alternative Splicing and Differently Expressed Genes.

## Installation

### Using pre-built Library from PYPA

You need Python interpreter (CPython implementation) >= 3.7 and latest [`pip`](https://pip.pypa.io/) to install this software from PYPA. Command:

```shell
pip install yasim==1.0.0
```

You are recommended to use this application inside a virtual environment.

### Build from Source

You need Python interpreter (CPython implementation) >= 3.7 and latest [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build the simulator using:

```shell
python3 setup.py sdist
pip install dist/yasim-1.0.0.tar.gz
```

If you wish to contribute in development of the package, see Conda environment at `env/yasim_dev.yml` for tools that builds documentation.

## `yasim` Quickstart

List all available subcommands: `python3 -m yasim lscmd`. You would see some output like:

```text
2023-01-09 03:28:38,268 [INFO] yasim -- Yet Another SIMulator for Alternative Splicing and Differentially Expressed Gene ver. 1.0.0
2023-01-09 03:28:38,268 [INFO] Called by: /home/yuzj/Documents/yasim/src/yasim/__main__.py lscmd
2023-01-09 03:28:38,268 [INFO] Listing modules...
badread
dwgsim
generate_as_events
generate_depth
pbsim
pbsim2
transcribe
```

### Generate Alternative Splicing Events: `generate_as_events`

usage: `python3 -m yasim generate_as_events [-h] -f [FASTA] -g [GTF] -o [OUT]`

optional arguments:

- `-h`, `--help` show this help message and exit
- `-f [FASTA]`, `--fasta [FASTA]` Reference genome, in FASTA format
- `-g [GTF]`, `--gtf [GTF]` Reference genome, in GTF format
- `-o [OUT]`, `--out [OUT]` Output GTF

### Generate Sequencing Depth of Gene: `generate_depth`

usage: `python3 -m yasim generate_depth [-h] -f [FASTA] -g [GTF] -o [OUT]`

Arguments:

- `-h`, `--help` show this help message and exit
- `-g [GTF]`, `--gtf [GTF]` Reference genome, in GTF format
- `-o [OUT]`, `--out [OUT]` Output TSV
- `-d [MU]`, `--mu [MU]` Average depth

### Transcribe GTF to FASTA: `transcribe`

- `-h`, `--help` show this help message and exit
- `-f [FASTA]`, `--fasta [FASTA]` Reference genome, in FASTA format
- `-g [GTF]`, `--gtf [GTF]` GTF to be transcribed
- `-o [OUT]`, `--out [OUT]` Output FASTA

### Simulate: Use `badread` for Example

- `-h`, `--help` show this help message and exit
- `-F [FASTAS]`, `--fastas [FASTAS]` Directory of transcribed DGE FASTAs from `transcribe` step
- `-g [GTF]`, `--gtf [GTF]` GTF to be transcribed
- `-o [OUT]`, `--out [OUT]` Output transcript prefix
- `-d [DEPTH]`, `--depth [DEPTH]` Depth generated by `generate_depth` step
- `-m [MODEL_NAME]`, `--model_name [MODEL_NAME]`. Used model. Should be one of  -m `nanopore2018`,`nanopore2020`,`pacbio2016`,`verybad`,`verynice`
- `-e [EXENAME]`, `--exename [EXENAME]` Executable name or absolute path.
- `-j [JOBS]`, `--jobs [JOBS]` Number of threads used.

### Put them Together

This is an example of simulating AS events using hg38 genome with PBSIM2 as LLRG.

```shell
python3 -m yasim generate_as_events -f hg38.fa -g hg38.gtf -o hg38_as.gtf
python3 -m yasim generate_depth -g hg38_as.gtf -o hg38_depth.tsv -d 100
python3 -m yasim transcribe -f hg38.fa -g hg38_as.gtf -o hg38_trans.fa
python3 -m yasim pbsim2 \
    -F hg38_trans.fa.d \
    -o hg38_pbsim \
    -d hg38_depth.tsv \
    --hmm_model P4C2 \
    --exename PATH_TO_PBSIM2 \
    -j 40
```

## Hints

To simulate using simulators like `pbsim`, `pbsim2`, etc., you are recommended to provide the **ABSOLUTE PATH** to the executable using `-e` argument.

````{hint}
You may use wrapper scripts for LLRGs that requires a large environment to run.

Provided wrapper scripts:

Following is a wrapper for `badread`. This script would:

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
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^badread &>> /dev/null; then
    conda create -y -n badread -c bioconda badread==0.2.0 python-edlib
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate badread
exec badread "${@}"
```

Similiar script for `pbsim`:

```shell
#!/usr/bin/env bash
set -e
if which pbsim &>> /dev/null; then
    exec pbsim "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^pbsim &>> /dev/null; then
    conda create -y -n pbsim -c bioconda pbsim==1.0.3
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pbsim
exec pbsim "${@}"
```

Similiar script for `pbsim2`. Mind that in official distribution, the executable name of `pbsim` and `pbsim2` is the same. So this wrapper would search for an executable named `pbsim2`.

```shell
#!/usr/bin/env bash
set -e
if which pbsim2 &>> /dev/null; then
    exec pbsim2 "${@}"
fi
if ! which conda &>> /dev/null; then
    echo "conda not found!"
    exit 127
fi

if ! conda env list | grep ^pbsim2 &>> /dev/null; then
    conda create -y -n pbsim2 -c bioconda pbsim2==2.0.1
fi

eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pbsim2
exec pbsim "${@}"
```

````

## NEWS

- 2022-11-09: The `yasim.commonutils` and `yasim.bioutils` package would be separated into a new package (still under internal development) and `yasim` of later versions would depend on it.

## Steps

### Generating Differentially Expressed Genes and Alternative Splicing



### Call Low-Level Simulators (aka. Sequence Generators)

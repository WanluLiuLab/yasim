# YASIM -- Yet Another SIMulator for Alternative Splicing Events and Realistic Gene Expression Profile

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

With the development of Third-Generation Sequencing (TGS) and related technologies, accurate quantification of transcripts in the isoform level with precise detection of novel isoforms from Alternative Splicing (AS) events or relocation of Transposable Elements (TEs) had become possible. YASIM is the tool that simulates Next- or Third-Generation bulk RNA-Sequencing raw FASTQ reads with ground truth genome annotation and realistic gene expression profile (GEP). It can be used to benchmark tools that are claimed to be able to detect isoforms (e.g., [StringTie](https://ccb.jhu.edu/software/stringtie/)) or quantify reads on an isoform level (e.g., [featureCounts](https://subread.sourceforge.net/featureCounts.html)).

YASIM serves different simulation purposes. For example, it can be used to simulate a count matrix from reference genome annotation or to simulate raw FASTQ reads from a user-provided count matrix. When combined with other tools, the user can also simulate reads from the genome with Single Nucleotide Polymorphism (SNP), Insertions \& Deletions (InDels), Structural Variations (SVs), and other genomic variations.

YASIM cannot generate machine noises for each sequencer, and third-party DNA- or RNA-Seq simulators (Referred to as Low-Level Read generators, LLRGs) are needed to convert cDNA sequences to reads with machine errors and quality information (Except PacBio Sequel/Sequel II model). This gives YASIM extreme flexibility over sequencer models. Till now, YASIM can simulate most Illumina NGS sequencers and most PacBio/ONT TGS sequencing platforms.

YASIM is designed to be modularized, as some of the modules are general-purpose and can be used in other simulation tasks. Implemented in Python 3, YASIM follows Object-Oriented Programming (OOP) styles and can be easily extended. Theoretically, YASIM can run on any platform that supports Python3. However, most LLRG are POSIX-only (i.e., work on GNU/Linux, MacOS, and friends). So it is recommended to deploy this tool inside major GNU/Linux distributions like Ubuntu, Debian, CentOS, Fedora, etc. Using YASIM on Microsoft Windows Subsystem of Linux (WSL), version 1 or 2, is **NOT** recommended -- It would lead to impaired performance and may cause other problems due to LLRG incompatibilities. Using YASIM on other platforms (e.g., Oracle Solaris) is neither tested nor recommended.

## Installation

### Using the Pre-Built Version from PYPI

You need a working [Python](https://www.python.org) interpreter (CPython implementation) >= 3.7 (**recommended 3.8**) and the latest [`pip`](https://pip.pypa.io/) to install this software from [PYPI](https://pypi.org). Command:

```shell
pip install yasim==3.1.6
```

You are recommended to use this application inside a virtual environment like [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io), or [`poetry`](https://python-poetry.org).

### Build from Source

Before building from the source, get a copy of the latest source code from <https://github.com/WanluLiuLab/yasim> using [Git](https://git-scm.com):

```shell
git clone https://github.com/WanluLiuLab/yasim
```

Or, if you prefer to use [GNU Wget](https://www.gnu.org/software/wget).

```shell
wget -o yasim-master.zip https://github.com/WanluLiuLab/yasim/archive/refs/heads/master.zip
unzip yasim-master.zip
```

You need Python interpreter (CPython implementation) >= 3.7, latest PYPA [`build`](https://pypa-build.readthedocs.io), and [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build and install the simulator using:

```shell
cd yasim
python3 -m build
pip install dist/yasim-3.1.6-py3-none-any.whl
```

### Installation of Third-Party Programs

For TGS LLRGs:

- [PBSIM](https://github.com/yukiteruono/pbsim), which simulates PacBio RS C1 and C2 chemistry, with CCS support.
- [PBSIM2](https://github.com/yukiteruono/pbsim2), which simulates PacBio RS II P4C2, P5C3, and P6C4 chemistry, CLR only; ONT R9.4, R9.5, and R10.3 pore.
- [PBSIM3](https://github.com/yukiteruono/pbsim3), which simulates PacBio RS II and Sequel model, with CCS support.
- [BadRead](https://github.com/rrwick/Badread), which simulates arbitrary PacBio and ONT models.

For NGS LLRGs:

- [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm), which simulates various Illumina NGS platforms. Namely, GenomeAnalyzer I, GenomeAnalyzer II, HiSeq 1000, HiSeq 2000, HiSeq 2500, HiSeqX PCR free, HiSeqX TruSeq, MiniSeq TruSeq, MiSeq v3, NextSeq500 v2.
- [DWGSIM](https://github.com/nh13/DWGSIM), which simulates arbitrary Illumina models.

You may refer to the LLRG tutorial for detailed guidance on the utilization of these pieces of software.

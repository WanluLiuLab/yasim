# `yasim` -- Yet Another SIMulator

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

YASIM is a read simulator for Next- and Third-Generation bulk RNA Sequencing with Alternative Splicing and realistic Gene Expression Profile. It can be used to benchmark various of tools that claimed to be able to detect or quantify isoforms.

## Installation

### Using pre-built Library from PYPA

You need Python interpreter (CPython implementation) >= 3.7 (recommended 3.8) and latest [`pip`](https://pip.pypa.io/) to install this software from PYPA. Command:

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

### Installation of Third-Party Programs

YASIM does not have the ability to generate machine noises for each sequencer, and third-party DNA- or RNA-Seq simulators (Referred to as Low-Level Read generators, or LLRGs) are needed to convert cDNA sequences to reads with machine errors and quality information [^qual]. A list of supported LLRGs:

For TGS LLRGs:

- `pbsim`, which simulates PacBio RS C1 and C2 chemistry, with CCS support.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/pbsim) [Debian](https://packages.debian.org/stable/pbsim)
  - Y. Ono, K. Asai, and M. Hamada, "Pbsim: Pacbio reads simulator–toward accurate genome assembly.," _Bioinformatics (Oxford, England)_, vol. 29, pp. 119–121, 1 Jan. 2013, ISSN : 1367-4811. DOI: [10.1093/bioinformatics/bts649](https://doi.org/10.1093/bioinformatics/bts649)
- `pbsim2`, which simulates PacBio RS II P4C2, P5C3 and P6C4 chemistry, CLR only; ONT R94, R95 and R103 chemistry.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim2)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/pbsim2)
  - Y. Ono, K. Asai, and M. Hamada, "PBSIM2: A simulator for long-read sequencers with a novel generative model of quality scores," _Bioinformatics (Oxford, England)_, vol. 37, no. 5, pp. 589–595, May 5, 2021, Number: 5, ISSN: 1367-4811. DOI: [10.1093/bioinformatics/btaa835](https://doi.org/10.1093/bioinformatics/btaa835)
- `pbsim3`, which simulates PacBio RS II and Sequel model, with CCS support.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim3)
  -  Y. Ono, M. Hamada, and K. Asai, "Pbsim3: A simulator for all types of pacbio and ont long reads.," _NAR genomics and bioinformatics_, vol. 4, lqac092, 4 Dec. 2022, ISSN: 2631-9268. DOI: [10.1093/nargab/lqac092](https://doi.org/10.1093/nargab/lqac092)
- `badread`, which simulates arbitrary PacBio and ONT models.
  - Official site: [GitHub](https://github.com/rrwick/Badread)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/badread)
  - R. Wick, "Badread: Simulation of error-prone long reads," Journal of Open Source Software, vol. 4, no. 36, p. 1316, Apr. 2019. DOI: [10.21105/joss.01316](https://doi.org/10.21105/joss.01316)

For NGS LLRGs:

- `art`, which simulates Illumina GenomeAnalyzer I, GenomeAnalyzer II, HiSeq 1000, HiSeq 2000, HiSeq 2500, HiSeqX PCR free, HiSeqX TruSeq, MiniSeq TruSeq, MiSeq v3, NextSeq500 v2.
  - Official site: [NIEHS](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/art) [Debian](https://packages.debian.org/stable/art-nextgen-simulation-tools)
  - W. Huang, L. Li, J. R. Myers, and G. T. Marth, "Art: A next-generation sequencing read simulator.," _Bioinformatics (Oxford, England)_, vol. 28, pp. 593–594, 4 Feb. 2012, ISSN: 1367-4811. DOI: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708)
- `dwgsim`, which simulates arbitrary Illumina models.
  - Official site: [GitHub](https://github.com/nh13/DWGSIM)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/dwgsim) [Debian](https://packages.debian.org/stable/dwgsim)
  - NO PUB

You may refer to LLRG tutorial for detailed guidance on utilization of these software.

[^qual]: Except PacBio Sequel model.

## NEWS

- 2022-11-09: The `yasim.commonutils` and `yasim.bioutils` package would be separated into a new package (still under internal development) and `yasim` of later versions would depend on it.

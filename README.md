[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![GitHub release](https://img.shields.io/github/release/ngs-fzb/mtbseq_source.svg)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mtbseq/README.html)
![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/mtbseq.svg)
![Conda](https://img.shields.io/conda/dn/bioconda/mtbseq.svg)

[![Docker](https://quay.io/repository/biocontainers/mtbseq/status)](https://quay.io/repository/biocontainers/mtbseq/)

# MTBseq

MTBseq is an automated pipeline for mapping, variant calling and detection of resistance mediating and phylogenetic variants from Illumina whole genome sequence data of *Mycobacterium tuberculosis* complex isolates.

## Getting Started

For complete installation instructions, description and usage examples please read the [MANUAL.md](https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md).

# Installation

## Conda
Install [Conda](https://conda.io/docs/) or [Miniconda](https://conda.io/miniconda.html) hereafter install MTBseq with:
```
conda install -c bioconda mtbseq
```

## Source
Please see the [MANUAL.md](https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md) for installation from source.

### Requirements

```
* Perl: Perl 5, version 22, subversion 1 (v5.22.1)
* Java: Oracle Java 8 or OpenJDK 8 (no other version work with the GenomeAnalysisTK 3.8)

** MTBseq uses the following CPAN and core modules: **
* MCE                 (v1.833)
* Statistics::Basic   (v1.6611)

* FindBin             (v1.51)core
* Cwd                 (v3.62)core
* Getopt::Long        (v2.5)core
* File::Copy          (v2.30)core
* List::Util          (v1.49)core
* Exporter            (v5.72)core
* vars                (v1.03)core
* lib                 (v0.63)core
* strict              (v1.09)core
* warnings            (v1.34)core

** MTBseq uses the following third party software: **
** Binaries (compiled on Ubuntu 16.04) are included, except for GenomeAnalysisTK 3.8**
* bwa                 (v0.7.17)
* GenomeAnalysisTK    (v3.8)
* picard              (v2.17.0)
* samtools            (v1.6)
```

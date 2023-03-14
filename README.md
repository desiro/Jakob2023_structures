# [<samp>Jakob2023 structures</samp>](https://github.com/desiro/Jakob2023_structures)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

RNA structure prediction script for <samp>Jakob2023</samp>.

### Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![ViennaRNA v2.5.0](https://img.shields.io/badge/ViennaRNA_v2.5.0-006795.svg)](https://www.tbi.univie.ac.at/RNA/)

### Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Installation

I recommend using Miniconda and following the steps below. If this is the first time using conda, you should probably restart your shell after the installation of Miniconda. The following will demonstrate the installation and set up of Miniconda on Linux, which should be similar on other platforms. For Windows 10 users, I advise using the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem. More information can be found on the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [Bioconda](https://bioconda.github.io/) pages.

### Conda Installation

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Updating Miniconda and setting channels:
```
conda update conda
conda update python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
conda create --name Jakob2023_structures python=3.9.7
conda activate Jakob2023_structures
conda install -c bioconda viennarna=2.5.0
git clone https://github.com/desiro/Jakob2023_structures.git
cd Jakob2023_structures
```

***

## Execution

```
python SPLASHPeaksStructures.py
```

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite <samp>Jakob2023</samp>.

```
C. Jakob, G.L. Lovate, D. Desirò, L. Gießler, R. Smyth, R. Marquet, K. Lamkiewicz, M. Marz, M. Schwemmle and H. Bolte.
"Sequential disruption of SPLASH-identified vRNA-vRNA interactions challenges their role in influenza A virus genome packaging."
In review, 2023.
```

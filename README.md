# VecTrap

> **Don't let synthetic vector contamination slip into your pathogen assemblies. Catch them with VecTrap.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()

---

## Overview

VecTrap is a high-throughput Python pipeline for the detection and classification of synthetic laboratory vector sequences within raw bacterial genome assemblies and plasmid datasets.

VecTrap uses a catalog-driven homology scanning strategy backed by an empirically derived sequence database built from 120,684 real Addgene plasmids. Detection covers replication origins, resistance markers, promoters, terminators, primer binding sites, recombination sites, and other synthetic vector elements.

---

## Architecture

```
vectrap/
    catalogs/          FASTA sequence catalogs + search indexes
    modules/           Core scanning and scoring modules
db/
    build_db.py        One-time catalog index preparation
vectrap.py             Unified CLI entry point
requirements.txt
```

---

## Installation

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
pip install -r requirements.txt
```

---

## Citation

If you use VecTrap in your research, please cite:

> *VecTrap: A catalog-driven pipeline for detecting synthetic vector contamination in bacterial genome assemblies.* (Manuscript in preparation)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

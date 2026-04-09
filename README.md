# VecTrap

> **Don't let synthetic vector contamination slip into your pathogen assemblies. Catch them with VecTrap.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()

---

## Overview

VecTrap is a high-throughput, computationally lightweight Python pipeline for the detection and classification of synthetic laboratory vector sequences within raw bacterial genome assemblies and plasmid datasets.

Synthetic vector contamination in public databases (NCBI RefSeq/GenBank) is a well-documented and pervasive problem. Researchers often clone fragments of pathogen genomes into pUC-, pET-, or BioBrick-based vectors, then accidentally deposit the complete vector backbone sequence alongside the biological target sequence. VecTrap is designed to systematically flag these contaminants before they can corrupt downstream genomic, phylogenetic, or epidemiological analyses.

VecTrap uses a multi-module, evidence-tiered scanning strategy to scan both DNA strands in a single pass, avoiding expensive protein translation or BLAST alignment for its primary candidate selection stage.

---

## Key Features

- **Blazing-fast regex-based scanning** of raw nucleotide sequences (no 6-frame translation required)
- **Evidence-tiered classification** — every hit is classified as `PRIMARY` (standalone proof) or `SUPPORTIVE` (requires genomic context)
- **Context-Aware Filtering** — short assembly markers like BioBrick scars are only flagged when their required flanking elements are also present on the same contig
- **Bidirectional strand scanning** with precise 0-based coordinate mapping back to the forward sense strand
- **Modular architecture** — each marker class is implemented as an independent, swappable Python module
- **Pure Python 3** — zero external dependencies for all candidate selection modules

---

## Modules

| Module | Description | Dependencies |
| :--- | :--- | :--- |
| `origin_scanner_module.py` | Detects synthetic replication origins (ColE1, p15A, RK2, pBBR1, f1, SV40, 2μ) | None |
| `regulatory_scanner_module.py` | Detects engineered promoters (tac/trc/T5/T7/SP6/T3), sequencing primers (M13/VF2/VR), BioBrick assembly elements (RFC10), and terminators | None |
| `mcs_module.py` | Identifies Multiple Cloning Site clusters by restriction site density peak-finding | None |
| `peptide_tag_module.py` | Detects common fusion protein tags (His, FLAG, Strep-II, SUMO, MBP, GST, thrombin/TEV cleavage sites) | None |

---

## Installation

VecTrap requires Python 3.8 or higher. All candidate selection modules are written in pure Python and require no external dependencies.

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
```

For the planned homology-based validation modules (Reporter Proteins, Selectable Markers), additional dependencies will be required:

```bash
pip install -r requirements.txt
```

---

## Quick Start

Run each module independently against a FASTA or gzipped FASTA file:

```bash
# Detect synthetic replication origins
python origin_scanner_module.py -i your_plasmids.fasta -o results/

# Detect engineered promoters, primers, and BioBrick elements
python regulatory_scanner_module.py -i your_plasmids.fasta -o results/

# Detect Multiple Cloning Sites
python mcs_module.py -i your_plasmids.fasta -o results/

# Detect peptide fusion tags
python peptide_tag_module.py -i your_plasmids.fasta -o results/
```

---

## Output Format

All modules write a tab-separated values (TSV) file to the specified output directory. The regulatory scanner output includes the following columns:

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header of the originating contig |
| `start_0based` | 0-based start coordinate of the hit on the forward strand |
| `end_0based` | 0-based end coordinate of the hit on the forward strand |
| `strand` | `+` for forward, `-` for reverse complement |
| `marker_name` | Identifier of the detected regulatory element |
| `evidence_class` | `PRIMARY` or `SUPPORTIVE` |
| `architecture` | `monolithic` or `bipartite_laco` |
| `matched_sequence` | Exact nucleotide sequence matched by the regex engine |

---

## Evidence Classification

VecTrap partitions every hit into one of two evidence tiers:

- **`PRIMARY`** — The matched sequence is a fully synthetic, laboratory-designed element (e.g., a tac hybrid promoter, M13 sequencing primer, RFC10 BioBrick prefix). Finding even a single `PRIMARY` hit on a contig is sufficient evidence to classify it as an engineered vector.
- **`SUPPORTIVE`** — The matched sequence can also occur in wild-type phages or natural bacterial genomes (e.g., a T7 promoter core, lacWT promoter, lambda pL). These hits require additional genomic context to confirm engineering and are used as supporting evidence only.

---

## Citation

If you use VecTrap in your research, please cite:

> *VecTrap: A high-throughput pipeline for detecting synthetic vector contamination in bacterial genome assemblies.* (Manuscript in preparation)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

# VecTrap


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()

---

## Overview

VecTrap is a high-throughput Python pipeline for the detection and classification of synthetic laboratory vector sequences within raw bacterial genome assemblies and plasmid datasets.

Synthetic vector contamination in public databases (NCBI RefSeq/GenBank) is a well-documented and pervasive problem. VecTrap is designed to systematically flag these contaminants before they can corrupt downstream genomic, phylogenetic, or epidemiological analyses.

VecTrap uses a catalog-driven homology scanning strategy backed by an empirically derived sequence database built from more than 100,000 experimentally validated synthetic constructs. The catalog covers the full diversity of engineered vector elements including replication origins, selectable markers, regulatory sequences, recombination sites, and sequencing primer binding sites found in prokaryotic, broad-host-range, and mammalian expression vectors.

---

## Key Features

- **Catalog-driven homology detection** against a curated database of real synthetic vector sequences
- **Evidence-tiered classification** -- every hit is classified as `PRIMARY` (standalone proof) or `SUPPORTIVE` (requires genomic context)
- **Broad vector coverage** -- prokaryotic, broad-host-range, lentiviral, and mammalian expression vector backbones
- **Bidirectional strand scanning** with precise coordinate mapping back to the forward sense strand
- **Two output levels** -- per-hit evidence table and per-contig classification (VECTOR / SUSPECTED / CLEAN)
- **Modular architecture** -- scanning and scoring are fully decoupled

---

## Detected Vector Elements

| Category | Evidence Tier | Examples |
| :--- | :--- | :--- |
| Replication origins | PRIMARY | ColE1, p15A, RK2, pBBR1, f1, SV40, 2-micron |
| Transfer origins | PRIMARY | oriT IncP, IncQ, IncW |
| Recombination sites | PRIMARY | attB/attP, loxP, FRT |
| Ribosome binding sites | PRIMARY | Shine-Dalgarno variants, synthetic RBS |
| Selectable markers | SUPPORTIVE | KanR, AmpR, CmR, HygR, PuroR, BlastR |
| Engineered promoters | SUPPORTIVE | T7, T5, tac, trc, CMV, EF1a, PGK, CAG |
| Terminators | SUPPORTIVE | T7Te, rrnB, BGH, SV40 polyA |
| Primer binding sites | SUPPORTIVE | M13, T7, SP6, pUC, T3 sequencing primers |
| Protein binding sites | SUPPORTIVE | lacO, tetO, cumate operator |
| Retroviral elements | SUPPORTIVE | LTR sequences, lentiviral packaging signals |
| Regulatory elements | SUPPORTIVE | IRES, Kozak, splice donor/acceptor |
| Mobile elements | SUPPORTIVE | Transposon remnants in vector backbones |

---

## Architecture

```
vectrap/
    catalogs/      FASTA sequence catalogs and search indexes
    modules/       Core scanning and scoring modules
db/
    build_db.py    One-time catalog index preparation
vectrap.py         Unified CLI entry point
requirements.txt
```

---

## Installation

VecTrap requires Python 3.8 or higher.

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
pip install -r requirements.txt
```

---

## Quick Start

```bash
# Prepare catalog indexes (run once after cloning)
python db/build_db.py --catalog-dir vectrap/catalogs/

# Run full pipeline
python vectrap.py -i assembly.fasta -o results/ -c vectrap/catalogs/
```

---

## Output Format

VecTrap writes two output files per run to the specified output directory.

**Per-hit evidence table** (`*.all_hits.tsv`):

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header of the originating contig |
| `start_0based` | 0-based start coordinate on the forward strand |
| `end_0based` | 0-based end coordinate on the forward strand |
| `strand` | `+` for forward, `-` for reverse complement |
| `evidence_class` | `PRIMARY` or `SUPPORTIVE` |
| `marker` | Catalog category of the detected element |
| `detail` | Matched catalog entry identifier |

**Per-contig verdict table** (`*.verdicts.tsv`):

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header |
| `classification` | `VECTOR`, `SUSPECTED`, or `CLEAN` |
| `primary_hits` | Number of PRIMARY-tier hits |
| `supportive_score` | Weighted sum of SUPPORTIVE-tier hits |
| `evidence_summary` | Human-readable evidence breakdown |

---

## Evidence Classification

VecTrap partitions every hit into one of two evidence tiers:

- **`PRIMARY`** -- The matched sequence is an unambiguously synthetic, laboratory-designed element (replication origin, transfer origin, recombination site, or synthetic RBS). A single `PRIMARY` hit on a contig is sufficient to classify it as an engineered vector.
- **`SUPPORTIVE`** -- The matched sequence is strongly associated with synthetic vectors but may also occur in natural contexts (promoter cores, terminator sequences, primer binding sites). These hits contribute to a weighted evidence score and require corroboration to confirm engineering.

---

## Citation

If you use VecTrap in your research, please cite:

> *VecTrap: A catalog-driven pipeline for detecting synthetic vector contamination in bacterial genome assemblies.* (Manuscript in preparation)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

# VecTrap

> **Don't let synthetic vector contamination slip into your pathogen assemblies. Catch them with VecTrap.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844271.svg)](https://doi.org/10.5281/zenodo.20844271)

---

## Overview

VecTrap is a high-throughput Python pipeline for the detection and classification of synthetic laboratory vector sequences within raw bacterial genome assemblies and plasmid datasets.

Synthetic vector contamination in public databases (NCBI RefSeq/GenBank) is a well-documented and pervasive problem. Researchers often clone fragments of pathogen genomes into engineered cloning and expression vectors, then accidentally deposit the complete vector backbone sequence alongside the biological target sequence. VecTrap is designed to systematically flag these contaminants before they can corrupt downstream genomic, phylogenetic, or epidemiological analyses.

VecTrap uses a catalog-driven homology scanning strategy backed by an empirically derived sequence database built from more than 100,000 experimentally validated synthetic constructs. The catalog covers the full diversity of engineered vector elements including replication origins, selectable markers, regulatory sequences, recombination sites, and sequencing primer binding sites found in prokaryotic, broad-host-range, and mammalian expression vectors.

---

## Sequence Catalog

The VecTrap sequence catalog is deposited on Zenodo and must be downloaded separately before running the pipeline.

**DOI: [10.5281/zenodo.20844271](https://doi.org/10.5281/zenodo.20844271)**

The catalog consists of 14 FASTA files, one per feature type category, plus a manifest file. Download and index preparation is handled automatically by the `vectrap-build-db` command (see Quick Start below).

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
    __init__.py
    cli/
        build_db.py    vectrap-build-db entry point
        run.py         vectrap entry point
    modules/
        utils.py
        homology_scanner.py
        scorer.py
    catalogs/
        README.md      how to obtain catalog files
pyproject.toml
requirements.txt
```

---

## Installation

VecTrap requires Python 3.8 or higher and [minimap2](https://github.com/lh3/minimap2) available in `PATH`.

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
pip install .
```

Or once published on PyPI:

```bash
pip install vectrap
```

---

## Quick Start

```bash
# Download catalog files from Zenodo and prepare indexes (run once after installation)
vectrap-build-db --download

# Or use locally provided catalog files
vectrap-build-db --catalog-dir /path/to/catalogs/

# Run the pipeline
vectrap -i assembly.fasta -o results/ -c vectrap/catalogs/
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

For the sequence catalog:

> VecTrap Sequence Catalog v1.0. Zenodo. https://doi.org/10.5281/zenodo.20844271

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

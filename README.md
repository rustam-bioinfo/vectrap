# VecTrap

> **Identify sequences of synthetic origin in bacterial genome assemblies and plasmid datasets.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844271.svg)](https://doi.org/10.5281/zenodo.20844271)

---

## Overview

VecTrap is a high-throughput Python pipeline for detecting and classifying sequences of synthetic origin in bacterial genome assemblies and plasmid datasets. Its primary use case is identifying assembled contigs that derive fully or partially from engineered laboratory constructs — cloning vectors, expression vectors, lentiviral transfer plasmids, and similar — rather than from the biological sample under study.

Synthetic sequence contamination in public databases (NCBI RefSeq/GenBank) is a well-documented problem. Researchers clone fragments of pathogen or environmental genomes into engineered vector backbones, then accidentally deposit the complete construct alongside the intended biological sequence. VecTrap is designed to flag these cases systematically before they corrupt downstream genomic, phylogenetic, or epidemiological analyses.

VecTrap uses a catalog-driven homology scanning strategy backed by a curated sequence database of more than 100,000 experimentally validated synthetic construct elements. The catalog covers the full diversity of engineered vector components found in prokaryotic, broad-host-range, and mammalian expression systems.

---

## Sequence Catalog

The VecTrap sequence catalog is deposited on Zenodo and must be downloaded separately before running the pipeline.

**DOI: [10.5281/zenodo.20844271](https://doi.org/10.5281/zenodo.20844271)**

The catalog consists of 14 FASTA files, one per feature type category, plus a manifest file. Download and index preparation is handled automatically by the `vectrap-build-db` command (see Quick Start below).

---

## Key Features

- **Catalog-driven homology detection** against a curated database of real synthetic vector sequences
- **Evidence-tiered classification** — every hit is assigned to `ENGINEERED` (no natural analog; standalone proof of synthetic origin) or `RECRUITED` (natural sequence recruited into vectors; interpreted in genomic context)
- **Broad vector coverage** — prokaryotic, broad-host-range, lentiviral, and mammalian expression vector backbones
- **Bidirectional strand scanning** with precise coordinate mapping back to the forward sense strand
- **Two output levels** — per-hit evidence table and per-contig classification (`VECTOR` / `SUSPECTED` / `CLEAN`)
- **Modular architecture** — scanning and scoring are fully decoupled

---

## Detected Vector Elements

| Category | Evidence Tier | Examples |
| :--- | :--- | :--- |
| Synthetic ribosome binding sites | ENGINEERED | BBa-series RBS, optimised Shine-Dalgarno variants |
| Site-specific recombination sites | ENGINEERED | attB/attP (Lambda, PhiC31), loxP, FRT |
| Fully synthetic promoters | ENGINEERED | T7, T5, tac, trc, CMV, EF1α, PGK, CAG |
| Synthetic terminators | ENGINEERED | T7Te, BGH polyA, SV40 polyA |
| Sequencing primer binding sites | ENGINEERED | M13, T7, SP6, pUC, T3 |
| Retroviral LTR sequences | ENGINEERED | LTR sequences, lentiviral packaging signals |
| Synthetic regulatory elements | ENGINEERED | IRES, synthetic Kozak, splice donor/acceptor |
| Replication origins | RECRUITED | ColE1, p15A, RK2, pBBR1, f1 (natural plasmid/phage ancestors) |
| Transfer origins | RECRUITED | oriT IncP, IncQ, IncW |
| Selectable markers | RECRUITED | KanR, AmpR, CmR, HygR, PuroR, BlastR (Tn-derived) |
| Natural-sequence promoters in vectors | RECRUITED | lac, ara, trp, phoA, tet |
| Natural terminators in vectors | RECRUITED | rrnB T1/T2, lambda t0 |
| Protein binding sites | RECRUITED | lacO, tetO, cumate operator |
| Mobile element remnants | RECRUITED | Tn3, Tn5, Tn10, Tn903, IS borders |

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
        README.md      catalog acquisition and tier reference
pyproject.toml
requirements.txt
```

---

## Installation

VecTrap requires Python 3.8 or higher.

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
| `evidence_tier` | `ENGINEERED` or `RECRUITED` |
| `feature_class` | Functional class of the detected element |
| `catalog_id` | Matched catalog entry identifier |

**Per-contig verdict table** (`*.verdicts.tsv`):

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header |
| `classification` | `VECTOR`, `SUSPECTED`, or `CLEAN` |
| `engineered_hits` | Number of ENGINEERED-tier hits |
| `recruited_classes` | Number of distinct functional classes with RECRUITED hits |
| `evidence_summary` | Human-readable evidence breakdown |

---

## Evidence Classification

VecTrap assigns every catalog hit to one of two evidence tiers:

- **`ENGINEERED`** — The matched sequence has no natural biological analog. It was designed in a laboratory and does not exist in wild-type organisms. A single high-confidence ENGINEERED hit is sufficient to classify a contig as synthetic. Examples: T7 promoter, CMV enhancer, BGH polyA signal, attB/FRT recombination sites, synthetic RBS, sequencing primer binding sites, retroviral LTR sequences in non-viral assemblies.

- **`RECRUITED`** — The matched sequence derives from a natural biological source that was recruited into synthetic vectors because it functions well in engineered contexts. These sequences also occur in wild genomes and cannot individually prove synthetic origin. Their evidential weight comes from **functional class diversity**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter is almost certainly synthetic, whereas a contig with only one of these is ambiguous. Examples: ColE1 origin, KanR/AmpR resistance genes, rrnB terminator, conjugative oriT, lac promoter.

---

## Repository Docs

- [`INFO.md`](INFO.md) — technical internals: algorithms, data structures, coordinate conventions, evidence model, scoring logic, performance notes
- [`vectrap/catalogs/README.md`](vectrap/catalogs/README.md) — full catalog tier reference with per-category rationale

---

## Citation

If you use VecTrap in your research, please cite:

> *VecTrap: A catalog-driven pipeline for detecting sequences of synthetic origin in bacterial genome assemblies.* (Manuscript in preparation)

For the sequence catalog:

> VecTrap Sequence Catalog v1.0. Zenodo. https://doi.org/10.5281/zenodo.20844271

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

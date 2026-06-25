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
- **Evidence-tiered classification** — every hit is assigned to `ENGINEERED` (no natural analog in wild-type bacteria; standalone proof of synthetic origin), `CONTEXT_DEPENDENT` (natural in bacteria but anomalous depending on host organism and genomic context), or `WEAK` (widespread in wild-type bacteria; contributes only in combination with other evidence)
- **Broad vector coverage** — prokaryotic, broad-host-range, lentiviral, and mammalian expression vector backbones
- **Bidirectional strand scanning** with precise coordinate mapping back to the forward sense strand
- **Two output levels** — per-hit evidence table and per-contig classification (`VECTOR` / `SUSPECTED` / `CLEAN`)
- **Modular architecture** — scanning and scoring are fully decoupled

---

## Detected Vector Elements

All tier assignments reflect the bacterial genome context — the question asked is whether a given sequence could plausibly appear in a wild-type bacterial chromosome or natural bacterial plasmid.

| Category | Evidence Tier | Examples |
| :--- | :--- | :--- |
| Eukaryotic/viral replication origins | ENGINEERED | SV40 ori, EBV oriP |
| Filamentous phage origins used as cloning tools | ENGINEERED | f1 ori, M13 ori |
| Eukaryotic and viral promoters | ENGINEERED | CMV, EF1α, PGK, CAG, CaMV35S |
| Phage RNA polymerase promoters | ENGINEERED | T7, SP6, T3 |
| Eukaryotic transcription terminators / polyA signals | ENGINEERED | BGH polyA, SV40 polyA |
| Mammalian selectable markers | ENGINEERED | Puromycin, Blasticidin, Hygromycin (mammalian context) |
| Fluorescent and reporter proteins | ENGINEERED | EGFP, mCherry, mTurquoise, Luciferase |
| CRISPR/genome-editing elements | ENGINEERED | Cas9, Cas12, sgRNA scaffolds |
| Recombinase systems | ENGINEERED | Cre, Flp, PhiC31 |
| Site-specific recombination sites | ENGINEERED | attB/attP (Lambda, PhiC31), loxP, FRT |
| Synthetic regulatory elements | ENGINEERED | IRES, WPRE, P2A/T2A self-cleaving peptides |
| Sequencing primer binding sites | ENGINEERED | M13, T7, SP6, pUC, T3 |
| Retroviral/lentiviral elements | ENGINEERED | LTR, lentiviral packaging signals, HIV PSI |
| Bacterial replication origins | CONTEXT_DEPENDENT | ColE1, pMB1, pUC ori, p15A, pSC101, RK2, pBBR1 |
| Bacterial transfer origins | CONTEXT_DEPENDENT | oriT IncP, IncQ, IncW |
| Bacterial selectable markers | CONTEXT_DEPENDENT | KanR, AmpR, CmR, TetR, GenR (Tn-derived, endemic via HGT) |
| T7 expression system elements | CONTEXT_DEPENDENT | T7 terminator, T7 tag (T7 RNAP not endogenous in wild-type bacteria) |
| Inducible expression system components | CONTEXT_DEPENDENT | lacI, lacO, IPTG-responsive elements, tac/trc promoters |
| Natural terminators used in vectors | WEAK | rrnB T1/T2, lambda t0, tonB terminator |
| Natural sigma-factor promoters | WEAK | lac, ara, trp, phoA (widespread in enterobacteria) |
| Mobile element remnants | CONTEXT_DEPENDENT | Tn3, Tn5, Tn10, Tn903, IS borders |

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
| `evidence_tier` | `ENGINEERED`, `CONTEXT_DEPENDENT`, or `WEAK` |
| `feature_class` | Functional class of the detected element |
| `catalog_id` | Matched catalog entry identifier |

**Per-contig verdict table** (`*.verdicts.tsv`):

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header |
| `classification` | `VECTOR`, `SUSPECTED`, or `CLEAN` |
| `engineered_hits` | Number of ENGINEERED-tier hits |
| `recruited_classes` | Number of distinct functional classes with CONTEXT_DEPENDENT hits |
| `evidence_summary` | Human-readable evidence breakdown |

---

## Evidence Classification

VecTrap is designed for **bacterial genome assemblies**. All tier assignments reflect whether a sequence could plausibly occur in a wild-type bacterial chromosome or natural bacterial plasmid.

Every catalog hit is assigned to one of three evidence tiers:

- **`ENGINEERED`** — The matched sequence has no natural occurrence in any wild-type bacterium. It was designed in a laboratory, is derived from eukaryotic viruses or mammalian systems, or originates from contexts (phage RNA polymerase, mammalian cell culture) that have no plausible route into a bacterial genome without deliberate cloning. A single high-confidence ENGINEERED hit is sufficient to classify a contig as synthetic. Examples: T7 promoter, CMV enhancer, SV40 ori, f1/M13 ori, BGH polyA signal, EGFP, Cas9, attB/FRT recombination sites, IRES, WPRE, retroviral LTR.

- **`CONTEXT_DEPENDENT`** — The matched sequence exists naturally in bacteria or natural bacterial plasmids but its presence in an assembly may indicate synthetic contamination depending on the host organism and genomic context. A ColE1 origin in an *E. coli* assembly is unremarkable; the same origin in a *Salmonella* or *Bacillus* assembly warrants scrutiny. Antibiotic resistance genes (KanR, AmpR, CmR) are endemic in environmental and clinical bacteria via horizontal gene transfer and are not diagnostic alone. Their evidential weight comes from **functional class co-occurrence**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter is almost certainly synthetic. Examples: ColE1 ori, p15A ori, KanR, AmpR, CmR, tac/trc promoters, T7 terminator, lacI/lacO.

- **`WEAK`** — The matched sequence is widespread across wild-type bacteria with low specificity for synthetic origin. These elements contribute to a verdict only when co-occurring with stronger evidence from other tiers. Examples: rrnB T1/T2 terminator, tonB terminator, natural sigma70 promoters (lac, ara, trp).

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

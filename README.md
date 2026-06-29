# VecTrap

> **Identify sequences of synthetic origin in bacterial genome assemblies and plasmid datasets.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Status: Active Development](https://img.shields.io/badge/status-active%20development-green.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844271.svg)](https://doi.org/10.5281/zenodo.20844271)

---

## Overview

VecTrap is a high-throughput Python pipeline for detecting and classifying sequences of synthetic origin in bacterial genome assemblies and plasmid datasets. Its primary use case is identifying assembled contigs that derive fully or partially from engineered laboratory constructs — cloning vectors, expression vectors, lentiviral transfer plasmids, and similar — rather than from the biological sample under study.

VecTrap uses a catalog-driven homology scanning strategy backed by a curated sequence database of more than 100,000 experimentally validated synthetic construct elements. The catalog covers the full diversity of engineered vector components found in prokaryotic, broad-host-range, and mammalian expression systems.

---

## Sequence Catalog

The VecTrap sequence catalog is deposited on Zenodo and must be downloaded separately before running the pipeline.

**DOI: [10.5281/zenodo.20844271](https://doi.org/10.5281/zenodo.20844271)**

The catalog consists of 14 FASTA files, one per feature type category, plus a manifest file. Download and index preparation is handled automatically by the `vectrap-build-db` command (see [Quick Start](#quick-start) below).

---

## Key Features

- **Catalog-driven homology detection** against a curated database of real synthetic vector sequences
- **Evidence-tiered classification** — every hit is assigned to `ENGINEERED`, `CONTEXT_DEPENDENT`, or `WEAK` with a `confidence` level (`HIGH` / `MEDIUM` / `LOW`) and a one-sentence human-readable `reasoning` string
- **Broad vector coverage** — prokaryotic, broad-host-range, lentiviral, and mammalian expression vector backbones
- **Bidirectional strand scanning** with precise coordinate mapping back to the forward sense strand
- **Three output levels** — per-hit evidence table, per-contig verdict table, and a per-sample summary row in `vectrap_summary.tsv`
- **Per-sample error isolation** — a failure on one sample prints a warning and writes an `ERROR` row to the summary; the rest of the batch continues
- **Modular architecture** — scanning, scoring, and reporting are fully decoupled

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
| T7 expression system elements | CONTEXT_DEPENDENT | T7 terminator, T7 tag |
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
        reporter.py
    catalogs/
        README.md      catalog acquisition and tier reference
pyproject.toml
requirements.txt
```

### Scanning pipeline

```
catalog indexes (built once)
    combined_long.mmi        mappy minimap2 index — sequences >= 50 bp
    short_index.pkl          Aho-Corasick k-mer automaton — sequences < 50 bp
    catalog_metadata.pkl     tier / confidence / reasoning per catalog entry

per-sample run
    1. read_fasta (single pass) → contigs + contig_lengths
    2. mappy aligner           → HomologyHit list (long sequences)
    3. Aho-Corasick scanner    → HomologyHit list (short sequences)
    4. metadata enrichment     → adds tier, confidence, reasoning, label to each hit
    5. scorer.summarize()      → ContigSummary list
    6. reporter.build_sample_row() → summary dict
    7. write TSVs
```

In batch mode the mappy index, k-mer automaton, and metadata dict are loaded **once** and reused across all samples.

---

## Installation

VecTrap requires Python 3.8 or higher.

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
pip install .
```

**Install from a specific branch:**

```bash
pip install "git+https://github.com/rustam-bioinfo/vectrap.git@<branch-name>"
```

**Editable install for development:**

```bash
git clone https://github.com/rustam-bioinfo/vectrap.git
cd vectrap
pip install -e .
```

Once published on PyPI:

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

# Run on a single assembly
vectrap -i assembly.fasta -o results/

# Run on a directory of assemblies (batch mode)
vectrap -i genomes/ -o results/

# Verbose output (per-stage progress + full tracebacks on error)
vectrap -i assembly.fasta -o results/ -v
```

---

## CLI Reference

```
vectrap -i FASTA|DIR -o DIR [options]

Required:
  -i, --input FASTA|DIR         Input assembly FASTA or directory of FASTA files.
  -o, --output DIR              Output directory. Created if it does not exist.

Optional:
  -c, --catalog-dir DIR         Catalog index directory (default: bundled vectrap/catalogs/).
  --min-identity FLOAT          Minimum sequence identity for mappy hits (default: 0.90).
  --min-coverage FLOAT          Minimum fraction of catalog sequence covered (default: 0.80).
  --min-len INT                 Length threshold (bp) separating mappy from k-mer scanner (default: 50).
  --best-n INT                  Max mappy alignments per query sequence (default: 10).
  --min-engineered-contamination INT
                                Engineered-hit count threshold for CONTAMINATION verdict (default: 3).
  --min-context-suspected INT   Context-dependent-hit count for SUSPECTED verdict (default: 1).
  -t, --threads INT             Aligner threads (default: 4).
  -v, --verbose                 Detailed per-stage progress and full tracebacks to stderr.
  --version                     Show version and exit.

Accepted FASTA extensions: .fasta  .fa  .fna  .fasta.gz  .fa.gz  .fna.gz
```

---

## Output Files

VecTrap writes three file types to the output directory.

### `<sample>_hits.tsv` — per-hit evidence table

One row per homology hit.

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header of the originating contig |
| `start` | 0-based start coordinate on the forward strand |
| `end` | 0-based end coordinate on the forward strand |
| `strand` | `+` for forward, `-` for reverse complement |
| `length` | Hit length in bp |
| `identity` | Per-base sequence identity (mappy hits only; `1.0` for k-mer hits) |
| `coverage` | Fraction of catalog sequence covered by the hit |
| `catalog_id` | Matched catalog entry identifier |
| `source` | `mappy` or `kmer` |
| `feature_type` | Functional category (e.g. `promoter`, `LTR`, `misc_feature`) |
| `label` | Human-readable element name (e.g. `CMV enhancer`, `HIV-1 Psi`) |
| `tier` | `ENGINEERED`, `CONTEXT_DEPENDENT`, or `WEAK` |
| `confidence` | `HIGH`, `MEDIUM`, or `LOW` |
| `reasoning` | One-sentence rationale for the tier assignment |

### `<sample>_contig_verdicts.tsv` — per-contig summary

One row per contig that has at least one hit.

| Column | Description |
| :--- | :--- |
| `contig` | FASTA header |
| `contig_length` | Length in bp |
| `total_hits` | Total hits on this contig |
| `engineered_hits` | Hits with `tier = ENGINEERED` |
| `context_dependent_hits` | Hits with `tier = CONTEXT_DEPENDENT` |
| `weak_hits` | Hits with `tier = WEAK` |
| `unannotated_hits` | Hits with no metadata entry |
| `unique_feature_types` | Number of distinct feature type categories |
| `unique_labels` | Number of distinct element labels |
| `covered_bp` | Merged non-overlapping bp covered by hits |
| `covered_fraction` | `covered_bp / contig_length` |
| `top_label` | Label of the hit with the highest identity |
| `top_tier` | Tier of the top hit |
| `top_confidence` | Confidence of the top hit |
| `evidence_summary` | Human-readable breakdown (e.g. `ENGINEERED[3]: CMV enhancer, HIV-1 Psi`) |

### `vectrap_summary.tsv` — per-sample summary

One row per input sample, written at the end of the run. Samples that fail produce an `ERROR` row with the error message in `top_labels`; all other samples in the batch are unaffected.

| Column | Description |
| :--- | :--- |
| `sample` | Sample stem name |
| `total_contigs` | Total contigs in assembly |
| `total_bp` | Total assembly length in bp |
| `total_hits` | Total hits across all contigs |
| `contigs_with_hits` | Number of contigs with at least one hit |
| `engineered_hits` | Total ENGINEERED hits |
| `context_dependent_hits` | Total CONTEXT_DEPENDENT hits |
| `weak_hits` | Total WEAK hits |
| `sample_verdict` | `CONTAMINATION`, `SUSPECTED`, `CLEAN`, or `ERROR` |
| `top_labels` | Top 5 matched element labels by hit count (or error message) |
| `unique_feature_types` | Number of distinct feature categories hit |
| `unique_labels` | Number of distinct element labels hit |
| `covered_bp` | Total merged covered bp across all contigs |
| `covered_fraction` | `covered_bp / total_bp` |
| `top_contig` | Contig with highest `covered_fraction` |
| `top_contig_length` | Length of `top_contig` |
| `top_contig_covered_fraction` | `covered_fraction` of `top_contig` |

---

## Evidence Classification

VecTrap is designed for **bacterial genome assemblies**. All tier assignments reflect whether a sequence could plausibly occur in a wild-type bacterial chromosome or natural bacterial plasmid.

Every catalog hit is assigned to one of three evidence tiers:

- **`ENGINEERED`** — The matched sequence has no natural occurrence in any wild-type bacterium. It was designed in a laboratory, is derived from eukaryotic viruses or mammalian systems, or originates from contexts (phage RNA polymerase, mammalian cell culture) that have no plausible route into a bacterial genome without deliberate cloning. A single high-confidence ENGINEERED hit is sufficient to flag a contig as synthetic. Examples: T7 promoter, CMV enhancer, SV40 ori, f1/M13 ori, BGH polyA signal, EGFP, Cas9, attB/FRT recombination sites, IRES, WPRE, retroviral LTR.

- **`CONTEXT_DEPENDENT`** — The matched sequence exists naturally in bacteria or natural bacterial plasmids, but its presence in an assembly may indicate synthetic contamination depending on the host organism and genomic context. A ColE1 origin in an *E. coli* assembly is unremarkable; the same origin in a *Salmonella* or *Bacillus* assembly warrants scrutiny. Antibiotic resistance genes (KanR, AmpR, CmR) are endemic in environmental and clinical bacteria via horizontal gene transfer and are not diagnostic alone. Their evidential weight comes from **functional class co-occurrence**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter is almost certainly synthetic. Examples: ColE1 ori, p15A ori, KanR, AmpR, CmR, tac/trc promoters, T7 terminator, lacI/lacO.

- **`WEAK`** — The matched sequence is widespread across wild-type bacteria with low specificity for synthetic origin. These elements contribute to a verdict only when co-occurring with stronger evidence from other tiers. Examples: rrnB T1/T2 terminator, tonB terminator, natural sigma70 promoters (lac, ara, trp).

### Sample verdict thresholds

| Condition | Verdict |
| :--- | :--- |
| `engineered_hits >= --min-engineered-contamination` (default 3) | `CONTAMINATION` |
| `engineered_hits > 0` OR `context_dependent_hits >= --min-context-suspected` (default 1) | `SUSPECTED` |
| Neither of the above | `CLEAN` |
| Unhandled exception during processing | `ERROR` |

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

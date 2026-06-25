# VecTrap — Coding Session Context Document

Use this document to onboard a coding model to the VecTrap project without repeating
background. Paste the entire document at the start of a session, then state the task.

---

## What VecTrap is

VecTrap is a Python command-line pipeline that detects sequences of synthetic origin
(vector contamination) in bacterial genome assemblies and plasmid datasets. It is
catalog-driven: it scans input FASTA contigs against a curated database of >100,000
synthetic construct sequences and classifies each contig as `VECTOR`, `SUSPECTED`, or
`CLEAN`.

Repository: https://github.com/rustam-bioinfo/vectrap  
Sequence catalog (Zenodo): https://doi.org/10.5281/zenodo.20844271  
License: MIT  
Python requirement: 3.8+

---

## Repository layout

```
vectrap/
    __init__.py
    cli/
        build_db.py       # vectrap-build-db entry point — COMPLETE
        run.py            # vectrap entry point — partially complete (arg parser done, wiring TODO)
    modules/
        utils.py          # COMPLETE
        homology_scanner.py   # COMPLETE
        scorer.py         # NOT YET WRITTEN
    catalogs/
        README.md         # catalog format + tier reference
pyproject.toml
requirements.txt
TODO.md
README.md
INFO.md
CONTEXT.md               # this file
```

---

## Module status

### `vectrap/modules/utils.py` — COMPLETE

Provides:
- `read_fasta(path)` — plain and gzipped FASTA reader, yields `(header, seq)` tuples
- `rev_comp(seq)` — reverse complement
- `open_text(path)` — transparent gz/plain file opener
- `write_fasta(records, path)` — plain and gzipped FASTA writer

### `vectrap/cli/build_db.py` — COMPLETE

- `--download` mode: fetches catalog FASTA files from Zenodo DOI, verifies MD5 checksums
  against `catalog_manifest.tsv`, builds indexes
- `--catalog-dir` mode: uses locally provided FASTA files
- Length-based routing: sequences >= 50 bp → minimap2 index (`combined_long.mmi`);
  sequences < 50 bp → k-mer hash pickle (`short_index.pkl`)

### `vectrap/modules/homology_scanner.py` — COMPLETE

Public interface: `scan(query_fasta, catalog_dir, min_identity, min_coverage, threads, verbose) -> List[HomologyHit]`

**`HomologyHit` dataclass fields:**
- `contig` (str) — contig name from FASTA header
- `start` (int) — 0-based start, always on forward strand
- `end` (int) — 0-based exclusive end, always on forward strand
- `strand` (str) — `"+"` or `"-"`
- `identity` (float) — fraction in [0, 1]
- `coverage` (float) — fraction in [0, 1], measured on the **catalog** (target) sequence
- `catalog_id` (str) — matched catalog entry identifier
- `source` (str) — `"mappy"` or `"kmer"`
- `.length` property — `end - start`

**Two internal scanners:**

1. `_run_minimap2(query_fasta, catalog_dir, min_identity, min_coverage, threads)` — uses
   `minimap2 -c --cs=short` against `combined_long.mmi`; parses PAF; identity derived
   from `de:f:` tag → `dv:f:` → NM/aln_len fallback; coverage = `(t_end - t_start) / t_len`

2. `_run_kmer_scanner(query_fasta, catalog_dir)` — loads `short_index.pkl` (a
   `{seq: [catalog_id, ...]}` dict); slides over each contig forward and RC; maps RC
   coordinates back to forward strand; deduplicates palindromes; returns identity=1.0,
   coverage=1.0 for all exact hits

**Coordinate convention (strictly enforced):** All coordinates are 0-based, half-open,
on the forward strand of the original contig, regardless of match strand.

### `vectrap/modules/scorer.py` — NOT YET WRITTEN

This is the next module to implement. See the Scorer specification section below.

### `vectrap/cli/run.py` — PARTIALLY COMPLETE

Argument parser is done (`-i`, `-o`, `-c`, `--min-identity`, `--min-coverage`,
`--version`, `-v`). Scanner and scorer are not yet wired. Output writing is not yet
implemented.

---

## Evidence tier model

All tiers are defined relative to the **bacterial genome context**. The question asked
for each catalog entry is: could this sequence plausibly appear in a wild-type bacterial
chromosome or natural bacterial plasmid without deliberate cloning?

### ENGINEERED

No natural occurrence in any wild-type bacterium. A single high-confidence ENGINEERED
hit is sufficient for a `VECTOR` verdict.

Examples: T7/SP6/T3 promoters, CMV/EF1α/CAG promoters, SV40 ori, f1/M13 ori, BGH/SV40
polyA, EGFP/mCherry and all fluorescent proteins, Cas9/Cas12/sgRNA scaffolds, Cre/Flp
recombinases, attB/attP/loxP/FRT sites, IRES/WPRE/P2A/T2A, sequencing primer binding
sites (M13/T7/pUC), retroviral LTR sequences, HIV PSI.

### CONTEXT_DEPENDENT

Exists naturally in bacteria or natural bacterial plasmids. Whether it signals synthetic
contamination depends on the host organism and genomic context. Evidential weight comes
from **functional class co-occurrence**: a contig carrying a replication origin + a
resistance marker + a promoter is almost certainly synthetic; a contig with only one of
these is ambiguous.

Examples: ColE1/p15A/pSC101/RK2/pBBR1 ori, oriT (IncP/IncQ/IncW), KanR/AmpR/CmR/TetR
(Tn-derived), lacI/lacO arrays, tac/trc promoters, T7 terminator, T7 tag, transposon
borders (Tn3/Tn5/Tn10/Tn903), IS-element borders, protein-binding operators in natural
form.

### WEAK

Widespread in wild-type bacterial genomes; low specificity for synthetic origin. A WEAK
hit alone is **never** sufficient for `VECTOR` or `SUSPECTED`. It contributes only when
ENGINEERED or CONTEXT_DEPENDENT evidence is already present on the same contig.

Examples: rrnB T1/T2 terminator, tonB terminator, lambda t0, unmodified lac/ara/trp/phoA
promoters, rRNA/tRNA genes.

---

## Scorer specification (next module to write)

File: `vectrap/modules/scorer.py`

### Inputs

A `List[HomologyHit]` from `homology_scanner.scan()` plus a catalog metadata source
(either a dict pre-loaded by the caller or derived from the catalog directory) that maps
`catalog_id → (evidence_tier, feature_class)`.

The two key per-hit metadata fields needed for scoring:
- `evidence_tier`: one of `"ENGINEERED"`, `"CONTEXT_DEPENDENT"`, `"WEAK"`
- `feature_class`: functional category string, e.g. `"rep_origin"`, `"resistance_marker"`,
  `"promoter"`, `"terminator"`, `"recombination_site"`, `"primer_bind"`, `"regulatory"`,
  `"LTR"`, `"mobile_element"`, `"oriT"`, `"CDS"`, `"protein_bind"`, `"enhancer"`,
  `"misc_feature"`

These come from the catalog manifest (`catalog_manifest.tsv`) or from the FASTA header
metadata. The scorer must load this mapping from the catalog directory.

### Output dataclass `ContigVerdict`

```python
@dataclass
class ContigVerdict:
    contig: str
    classification: str          # "VECTOR", "SUSPECTED", or "CLEAN"
    engineered_hits: int         # count of ENGINEERED-tier hits
    context_dependent_classes: int  # count of distinct feature_class values among CONTEXT_DEPENDENT hits
    weak_hits: int               # count of WEAK-tier hits
    evidence_summary: str        # human-readable one-line summary
```

### Classification logic

```
VECTOR    if:  engineered_hits >= 1
          OR   context_dependent_classes >= 2

SUSPECTED if:  context_dependent_classes == 1
          OR   weak_hits >= 1 AND context_dependent_classes >= 1

CLEAN     otherwise
```

The threshold of `context_dependent_classes >= 2` for `VECTOR` (rather than a simple
hit count) is the core design decision: functional class diversity, not raw hit count,
is the signal for CONTEXT_DEPENDENT evidence.

### Public interface

```python
def score(
    hits: List[HomologyHit],
    catalog_dir: str | Path,
) -> List[ContigVerdict]:
    ...
```

Returns one `ContigVerdict` per contig present in the input hits. Contigs with zero hits
are not included (they are CLEAN by absence).

### Catalog metadata loading

The scorer needs to map `catalog_id → (evidence_tier, feature_class)`. The catalog
manifest at `{catalog_dir}/catalog_manifest.tsv` contains this mapping. If a
`catalog_id` is not found in the manifest, log a warning and skip the hit rather than
crashing.

---

## Output format (for `run.py` wiring)

### Per-hit table: `{output_dir}/{stem}.all_hits.tsv`

Tab-separated, one row per `HomologyHit`, columns:

```
contig  start_0based  end_0based  strand  evidence_tier  feature_class  catalog_id  identity  coverage  source
```

### Per-contig verdict table: `{output_dir}/{stem}.verdicts.tsv`

Tab-separated, one row per `ContigVerdict`, columns:

```
contig  classification  engineered_hits  context_dependent_classes  weak_hits  evidence_summary
```

Contigs that are absent from all_hits (no hits at all) should appear in verdicts.tsv
as `CLEAN` with all numeric columns = 0 and evidence_summary = `"no hits"`.

---

## Coding conventions used in this project

- Python 3.8+; no walrus operator, no `match` statements
- `dataclasses.dataclass` for data objects (no Pydantic)
- No third-party logging framework; verbose output goes to `sys.stderr` with
  `print(f"[vectrap] {msg}", file=sys.stderr)` when `verbose=True`
- Type hints on all public functions; `from __future__ import annotations` at top of
  each file
- Pathlib (`Path`) used throughout; `str | Path` accepted in function signatures,
  converted to `Path` internally at the top of the function
- No use of `assert` for runtime validation; use explicit `raise ValueError(...)` or
  `raise FileNotFoundError(...)`
- Tests will go in `tests/` using pytest; not written yet

---

## What is still TODO (from TODO.md)

**Immediate next steps:**
1. Write `vectrap/modules/scorer.py` (spec above)
2. Wire scanner + scorer into `vectrap/cli/run.py`, implement output writing
3. Update `requirements.txt` with finalized dependencies (mappy, pyahocorasick)
4. Add `.gitignore` entries for `*.mmi` and `*.pkl`

**After that:**
5. `tests/test_scanner.py` — unit tests on positive/negative control sequences
6. `tests/test_scorer.py` — unit tests for classification logic
7. `tests/test_cli.py` — integration test for full pipeline run
8. Finalize `README.md` once CLI is stable
9. Add `CHANGELOG.md` and `CONTRIBUTING.md`
10. Write methods section for manuscript
11. Publish to PyPI

---

## Key design decisions (do not change without understanding the rationale)

| Decision | Rationale |
| :--- | :--- |
| Three tiers (ENGINEERED / CONTEXT_DEPENDENT / WEAK) not two | RECRUITED was too coarse; bacterial context matters for natural-origin sequences |
| f1/M13 ori are ENGINEERED (not CONTEXT_DEPENDENT) | Phagemid origins; do not occur in bacterial chromosomes or natural conjugative plasmids |
| T7 promoter is ENGINEERED; T7 terminator is CONTEXT_DEPENDENT | Promoter requires non-endogenous T7 RNAP; terminator occurs in natural phage genomes |
| rrnB / tonB terminators are WEAK (not CONTEXT_DEPENDENT) | Too widespread in wild-type bacteria for meaningful signal |
| Coverage measured on catalog sequence, not query contig | Answers "how much of the known element is present" — the biologically useful question |
| All coordinates in forward-strand 0-based half-open | Single coordinate system across both scanners; avoids strand-specific interval logic |
| Functional class diversity (not hit count) drives CONTEXT_DEPENDENT verdict | A contig with 20 ColE1-related hits is not more suspicious than one with 1; two different functional classes is the real signal |
| WEAK hits alone never trigger VECTOR or SUSPECTED | Prevents false positives from sequences found in both lab vectors and wild-type genomes |

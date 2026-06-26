# VecTrap — Coding Session Context Document

Use this document to onboard a coding model to the VecTrap project without repeating
background. Paste the entire document at the start of a session, then state the task.

---

## What VecTrap is

VecTrap is a Python command-line pipeline that detects sequences of synthetic origin
(vector contamination) in bacterial genome assemblies and plasmid datasets. It is
catalog-driven: it scans input FASTA contigs against a curated database of >100,000
synthetic construct sequences and classifies each sample as `CONTAMINATION`, `SUSPECTED`,
or `CLEAN`.

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
        build_db.py           # vectrap-build-db entry point — COMPLETE
        run.py                # vectrap entry point — COMPLETE
    modules/
        utils.py              # COMPLETE
        homology_scanner.py   # COMPLETE
        scorer.py             # COMPLETE
        reporter.py           # COMPLETE
    catalogs/
        README.md             # catalog format + tier reference
pyproject.toml
requirements.txt
TODO.md
README.md
INFO.md
CONTEXT.md                    # this file
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
  against `catalog_manifest.tsv`, builds all three indexes
- `--catalog-dir` mode: uses locally provided FASTA files
- Length-based routing: sequences >= 50 bp → minimap2 index (`combined_long.mmi`);
  sequences < 50 bp → Aho-Corasick k-mer automaton (`short_index.pkl`)
- Also builds `catalog_metadata.pkl` from `catalog_manifest.tsv`

### `vectrap/modules/homology_scanner.py` — COMPLETE

Public interface:
```python
scan(
    query_fasta: str | Path,
    catalog_dir: str | Path,
    min_identity: float = 0.90,
    min_coverage: float = 0.80,
    threads: int = 4,
    verbose: bool = False,
) -> List[HomologyHit]
```

**`HomologyHit` dataclass fields:**

| Field | Type | Description |
| :--- | :--- | :--- |
| `contig` | `str` | Contig name from FASTA header (first token) |
| `start` | `int` | 0-based start, always on forward strand |
| `end` | `int` | 0-based exclusive end, always on forward strand |
| `strand` | `str` | `"+"` or `"-"` |
| `identity` | `float` | Fraction in `[0, 1]` |
| `coverage` | `float` | Fraction in `[0, 1]`, measured on the **catalog** (target) sequence |
| `catalog_id` | `str` | Matched catalog entry identifier |
| `source` | `str` | `"mappy"` or `"kmer"` |
| `feature_type` | `str` | Functional category from catalog metadata; `""` if not found |
| `label` | `str` | Human-readable element name from catalog metadata; `""` if not found |
| `tier` | `str` | `"ENGINEERED"`, `"CONTEXT_DEPENDENT"`, `"WEAK"`, or `""` |
| `confidence` | `str` | `"HIGH"`, `"MEDIUM"`, `"LOW"`, or `""` |
| `reasoning` | `str` | One-sentence rationale from catalog metadata; `""` if not found |

Convenience property: `.length` = `end - start`

The `feature_type`, `label`, `tier`, `confidence`, and `reasoning` fields default to `""`
at scan time and are populated during a metadata enrichment step by looking up
`catalog_id` in `catalog_metadata.pkl`.

**Two internal scanners:**

1. `_run_minimap2(...)` — uses mappy against `combined_long.mmi`; filters by
   `min_identity` and `min_coverage`; identity = `mlen / blen`; coverage =
   `(r_en - r_st) / ctg_len`

2. `_run_kmer_scanner(...)` — loads `short_index.pkl`, builds Aho-Corasick automaton,
   single-pass per contig; RC patterns inserted at build time; all exact hits get
   `identity = 1.0`, `coverage = 1.0`

**Coordinate convention (strictly enforced):** All coordinates are 0-based, half-open,
on the forward strand of the original contig, regardless of match strand.

**Index loading:** The mappy aligner index and the Aho-Corasick automaton are loaded
once per batch session and passed into `scan()` — they are not reloaded per sample.

### `vectrap/modules/scorer.py` — COMPLETE

Public interface:
```python
summarize(
    hits: List[HomologyHit],
    contig_lengths: Dict[str, int] | None = None,
) -> List[ContigSummary]
```

**`ContigSummary` dataclass fields:**

| Field | Type | Description |
| :--- | :--- | :--- |
| `contig` | `str` | Contig name |
| `contig_length` | `int` | Length in bp; 0 if not supplied |
| `total_hits` | `int` | Total hit count (mappy + kmer) |
| `engineered_hits` | `int` | Hits with `tier == "ENGINEERED"` |
| `context_dependent_hits` | `int` | Hits with `tier == "CONTEXT_DEPENDENT"` |
| `weak_hits` | `int` | Hits with `tier == "WEAK"` |
| `unannotated_hits` | `int` | Hits with no metadata (`tier == ""`) |
| `unique_feature_types` | `int` | Distinct `feature_type` values |
| `unique_labels` | `int` | Distinct `label` values |
| `covered_bp` | `int` | Non-overlapping bp covered (interval-merged) |
| `covered_fraction` | `float` | `covered_bp / contig_length`; 0.0 if unknown |
| `top_label` | `str` | Label of hit with highest `identity × coverage` |
| `top_tier` | `str` | Tier of the top hit |
| `top_confidence` | `str` | Confidence of the top hit |
| `evidence_summary` | `str` | e.g. `ENG[3]: CMV enhancer, HIV-1 Psi; CD[1]: ColE1 ori` |

`covered_bp` is computed by merging all hit intervals with a sort-and-sweep algorithm.
Summaries are returned sorted by `(-engineered_hits, -context_dependent_hits, contig)`.

### `vectrap/modules/reporter.py` — COMPLETE

Aggregates `List[ContigSummary]` into a sample-level row dict for `vectrap_summary.tsv`.

**Sample verdict logic:**

| Condition | Verdict |
| :--- | :--- |
| `engineered_hits >= --min-engineered-contamination` (default 3) | `CONTAMINATION` |
| `engineered_hits > 0` OR `context_dependent_hits >= --min-context-suspected` (default 1) | `SUSPECTED` |
| Neither of the above | `CLEAN` |
| Unhandled exception | `ERROR` |

### `vectrap/cli/run.py` — COMPLETE

Full pipeline wiring: argument parsing → catalog index loading (once) → per-sample
loop → hit scanning → per-sample TSV writing → verdict/summary aggregation →
`vectrap_summary.tsv`.

**Per-sample error isolation:** each sample runs inside a try/except; on failure, partial
output files are deleted, an `ERROR` row is added to the summary, and the batch continues.

**Batch mode:** `-i <directory>` discovers all FASTA files (`.fasta`, `.fa`, `.fna`,
`.fasta.gz`, `.fa.gz`, `.fna.gz`) and processes them in sequence. Indexes are loaded once.

---

## Output files

### `{output_dir}/{stem}_hits.tsv`

One row per `HomologyHit`. Columns:

```
contig  start  end  strand  length  feature_type  label  tier  confidence  reasoning
catalog_id  identity  coverage  source
```

### `{output_dir}/{stem}_contig_verdicts.tsv`

One row per `ContigSummary`. Columns:

```
contig  contig_length  total_hits  engineered_hits  context_dependent_hits  weak_hits
unannotated_hits  unique_feature_types  unique_labels  covered_bp  covered_fraction
top_label  top_tier  top_confidence  evidence_summary
```

### `{output_dir}/vectrap_summary.tsv`

One row per sample. Includes: `sample`, `sample_verdict`, `total_hits`,
`engineered_hits`, `context_dependent_hits`, `weak_hits`, `contaminated_contigs`,
`covered_bp`, `top_labels`, and per-verdict contig counts.

---

## Evidence tier model

All tiers are defined relative to the **bacterial genome context**. The question asked
for each catalog entry is: could this sequence plausibly appear in a wild-type bacterial
chromosome or natural bacterial plasmid without deliberate cloning?

### ENGINEERED

No natural occurrence in any wild-type bacterium. A single HIGH-confidence ENGINEERED
hit triggers `SUSPECTED`; three or more trigger `CONTAMINATION` (configurable).

Examples: T7/SP6/T3 promoters, CMV/EF1α/CAG promoters, SV40 ori, f1/M13 ori, BGH/SV40
polyA, EGFP/mCherry and all fluorescent proteins, Cas9/Cas12/sgRNA scaffolds, Cre/Flp
recombinases, attB/attP/loxP/FRT sites, IRES/WPRE/P2A/T2A, universal sequencing primer
sites (M13/T7/pUC), retroviral LTR sequences (HIV-1, MMLV, HTLV-1, MSCV), HIV-1 PSI,
RRE, cPPT/CTS, CMV/SV40/hr5 enhancers.

### CONTEXT_DEPENDENT

Exists naturally in bacteria or natural bacterial plasmids. Evidential weight comes from
**functional class co-occurrence**: a contig carrying a replication origin + resistance
marker + promoter is almost certainly synthetic; a contig with only one of these is
ambiguous.

Examples: ColE1/p15A/pSC101/RK2/pBBR1 ori, oriT (IncP/IncQ/IncW), KanR/AmpR/CmR/TetR
(Tn-derived), lacI/lacO arrays, tac/trc promoters, T7 terminator, T7 tag, transposon
borders (Tn3/Tn5/Tn10/Tn903), IS-element borders, Agrobacterium T-DNA border repeats
(outside native context), pBR322 bom, protein-binding operators in natural form.

### WEAK

Widespread in wild-type bacterial genomes; low specificity for synthetic origin. A WEAK
hit alone is **never** sufficient for `CONTAMINATION` or `SUSPECTED`.

Examples: rrnB T1/T2 terminator, tonB terminator, lambda t0, unmodified lac/ara/trp/phoA
promoters, short regulatory motifs (human Pol II pause site).

### Confidence levels

Each catalog entry has a `confidence` field:
- `HIGH` — unambiguous evidence for its tier
- `MEDIUM` — alternative explanation possible in some contexts
- `LOW` — short, degenerate, or common enough to appear by chance

---

## Catalog metadata

`catalog_manifest.tsv` has one row per catalog entry with columns:
`catalog_id`, `feature_type`, `label`, `tier`, `confidence`, `reasoning`, `length`, `md5`

`catalog_metadata.pkl` is a compiled dict:
```python
{ "catalog_id": {"feature_type": ..., "label": ..., "tier": ..., "confidence": ..., "reasoning": ...} }
```
Built by `vectrap-build-db`; loaded once per batch run; used by `homology_scanner.scan()`
to enrich hits.

---

## Coding conventions

- Python 3.8+; no walrus operator, no `match` statements
- `dataclasses.dataclass` for data objects (no Pydantic)
- No third-party logging framework; verbose output goes to `sys.stderr` with timestamps
  when `verbose=True`
- Type hints on all public functions; `from __future__ import annotations` at top of each
  file
- `pathlib.Path` used throughout; `str | Path` accepted in function signatures, converted
  to `Path` internally
- No `assert` for runtime validation; use `raise ValueError(...)` or
  `raise FileNotFoundError(...)`
- Tests go in `tests/` using pytest

---

## Key design decisions

| Decision | Rationale |
| :--- | :--- |
| Three tiers (ENGINEERED / CONTEXT_DEPENDENT / WEAK) | Bacterial context matters; CONTEXT_DEPENDENT sequences are real positives in some organisms but not others |
| f1/M13 ori are ENGINEERED | Phagemid origins; do not occur in bacterial chromosomes or natural conjugative plasmids |
| T7 promoter is ENGINEERED; T7 terminator is CONTEXT_DEPENDENT | Promoter requires non-endogenous T7 RNAP; terminator occurs in natural phage genomes |
| rrnB / tonB terminators are WEAK | Too widespread in wild-type bacteria for meaningful signal |
| Coverage measured on catalog sequence, not query contig | Answers “how much of the known element is present” — the biologically useful question |
| All coordinates forward-strand 0-based half-open | Single coordinate system across both scanners; avoids strand-specific interval logic |
| Functional class diversity drives CONTEXT_DEPENDENT signal | A contig with 20 ColE1 hits is not more suspicious than one with 1; two different functional classes is the signal |
| WEAK hits alone never trigger a verdict | Prevents false positives from sequences in both lab vectors and wild-type genomes |
| Tier / confidence / reasoning come from catalog manifest, not heuristics | Evidence model is explicit and auditable; scanner looks up, does not infer |
| Indexes loaded once per batch run | Avoids redundant mappy index loads and automaton builds for large batches |
| Per-sample error isolation | One malformed FASTA does not abort a 1000-sample batch |
| CONTAMINATION / SUSPECTED / CLEAN (not VECTOR) | Clearer semantics; CONTAMINATION implies confirmed multi-hit evidence |

---

## What is still TODO

1. `tests/test_scanner.py` — unit tests on positive/negative control sequences
2. `tests/test_scorer.py` — unit tests for `summarize()` and verdict logic
3. `tests/test_cli.py` — integration test for full pipeline run
4. `CHANGELOG.md`
5. `CONTRIBUTING.md`
6. Methods section for manuscript
7. PyPI publication

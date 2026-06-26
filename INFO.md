# VecTrap — Technical Reference

This file is a technical reference for VecTrap internals. It covers algorithms, data structures, coordinate conventions, the evidence model, scoring logic, output formats, and performance characteristics. See `README.md` for installation and quick-start instructions.

---

## Purpose

VecTrap detects sequences of synthetic origin in assembled nucleotide sequences, with the current focus on bacterial genome assemblies and plasmid datasets. The core question it answers is: does any part of this assembled sequence derive from an engineered laboratory construct rather than from the biological sample under study?

VecTrap is designed for cases where a deposited or locally assembled sequence may contain engineered backbone fragments — replication origins, selectable markers, recombination sites, synthetic promoters, primer binding sites, terminators, or other laboratory-use elements. It does not rely on database annotations or assembly metadata; it works directly from sequence evidence.

---

## High-level workflow

```
catalog indexes (built once per session)
    combined_long.mmi        mappy minimap2 index — sequences >= 50 bp
    short_index.pkl          Aho-Corasick k-mer automaton — sequences < 50 bp
    catalog_metadata.pkl     tier / confidence / reasoning per catalog entry

per-sample run
    1. read_fasta (single pass) → contigs list + contig_lengths dict
    2. mappy aligner           → HomologyHit list (long sequences)
    3. Aho-Corasick scanner    → HomologyHit list (short sequences)
    4. metadata enrichment     → adds tier, confidence, reasoning, label to each hit
    5. scorer.summarize()      → ContigSummary list
    6. reporter.build_sample_row() → summary dict
    7. write TSVs
```

In batch mode (`-i <directory>`) the mappy index, k-mer automaton, and metadata dict are loaded **once** before the sample loop and reused across all samples.

---

## Inputs

VecTrap requires:

- An input FASTA file containing one or more contigs (plain or gzipped; `.fasta`, `.fa`, `.fna`, `.fasta.gz`, `.fa.gz`, `.fna.gz`).
- A prepared catalog directory built by `vectrap-build-db`.

The catalog directory must contain:

- `combined_long.mmi` — minimap2/mappy index for catalog entries >= 50 bp.
- `short_index.pkl` — serialised dict of short catalog sequences for k-mer matching.
- `catalog_metadata.pkl` — dict mapping `catalog_id` → `{feature_type, label, tier, confidence, reasoning}`.

---

## Why two scanners are used

Short exact motifs and longer homologous segments behave very differently computationally.

### Long sequences

Long catalog entries are best handled by an aligner that tolerates mismatches and partial divergence. VecTrap uses **mappy**, the Python binding to **minimap2**. This is appropriate for elements such as longer origins, resistance markers, or regulatory modules where exact matching would be too strict due to natural sequence variation across strains and vector generations.

### Short sequences

Short elements are a poor fit for standard alignment because they are abundant, exact matching is usually sufficient, and naive repeated substring search becomes slow. For these, VecTrap uses an **Aho-Corasick automaton**, which finds all short patterns simultaneously in a single pass through each contig.

This split design gives both sensitivity and speed.

---

## Algorithms

### 1. mappy / minimap2 for long catalog entries

For catalog entries at or above the internal length threshold, VecTrap queries the input contigs against a pre-built `.mmi` index using `mappy.Aligner`.

Hits are filtered by two configurable thresholds:

- **Minimum sequence identity** (`--min-identity`, default `0.90`)
- **Minimum catalog coverage** (`--min-coverage`, default `0.80`)

Identity is computed from the alignment object as:

```
identity = mlen / blen
```

where `mlen` = matching bases and `blen` = alignment block length.

Coverage is computed relative to the matched catalog sequence:

```
coverage = (r_en - r_st) / ctg_len
```

where `r_st`, `r_en` are start/end on the catalog entry and `ctg_len` is the total catalog entry length. This answers the biologically useful question: how much of the known synthetic element is covered by the hit?

### 2. Aho-Corasick automaton for short catalog entries

For short catalog entries, VecTrap uses the **Aho-Corasick** string matching algorithm.

Build phase (once per session):

1. All short catalog patterns are inserted into a trie, along with their reverse complements.
2. Failure links are added to turn the trie into a finite-state automaton.

Scan phase (per contig):

3. Each contig is streamed through the automaton in a single linear pass.
4. Every exact match to any stored pattern is emitted immediately with its `catalog_id` and strand.

This lets VecTrap find thousands of short motifs in one pass per contig rather than one pass per motif.

For exact short-pattern matches, `identity` and `coverage` are both reported as `1.0`.

---

## Strand handling

VecTrap scans both orientations.

For the short-pattern scanner, both the forward pattern and its reverse complement are inserted into the automaton. A single forward-strand pass over each contig therefore detects matches on either strand.

For long-sequence hits, strand is obtained directly from the mappy alignment output.

### Coordinate convention

All reported coordinates are:

- **0-based, half-open** (`start`, `end`) — following the Python/BED convention.
- **Mapped onto the forward strand of the original contig**, regardless of which strand the match biologically occurs on.

This unified coordinate system makes interval merging and downstream analysis straightforward because all hits live in the same reference frame.

---

## Internal hit representation

Both scanners produce `HomologyHit` objects with these fields:

| Field | Type | Description |
| :--- | :--- | :--- |
| `contig` | `str` | Contig name (FASTA header first token) |
| `start` | `int` | 0-based start on forward strand |
| `end` | `int` | 0-based exclusive end on forward strand |
| `strand` | `str` | `+` or `-` |
| `identity` | `float` | Sequence identity in `[0, 1]`; `1.0` for k-mer hits |
| `coverage` | `float` | Catalog entry coverage in `[0, 1]`; `1.0` for k-mer hits |
| `catalog_id` | `str` | Matched catalog entry identifier |
| `source` | `str` | `mappy` or `kmer` |
| `feature_type` | `str` | Functional category from catalog metadata |
| `label` | `str` | Human-readable element name from catalog metadata |
| `tier` | `str` | `ENGINEERED`, `CONTEXT_DEPENDENT`, `WEAK`, or `""` |
| `confidence` | `str` | `HIGH`, `MEDIUM`, `LOW`, or `""` |
| `reasoning` | `str` | One-sentence rationale for the tier assignment |

A convenience property `length = end - start` is also available.

The `feature_type`, `label`, `tier`, `confidence`, and `reasoning` fields default to `""` at scan time and are populated during the metadata enrichment step by looking up `catalog_id` in `catalog_metadata.pkl`.

---

## Length threshold split

VecTrap separates catalog entries at **50 bp**:

- `length < 50 bp` → short index (`short_index.pkl`, Aho-Corasick)
- `length >= 50 bp` → long index (`combined_long.mmi`, mappy)

This threshold is applied consistently during both database build (`vectrap-build-db`) and scan time. The two components must use the same threshold to avoid entries being missed or double-counted.

---

## Catalog metadata

Every catalog entry has a corresponding metadata record in `catalog_metadata.pkl`:

```python
{
    "catalog_id_1": {
        "feature_type": "LTR",
        "label":        "HIV-1 3' LTR",
        "tier":         "ENGINEERED",
        "confidence":   "HIGH",
        "reasoning":    "HIV-1 3' LTR is a retroviral regulatory sequence ...",
    },
    ...
}
```

This is the source of truth for tier, confidence, and reasoning — the scanner does not compute these; it looks them up. The metadata is assembled from `catalog_manifest.tsv` during `vectrap-build-db`.

---

## Evidence model

All tier assignments are made relative to the **bacterial genome context** — the defining question is whether a given sequence could plausibly appear in a wild-type bacterial chromosome or a natural bacterial plasmid without deliberate cloning.

Every catalog hit carries one of three evidence tiers.

### ENGINEERED tier

Sequences with no natural occurrence in any wild-type bacterium. They originate from eukaryotic viruses, mammalian cell biology, phage RNA polymerase systems, or are fully synthetic constructs with no plausible route into a bacterial genome except deliberate cloning. A single high-confidence ENGINEERED hit is sufficient to flag a contig.

Includes: eukaryotic and viral promoters (CMV, EF1α, PGK, CAG, CaMV35S), phage RNA polymerase promoters (T7, SP6, T3), eukaryotic replication origins (SV40 ori, EBV oriP), filamentous phage cloning origins (f1, M13), polyadenylation signals (BGH polyA, SV40 polyA), mammalian selectable markers, fluorescent proteins (EGFP, mCherry, Luciferase), CRISPR elements (Cas9, sgRNA scaffolds), recombinase systems (Cre, Flp, PhiC31), recombination sites (loxP, FRT, attB/attP), synthetic regulatory elements (IRES, WPRE, P2A/T2A), sequencing primer sites, retroviral/lentiviral elements (LTR, HIV PSI, RRE, cPPT/CTS).

### CONTEXT_DEPENDENT tier

Sequences that exist naturally in bacteria or natural bacterial plasmids. Their presence is not automatically diagnostic — it depends on the host organism and broader genomic context.

For example, a ColE1 origin in an *E. coli* assembly may reflect a natural endogenous plasmid. The same origin in a *Salmonella* or *Bacillus* assembly is anomalous. Antibiotic resistance genes (KanR, AmpR, CmR) are endemic via horizontal gene transfer and are not diagnostic alone. Their evidential weight comes from **functional class co-occurrence**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter from this tier is almost certainly synthetic.

Includes: bacterial replication origins (ColE1, pMB1, pUC ori, p15A, pSC101, RK2, pBBR1), transfer origins (oriT IncP/IncQ/IncW), antibiotic resistance markers (KanR, AmpR, CmR, TetR, GenR, SpecR), T7 expression system elements (T7 terminator, T7 tag), inducible expression components (lacI, lacO, tac, trc promoters), Agrobacterium T-DNA border repeats, mobile element remnants (Tn3, Tn5, Tn10, Tn903 ends).

### WEAK tier

Sequences widespread in wild-type bacterial genomes with low specificity for synthetic origin. Contribute to a verdict only when co-occurring with ENGINEERED or CONTEXT_DEPENDENT evidence.

Includes: rrnB T1/T2 terminator, tonB terminator, lambda t0, fd terminator; natural sigma70-dependent promoters (lac, ara, trp, phoA) in their unmodified form.

### Confidence levels

Each catalog entry also carries a `confidence` field (`HIGH`, `MEDIUM`, `LOW`) that reflects how unambiguous the tier assignment is:

- **HIGH** — the element's presence is unambiguous evidence for its tier (e.g., a viral LTR in a bacterial assembly is unambiguously ENGINEERED).
- **MEDIUM** — the element could have an alternative explanation in some genomic contexts.
- **LOW** — the element is short, degenerate, or common enough to appear by chance; weak evidence even within its tier.

---

## Per-contig summarisation (`scorer.summarize`)

After scanning, hits are grouped by contig and aggregated into `ContigSummary` objects.

### ContigSummary fields

| Field | Description |
| :--- | :--- |
| `contig` | Contig name |
| `contig_length` | Length in bp |
| `total_hits` | Total hits (mappy + kmer, before overlap merging) |
| `engineered_hits` | Hits with `tier == "ENGINEERED"` |
| `context_dependent_hits` | Hits with `tier == "CONTEXT_DEPENDENT"` |
| `weak_hits` | Hits with `tier == "WEAK"` |
| `unannotated_hits` | Hits with no metadata entry |
| `unique_feature_types` | Distinct `feature_type` values across all hits |
| `unique_labels` | Distinct `label` values across all hits |
| `covered_bp` | Merged non-overlapping bp covered by hits |
| `covered_fraction` | `covered_bp / contig_length` |
| `top_label` | Label of the hit with the highest `identity × coverage` score |
| `top_tier` | Tier of the top hit |
| `top_confidence` | Confidence of the top hit |
| `evidence_summary` | Compact semicolon-separated string, e.g. `ENG[3]: CMV enhancer, HIV-1 Psi; CD[1]: ColE1 ori` |

### Overlap merging

`covered_bp` is computed by merging all hit intervals on a contig using a sort-and-sweep algorithm, so overlapping or adjacent hits are counted only once. The `reporter` sums pre-computed `covered_bp` values across contigs rather than re-merging all intervals at the sample level.

### Contig sort order

Summaries are returned sorted by `(-engineered_hits, -context_dependent_hits, contig)` — the most strongly contaminated contigs appear first.

---

## Sample-level reporting (`reporter.build_sample_row`)

After per-contig summaries are produced, `reporter.build_sample_row` aggregates them into a single dict representing one row of `vectrap_summary.tsv`.

### Sample verdict logic

| Condition | Verdict |
| :--- | :--- |
| `engineered_hits >= --min-engineered-contamination` (default 3) | `CONTAMINATION` |
| `engineered_hits > 0` OR `context_dependent_hits >= --min-context-suspected` (default 1) | `SUSPECTED` |
| Neither of the above | `CLEAN` |
| Unhandled exception during processing | `ERROR` |

The `CONTAMINATION` verdict indicates strong, multi-hit evidence of synthetic contamination. `SUSPECTED` flags samples with weaker or ambiguous evidence for manual review. Thresholds are configurable via CLI flags.

---

## Per-sample error isolation

Each sample is processed inside an isolated try/except block. On failure:

1. Any partially written `_hits.tsv` or `_contig_verdicts.tsv` are deleted.
2. A one-line warning is printed to stderr: `[HH:MM:SS] [i/n] <sample>  ERROR (<ExcType>: <msg>)  elapsed=Xs`.
3. With `-v`, the full traceback is printed.
4. An `ERROR` row is inserted in `vectrap_summary.tsv` with `sample_verdict=ERROR` and the error message in `top_labels`.
5. The remaining samples continue normally.

A final warning is printed if any samples failed: `WARNING: k/n sample(s) failed — see ERROR rows in summary.`

---

## Verbose logging

When `-v` is passed, VecTrap prints timestamped progress messages to stderr at each stage of each sample run:

```
[14:31:02] Loading catalog indexes ...
[14:31:04] Indexes ready  elapsed=2.1s
[14:31:04] [1/3] sample1  hits=42  contigs_with_hits=3  verdict=CONTAMINATION  elapsed=1.8s
[14:31:05] [2/3] sample2  ERROR (UnicodeDecodeError: ...)  elapsed=0.1s
[14:31:07] [3/3] sample3  hits=0   contigs_with_hits=0  verdict=CLEAN  elapsed=2.0s
[14:31:07] Summary written to vectrap_summary.tsv
[14:31:07] WARNING: 1/3 sample(s) failed — see ERROR rows in summary.
[14:31:07] Done.
```

Per-stage verbose output also shows input path, catalog path, threshold settings, stage boundaries, hit counts, and elapsed time per stage.

---

## Performance characteristics

### Index loading

The mappy index (`combined_long.mmi`) and the Aho-Corasick automaton are loaded/built **once** per batch run and reused across all samples. For a batch of 100 samples this avoids 100 redundant index loads and automaton builds.

The catalog metadata dict (`catalog_metadata.pkl`) is likewise loaded once and passed to every `scan()` call.

### FASTA reading

Each assembly FASTA file is read in a **single pass** that simultaneously builds the contig sequence list and the contig-length dict. Both the mappy aligner and the k-mer scanner receive the pre-loaded contig list; the file is not opened again.

### Long-sequence scanning

Performance depends on assembly size, catalog size, and mappy thread count (`-t`, default 4). The mappy thread count controls the number of minimap2 alignment threads per sample.

### Short-sequence scanning

The Aho-Corasick automaton runs in O(contig length + number of matches) per contig after O(total pattern length) build time. This makes the short-motif phase fast regardless of catalog size — it is effectively linear in total assembly length.

### `covered_bp` computation

Overlap merging is done once per contig inside `scorer.summarize()`. The sample-level `covered_bp` in `vectrap_summary.tsv` is a simple sum over `ContigSummary.covered_bp` values and does not re-merge intervals.

---

## Design choices

### Unified coordinate system

All hits are projected onto forward-strand contig coordinates. This avoids error-prone strand-specific downstream merging and makes interval logic uniform across mappy and k-mer hits.

### Catalog-relative coverage

Coverage is measured on the catalog element, not the query contig. This answers: how much of the known synthetic element is present? A 90% covered hit is strong evidence even if the hit spans only a small fraction of a large contig.

### Scanning and scoring are decoupled

`homology_scanner.scan()` returns raw `HomologyHit` objects; `scorer.summarize()` aggregates them. Neither module knows about output formats. `reporter` and `cli/run.py` handle I/O. This makes each layer independently testable and replaceable.

### Tier and confidence come from the catalog, not computed heuristics

`tier`, `confidence`, and `reasoning` are static metadata fields attached to each catalog entry in `catalog_manifest.tsv` and compiled into `catalog_metadata.pkl` at build time. The scanner does not infer these fields — it looks them up. This keeps the evidence model explicit and auditable.

### Catalog-driven detection

VecTrap does not attempt ab initio prediction of synthetic origin. It is intentionally based on explicit sequence evidence from a curated catalog. This limits false positives at the cost of missing novel elements not yet in the catalog.

---

## Limitations

- Detection quality depends on catalog completeness. Novel or highly diverged synthetic elements not represented in the catalog will be missed.
- CONTEXT_DEPENDENT hits from a single functional class should not be taken as conclusive evidence. A confident synthetic verdict from this tier requires functional class diversity.
- WEAK hits are by definition low-specificity and should never be interpreted in isolation.
- The tier boundary between CONTEXT_DEPENDENT and ENGINEERED reflects current biological knowledge and is subject to revision as per-feature metadata is reviewed.
- Final interpretation should consider the biological context, the host organism being assembled, assembly quality, and the full pattern of evidence across the contig.

---

## Module layout

```
vectrap/cli/build_db.py         database preparation entry point
vectrap/cli/run.py              main pipeline entry point; batch loop, error isolation, TSV output
vectrap/modules/homology_scanner.py
                                long/short sequence scanning; HomologyHit dataclass;
                                load_mappy_index(), build_kmer_automaton(), scan()
vectrap/modules/scorer.py       hit aggregation per contig; ContigSummary dataclass; summarize()
vectrap/modules/reporter.py     sample-level aggregation and TSV writing; build_sample_row(), write_summary()
vectrap/modules/utils.py        rev_comp() and other helpers
```

---

## Example execution flow (single sample)

1. `load_mappy_index()` — open `combined_long.mmi` with mappy (once per session).
2. `build_kmer_automaton()` — load `short_index.pkl` and build automaton (once per session).
3. `_load_metadata()` — load `catalog_metadata.pkl` (once per session).
4. `scan()` — single-pass FASTA read → `contigs` list + `contig_lengths` dict.
5. `scan()` — mappy alignment of all contigs against the long index.
6. `scan()` — Aho-Corasick scan of all contigs for short patterns.
7. `scan()` — metadata enrichment: fill `tier`, `confidence`, `reasoning`, `label` on every hit.
8. `_write_hits_tsv()` — write `<sample>_hits.tsv`.
9. `scorer.summarize()` — aggregate hits into `ContigSummary` list.
10. `_write_verdicts_tsv()` — write `<sample>_contig_verdicts.tsv`.
11. `reporter.build_sample_row()` — produce one summary dict.
12. (After all samples) `reporter.write_summary()` — write `vectrap_summary.tsv`.

---

## Future extension points

- Additional evidence types or catalog categories
- Refined scoring weights per functional class
- Species- or context-aware false-positive suppression for CONTEXT_DEPENDENT hits
- Richer HTML or JSON reports
- Per-contig classification verdict (currently only the hit/summary tables are written; a `CONTAMINATION`/`SUSPECTED`/`CLEAN` label per contig is a natural addition)
- Multiprocessing across samples for very large batches

---

## Related files

- `README.md` — user-facing overview, quick start, CLI reference, output format tables
- `vectrap/catalogs/README.md` — full catalog tier reference with per-category rationale
- `vectrap/modules/homology_scanner.py` — scanning implementation
- `vectrap/modules/scorer.py` — per-contig summarisation
- `vectrap/modules/reporter.py` — sample-level aggregation and TSV writing

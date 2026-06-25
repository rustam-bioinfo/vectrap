# VecTrap INFO

This file is a technical reference for VecTrap. It explains what the tool does, how the pipeline works internally, what algorithms are used, what data structures are involved, and how to interpret the output.

## Purpose

VecTrap detects synthetic laboratory vector contamination in assembled nucleotide sequences, with the current focus on bacterial genome assemblies and plasmid datasets. The core idea is simple: scan each input contig against a curated catalog of known vector-associated sequence elements, then aggregate the resulting evidence into a contig-level verdict.

VecTrap is designed for cases where a deposited or locally assembled sequence may contain engineered backbone fragments such as replication origins, selectable markers, recombination sites, promoters, primer binding sites, terminators, or other laboratory-use elements.

## High-level workflow

The pipeline has two major phases:

1. Homology scanning
2. Evidence scoring and contig classification

At a high level:

- Long catalog sequences are searched with mappy/minimap2 indexing.
- Short catalog sequences are searched with an Aho-Corasick multi-pattern matcher.
- All raw hits are converted into a unified internal representation.
- Hits are then scored according to evidence strength.
- Each contig receives a final verdict such as `VECTOR`, `SUSPECTED`, or `CLEAN`.

## Inputs

VecTrap requires:

- An input FASTA file containing one or more contigs.
- A prepared sequence catalog directory.

The catalog directory is produced by `vectrap-build-db` and contains at least:

- `combined_long.mmi` — minimap2/mappy index for longer catalog entries.
- `short_index.pkl` — serialized dictionary of short catalog sequences.

## Why two scanners are used

VecTrap uses two different search strategies because short exact motifs and longer homologous segments behave very differently computationally.

### Long sequences

Long catalog entries are best handled by an aligner that tolerates mismatches and partial divergence. For this, VecTrap uses **mappy**, the Python binding to **minimap2**. This is appropriate for elements such as longer origins, markers, or regulatory modules where exact matching would be too strict.

### Short sequences

Short elements are a poor fit for standard long-sequence alignment because they are abundant, exact matching is usually sufficient, and naive repeated substring search becomes slow. For these, VecTrap uses an **Aho-Corasick automaton**, which can search for all short patterns simultaneously in a single pass through each contig.

This split design gives both sensitivity and speed.

## Algorithms used

### 1. mappy / minimap2 for long catalog entries

For catalog entries at or above the internal length threshold, VecTrap queries the input contigs against a pre-built `.mmi` index using `mappy.Aligner`.

Conceptually:

- minimap2 builds an efficient seed-and-chain index over the long catalog entries.
- Each contig is aligned against the indexed catalog.
- Candidate alignments are filtered by user-configurable thresholds.

Current filtering includes:

- Minimum sequence identity (`min_identity`)
- Minimum coverage of the catalog entry (`min_coverage`)

Identity is computed from the alignment object as:

- `identity = mlen / blen`

where:

- `mlen` = number of matching bases
- `blen` = alignment block length

Coverage is computed relative to the matched catalog sequence as:

- `coverage = (r_en - r_st) / ctg_len`

where:

- `r_st`, `r_en` = start and end on the catalog sequence
- `ctg_len` = total length of the catalog entry

This design answers the biologically useful question: how much of the known synthetic element is covered by the hit?

### 2. Aho-Corasick automaton for short catalog entries

For short catalog entries, VecTrap uses the **Aho-Corasick** string matching algorithm.

Why this matters:

- A naive approach would loop over every k-mer and search every contig separately.
- That scales roughly with number of patterns × contig length.
- Aho-Corasick instead builds a finite-state automaton once, then scans each contig in linear time.

Conceptually:

1. All short catalog patterns are inserted into a trie.
2. Failure links are added to turn the trie into an automaton.
3. Each contig is streamed through the automaton character by character.
4. Every exact match to any stored pattern is emitted immediately.

This lets VecTrap find thousands of short motifs in one pass per contig rather than one pass per motif.

## Strand handling

VecTrap scans both orientations.

For the short-pattern scanner, both the forward pattern and its reverse complement are inserted into the automaton. This means a single pass over the forward contig sequence can still detect matches on either strand.

For long-sequence hits, strand is obtained from the aligner output.

### Coordinate convention

All reported coordinates are:

- 0-based
- half-open (`start`, `end`)
- mapped onto the **forward strand of the original contig**

This is true even if the match biologically occurs on the reverse-complement strand.

That convention makes downstream interval logic much easier because all hits live in the same coordinate system.

## Internal hit representation

Both scanners return results as a shared `HomologyHit` object with these fields:

- `contig` — contig name
- `start` — 0-based start
- `end` — 0-based exclusive end
- `strand` — `+` or `-`
- `identity` — fraction in `[0, 1]`
- `coverage` — fraction in `[0, 1]`
- `catalog_id` — identifier of the matched catalog entry
- `source` — scanner origin, currently `mappy` or `kmer`

For exact short-pattern matches:

- `identity = 1.0`
- `coverage = 1.0`

A convenience property `length = end - start` is also available.

## Length threshold split

VecTrap internally separates catalog entries into:

- short sequences: length < 50 bp
- long sequences: length >= 50 bp

This threshold is chosen so that exact multi-pattern search handles the short motif-like entries efficiently, while longer biologically meaningful sequences are handled by alignment.

This same split is used during database build and scan time, so both components must stay synchronized.

## Catalog indexing

The `vectrap-build-db` command prepares the catalog in two forms.

### Long catalog index

Long entries are concatenated into a FASTA and indexed into `combined_long.mmi` using minimap2-compatible indexing through mappy/minimap2 infrastructure.

### Short catalog index

Short entries are stored in `short_index.pkl` as a Python dictionary roughly of the form:

```python
{
    "ACTG...": ["catalog_id_1", "catalog_id_2"],
    "TTGA...": ["catalog_id_3"]
}
```

This representation allows one short sequence to correspond to one or more catalog identifiers.

At runtime, this dictionary is loaded and transformed into the Aho-Corasick automaton.

## Evidence model

VecTrap does not treat all detected elements equally. The design distinguishes between stronger and weaker indicators of synthetic origin.

### PRIMARY evidence

These are elements considered highly specific to engineered constructs, such as:

- replication origins
- transfer origins
- recombination sites
- synthetic ribosome binding sites

A single strong PRIMARY hit may be sufficient to classify a contig as vector-derived.

### SUPPORTIVE evidence

These are elements associated with vectors but not always unique to them, such as:

- selectable markers
- promoters
- terminators
- primer binding sites
- regulatory motifs
- retroviral elements
- mobile-element-derived fragments

SUPPORTIVE hits are aggregated into a weighted score and interpreted in context.

## Contig-level classification

After scanning, hits are aggregated per contig.

While the exact logic lives in `scorer.py`, the intended design is:

- `VECTOR` — strong direct evidence of synthetic vector origin
- `SUSPECTED` — supportive but not definitive evidence
- `CLEAN` — no convincing vector evidence detected

This separation is important because some motifs commonly found in vectors can also occur in natural biological sequences.

## Verbose logging

When `-v` is used, VecTrap prints timestamped progress messages to stderr.

The log is intended to show:

- input path
- catalog path
- threshold settings
- stage boundaries
- hit counts
- elapsed time per stage
- total runtime

This is useful for large assemblies and for diagnosing slow stages or unexpected hit volumes.

## Performance notes

### Long-sequence performance

Long-sequence scanning performance depends mostly on:

- number of contigs
- total assembly size
- catalog size
- minimap2/mappy threading

### Short-sequence performance

The short-sequence scanner was intentionally implemented with Aho-Corasick instead of repeated `str.find()` loops.

This matters because:

- repeated substring search scales poorly with many patterns
- Aho-Corasick scans each contig once
- runtime becomes close to linear in total input length, after automaton construction

In practice, this makes the short-motif phase fast enough even when the catalog contains thousands of entries.

## Output interpretation

Typical raw hit output contains:

- where in a contig a feature was found
- which strand it matched on
- which catalog entry it matched
- which scanner found it
- identity and coverage for longer hits

These hits are then summarized into a per-contig verdict table.

A contig with multiple independent vector-associated features is much more suspicious than a contig with only one weak generic promoter-like motif.

## Design choices

A few implementation choices are deliberate:

### Unified coordinate system

All hits are projected onto forward-strand contig coordinates. This avoids error-prone strand-specific downstream merging.

### Catalog-relative coverage

Coverage is measured on the catalog element, not the query contig. This is more meaningful because the question is whether a known vector element is substantially represented.

### Separate scanning and scoring

Homology detection and biological interpretation are decoupled. This makes the pipeline easier to test, benchmark, and extend.

### Catalog-driven detection

VecTrap does not attempt ab initio prediction of synthetic origin. It is intentionally based on explicit sequence evidence from curated synthetic elements.

## Limitations

Users should be aware of several limitations:

- Detection quality depends on catalog completeness.
- Very short motifs can occur by chance, so context matters.
- Highly diverged or novel synthetic elements may be missed if they are insufficiently represented in the catalog.
- Naturally occurring sequences can resemble some supportive synthetic motifs.
- Final interpretation should consider biology, assembly quality, and neighboring evidence.

## Current module layout

The repository is organized around a small modular architecture:

- `vectrap/cli/build_db.py` — database preparation
- `vectrap/cli/run.py` — main command-line entry point
- `vectrap/modules/homology_scanner.py` — long/short sequence scanning
- `vectrap/modules/scorer.py` — evidence aggregation and contig classification
- `vectrap/modules/utils.py` — FASTA parsing and helper functions

## Example execution flow

A typical run looks like this:

1. Load the query FASTA.
2. Open `combined_long.mmi` with mappy.
3. Load `short_index.pkl`.
4. Build the Aho-Corasick automaton from all short patterns.
5. Scan all contigs with mappy for long-sequence hits.
6. Scan all contigs with the automaton for short exact hits.
7. Merge all hits into a single list.
8. Score evidence per contig.
9. Write hit-level and contig-level output tables.

## Future extension points

The current design allows several natural extensions:

- additional evidence types or catalog categories
- alternative scoring models
- interval merging of overlapping nearby hits
- reporting of clustered vector architectures
- species- or context-aware false-positive suppression
- richer HTML or JSON reports

## Suggested citation

If VecTrap is used in research, cite the software repository and the Zenodo-hosted sequence catalog described in the main `README.md`.

## Related files

- `README.md` — user-facing overview and quick start
- `vectrap/modules/homology_scanner.py` — technical scanning implementation
- `vectrap/modules/scorer.py` — classification logic
- `vectrap/catalogs/README.md` — catalog acquisition and preparation

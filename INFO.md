# VecTrap INFO

This file is a technical reference for VecTrap. It explains what the tool does, how the pipeline works internally, what algorithms are used, what data structures are involved, and how to interpret the output.

## Purpose

VecTrap detects sequences of synthetic origin in assembled nucleotide sequences, with the current focus on bacterial genome assemblies and plasmid datasets. The core question it answers is: does any part of this assembled sequence derive from an engineered laboratory construct rather than from the biological sample under study?

VecTrap is designed for cases where a deposited or locally assembled sequence may contain engineered backbone fragments — replication origins, selectable markers, recombination sites, synthetic promoters, primer binding sites, terminators, or other laboratory-use elements. It does not rely on database annotations or assembly metadata; it works directly from sequence evidence.

## High-level workflow

The pipeline has two major phases:

1. Homology scanning
2. Evidence scoring and contig classification

At a high level:

- Long catalog sequences are searched with mappy/minimap2 indexing.
- Short catalog sequences are searched with an Aho-Corasick multi-pattern matcher.
- All raw hits are converted into a unified internal representation.
- Hits are scored according to their evidence tier and functional class diversity.
- Each contig receives a final verdict: `VECTOR`, `SUSPECTED`, or `CLEAN`.

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

Long catalog entries are best handled by an aligner that tolerates mismatches and partial divergence. For this, VecTrap uses **mappy**, the Python binding to **minimap2**. This is appropriate for elements such as longer origins, resistance markers, or regulatory modules where exact matching would be too strict due to natural sequence variation across strains and vector generations.

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

This answers the biologically useful question: how much of the known synthetic element is covered by the hit?

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

VecTrap does not treat all detected elements equally. All tier assignments are made relative to the **bacterial genome context** — the defining question is whether a given sequence could plausibly appear in a wild-type bacterial chromosome or a natural bacterial plasmid without deliberate cloning.

Every catalog entry is assigned to one of three evidence tiers.

### ENGINEERED tier

Sequences in this tier have no natural occurrence in any wild-type bacterium. They originate from eukaryotic viruses, mammalian cell biology, phage RNA polymerase systems, or are fully synthetic constructs with no plausible route into a bacterial genome except deliberate cloning. A single high-confidence ENGINEERED hit is sufficient to classify a contig as synthetic.

This tier includes:

- **Eukaryotic and viral promoters**: CMV, EF1α, PGK, CAG, CaMV35S. These are active only in eukaryotic cells and have no function in bacteria; their presence in a bacterial assembly is unambiguous.
- **Phage RNA polymerase promoters**: T7, SP6, T3. These require T7/SP6/T3 RNA polymerase, which is never endogenous in wild-type bacteria. Their use is exclusively in engineered expression systems.
- **Eukaryotic and viral replication origins**: SV40 ori, EBV oriP. Never functional in bacteria.
- **Filamentous phage origins used as cloning tools**: f1 ori, M13 ori. These are phage-derived elements used in phagemid systems; they do not occur in bacterial chromosomes or natural conjugative plasmids.
- **Eukaryotic transcription terminators and polyadenylation signals**: BGH polyA, SV40 polyA. Eukaryotic-specific; have no role in bacterial transcription.
- **Mammalian selectable markers**: Puromycin resistance (pac), Blasticidin resistance (bsr), Hygromycin resistance (hph) when used in mammalian vector context. These are sometimes found in bacterial vectors for propagation, but their primary design context is mammalian cell selection.
- **Fluorescent and reporter proteins**: EGFP, mCherry, mTurquoise, mVenus, mRuby, tdTomato, Luciferase, and all GFP-family variants. Never encoded in wild-type bacterial genomes.
- **CRISPR and genome-editing elements**: Cas9, Cas12, Cpf1, sgRNA scaffolds. Engineered systems derived from bacterial immunity components but assembled into configurations that do not occur naturally.
- **Recombinase systems**: Cre recombinase, Flp recombinase, PhiC31 integrase. Derived from phage but used exclusively as engineered tools.
- **Site-specific recombination sites**: attB, attP, loxP, FRT, and designer variants. Distinguished from naturally occurring phage att sites by the synthetic flanking context in which they appear.
- **Synthetic regulatory elements**: IRES sequences, WPRE, P2A/T2A/E2A self-cleaving peptides, synthetic Kozak sequences, synthetic splice donor/acceptor sites. Eukaryotic-specific or fully synthetic.
- **Sequencing primer binding sites**: M13, T7, SP6, pUC, T3 primer sites. Designed for laboratory use; no natural occurrence.
- **Retroviral and lentiviral elements**: LTR sequences, lentiviral packaging signals (PSI), HIV-derived elements. Viral in origin; never present in bacterial genomes outside deliberate cloning.

### CONTEXT_DEPENDENT tier

Sequences in this tier exist naturally in bacteria or natural bacterial plasmids. Their presence in an assembly is not automatically diagnostic — it depends on which organism is being assembled and what the broader genomic context looks like.

For example, a ColE1 origin found in an *E. coli* assembly could reflect a natural endogenous plasmid. The same origin found in a *Salmonella enterica* or *Bacillus subtilis* assembly is anomalous and warrants scrutiny. Antibiotic resistance genes such as KanR and AmpR are endemic in environmental and clinical bacterial populations via horizontal gene transfer and cannot individually prove synthetic contamination.

Their evidential weight comes from **functional class co-occurrence**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter from this tier is almost certainly synthetic, whereas a contig with only a single hit from one functional class is ambiguous.

This tier includes:

- **Bacterial replication origins**: ColE1, pMB1, pUC ori, p15A, pSC101, RK2, pBBR1, pSa. Derived from natural *E. coli* and broad-host-range plasmids.
- **Transfer origins (oriT)**: IncP, IncQ, IncW oriT sequences. Found in natural conjugative plasmids.
- **Bacterial antibiotic resistance markers**: aphA (KanR), bla (AmpR), cat (CmR), tetA (TetR), aac(3) (GenR), aadA (SpecR). Originally derived from transposons; widely disseminated in natural populations via HGT.
- **T7 expression system elements**: T7 terminator, T7 tag. The T7 promoter itself is ENGINEERED; these associated elements are CONTEXT_DEPENDENT because the terminator sequence occurs in some natural phage genomes and the tag is a protein epitope.
- **Inducible expression system components**: lacI, lacO arrays, tac promoter, trc promoter, IPTG-responsive elements. Derived from the natural *E. coli* lac operon but assembled into configurations not found in wild-type chromosomes.
- **Mobile element remnants**: Transposon end sequences (Tn3, Tn5, Tn10, Tn903) and IS-element borders incorporated into vector backbones during historical cloning from natural sources.

### WEAK tier

Sequences in this tier are widespread in wild-type bacterial genomes and carry low specificity for synthetic origin on their own. A WEAK hit contributes to a verdict only when it co-occurs with hits from ENGINEERED or CONTEXT_DEPENDENT tiers.

This tier includes:

- **Common transcription terminators**: rrnB T1/T2, tonB terminator, lambda t0, fd terminator. These are standard terminators used in expression vectors precisely because they are efficient in many bacteria — and therefore also common in natural genomes.
- **Natural sigma-factor-dependent promoters**: lac, ara, trp, phoA promoters in their unmodified natural form. Widespread in Enterobacteriaceae and related organisms.
- **Ribosomal RNA and tRNA elements**: rRNA genes, tRNA genes. Present in all bacteria.

## Contig-level classification

After scanning, hits are aggregated per contig by the scorer.

The classification logic is:

- `VECTOR` — one or more ENGINEERED hits with sufficient identity and coverage, OR CONTEXT_DEPENDENT hits spanning two or more independent functional classes.
- `SUSPECTED` — CONTEXT_DEPENDENT hits present from only one functional class, OR low-confidence ENGINEERED hits below identity/coverage thresholds, OR WEAK hits co-occurring with a single CONTEXT_DEPENDENT hit.
- `CLEAN` — no convincing evidence of synthetic origin detected.

## Verbose logging

When `-v` is used, VecTrap prints timestamped progress messages to stderr.

The log shows:

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

The short-sequence scanner uses Aho-Corasick instead of repeated `str.find()` loops.

This matters because:

- repeated substring search scales poorly with many patterns
- Aho-Corasick scans each contig once regardless of how many patterns are in the catalog
- runtime becomes close to linear in total input length after automaton construction

In practice this makes the short-motif phase fast even when the catalog contains thousands of entries.

## Output interpretation

Typical raw hit output contains:

- where in a contig a feature was found
- which strand it matched on
- which catalog entry it matched
- which scanner found it
- identity and coverage for longer hits

These hits are summarized into a per-contig verdict table.

A single ENGINEERED hit is strong evidence on its own. For CONTEXT_DEPENDENT hits, the signal is co-occurrence of hits from multiple independent functional classes on the same contig — not the raw number of hits from any one class. WEAK hits alone are never sufficient for a `VECTOR` or `SUSPECTED` verdict.

## Design choices

A few implementation choices are deliberate:

### Unified coordinate system

All hits are projected onto forward-strand contig coordinates. This avoids error-prone strand-specific downstream merging.

### Catalog-relative coverage

Coverage is measured on the catalog element, not the query contig. This answers the biologically useful question: how much of the known synthetic element is present in the contig?

### Separate scanning and scoring

Homology detection and biological interpretation are decoupled. This makes the pipeline easier to test, benchmark, and extend independently.

### Catalog-driven detection

VecTrap does not attempt ab initio prediction of synthetic origin. It is intentionally based on explicit sequence evidence from curated synthetic and recruited elements.

## Limitations

Users should be aware of several limitations:

- Detection quality depends on catalog completeness. Novel or highly diverged synthetic elements not represented in the catalog will be missed.
- CONTEXT_DEPENDENT hits from a single functional class should not be taken as conclusive evidence. The scorer requires functional class diversity for a confident verdict.
- WEAK hits are by definition low-specificity. They should never be interpreted in isolation.
- The boundary between tiers reflects current biological knowledge and is subject to revision. Some elements currently assigned CONTEXT_DEPENDENT may be promoted to ENGINEERED as per-feature metadata is refined through systematic review.
- Final interpretation should consider the biological context, the host organism being assembled, assembly quality, and the full pattern of evidence across the contig.

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
8. Score evidence per contig using tier and functional class diversity.
9. Write hit-level and contig-level output tables.

## Future extension points

The current design allows several natural extensions:

- additional evidence types or catalog categories
- refined scoring weights per functional class
- interval merging of overlapping nearby hits
- reporting of clustered vector architectures
- species- or context-aware false-positive suppression for CONTEXT_DEPENDENT hits
- richer HTML or JSON reports

## Suggested citation

If VecTrap is used in research, cite the software repository and the Zenodo-hosted sequence catalog described in the main `README.md`.

## Related files

- `README.md` — user-facing overview and quick start
- `vectrap/catalogs/README.md` — full catalog tier reference with per-category rationale
- `vectrap/modules/homology_scanner.py` — technical scanning implementation
- `vectrap/modules/scorer.py` — classification logic

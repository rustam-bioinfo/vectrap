# VecTrap Sequence Catalogs

This directory holds the VecTrap sequence catalog files and search indexes.

Catalog files are **not bundled with the repository** due to their size. They must be downloaded separately from Zenodo before running the pipeline.

## Download

```bash
vectrap-build-db --download
```

This command downloads all catalog FASTA files from Zenodo, verifies their checksums, and builds the minimap2 index and k-mer hash index in this directory.

## Zenodo Deposit

DOI: [10.5281/zenodo.20844271](https://doi.org/10.5281/zenodo.20844271)

## Files After Download

| File | Description |
| :--- | :--- |
| `catalog_*.fasta.gz` | Compressed FASTA files, one per feature type category |
| `catalog_manifest.tsv` | File list with sequence counts and MD5 checksums |
| `combined_long.fasta` | All sequences >= 50 bp, input for minimap2 indexing |
| `combined_long.mmi` | minimap2 index for long-feature scanning |
| `short_index.pkl` | Exact k-mer hash index for short-feature scanning |

## Evidence Tiers

VecTrap assigns every catalog category to one of two evidence tiers: **ENGINEERED** or
**RECRUITED**. The distinction is based on whether the sequences in a category have any
plausible natural occurrence, not simply on how commonly they appear in vector maps.

### ENGINEERED

Sequences in this tier have **no natural biological origin**. They were designed in a
laboratory, do not exist in wild-type organisms, and their presence in an assembled
sequence is essentially unambiguous evidence of synthetic origin. A single high-confidence
ENGINEERED hit is sufficient to classify a contig as synthetic.

This tier includes:

- Synthetic ribosome binding sites (BBa-series, consensus Shine-Dalgarno variants optimised
  for expression, and other designer RBS sequences that differ substantially from any
  chromosomal RBS).
- Fully synthetic promoters with no direct natural counterpart: T7 class, T5, tac, trc,
  lac-T7 hybrids, CMV, EF1α, PGK, CAG, CaMV 35S, and similar.
- Synthetic transcription terminators: T7Te, designed rho-independent terminators,
  BGH polyA, SV40 polyA, and similar.
- Site-specific recombination sites introduced by engineering: attB, attP, attL, attR
  (Lambda and PhiC31), loxP, loxN, FLP/FRT, rox, and other designer variants. These are
  distinguished from naturally occurring phage att sites by their flanking sequence context
  stored in the catalog.
- Sequencing primer binding sites in canonical orientations: M13 forward/reverse, T7, SP6,
  T3, pUC, and similar universal sequencing primers. These short exact sequences were
  designed for laboratory use and have no meaningful natural occurrence at the frequencies
  seen in cloning vectors.
- Retroviral LTR sequences used in lentiviral and gamma-retroviral expression and
  packaging vectors. These are derived from viral genomes but appear in assembled bacterial
  or eukaryotic sequences only through deliberate cloning.
- Synthetic regulatory elements: IRES sequences, Kozak consensus inserts, synthetic splice
  donor and acceptor sites, tetO and cumate operator arrays used in inducible systems.

### RECRUITED

Sequences in this tier are **derived from natural biological sources** and are commonly
recruited into synthetic vectors precisely because they function well in engineered
constructs. They exist in natural genomes as well and cannot individually serve as proof
of synthetic origin. Their evidential weight depends on co-occurrence with other hits and
on genomic context.

Multiple independent RECRUITED hits from **different functional classes** on the same
contig strongly support a synthetic origin conclusion, even when no ENGINEERED hit is
present. The key signal is functional diversity: a natural sequence may contain a
terminator-like hairpin, but it will not simultaneously carry a resistance gene, a
replication origin, and a promoter.

This tier includes:

- Replication origins: ColE1, p15A, pBBR1, RK2, pSC101, f1, and similar. These are
  derived from natural plasmids and bacteriophage and occur in wild-type extrachromosomal
  elements. Engineered variants (rop-deleted ColE1, mutated copy-number alleles) are
  included but the category as a whole cannot be called ENGINEERED because the underlying
  sequences are naturally occurring.
- Transfer origins (oriT): IncP, IncQ, IncW, and similar conjugative transfer sequences.
  These are found in natural conjugative plasmids.
- Selectable marker coding sequences: resistance genes such as aphA (KanR), bla (AmpR),
  cat (CmR), aac(3) (GentR), hph (HygR), pac (PuroR), bsr (BlastR), and similar. These
  are derived from natural bacterial resistance determinants, many from transposons, and
  occur in wild-type organisms. Their presence in a contig is suspicious but not
  individually conclusive without additional context.
- Engineered promoters derived from natural sigma-factor-dependent sequences: lac, ara,
  trp, phoA, tet, and similar. The fully synthetic core promoters listed under ENGINEERED
  are distinct; this sub-category covers natural promoter sequences that appear in vector
  maps with minimal or no modification.
- Transcription terminators derived from natural sequences: rrnB T1/T2, lambda t0, fd,
  and similar. These occur in wild-type genomes.
- Protein-binding operator sequences: lacO, tetO arrays, and similar, when cataloged in
  their natural-length form rather than as synthetic multimers.
- Mobile-element remnants: transposon end sequences (Tn3, Tn5, Tn10, Tn903) and IS-element
  borders that appear in vector backbones due to historical cloning from natural sources.
  These also occur in natural chromosomes and plasmids.
- Miscellaneous synthetic features that do not fit the above categories cleanly but share
  the property of having a recognisable natural sequence ancestor.

## Catalog Categories

| Catalog file | Feature type | Evidence tier | Rationale |
| :--- | :--- | :--- | :--- |
| `catalog_rep_origin.fasta.gz` | Replication origins | RECRUITED | Derived from natural plasmid/phage replicons; engineered variants included but natural counterparts exist |
| `catalog_oriT.fasta.gz` | Transfer origins | RECRUITED | Found in natural conjugative plasmids |
| `catalog_misc_recomb.fasta.gz` | Recombination sites | ENGINEERED | Designed att/lox/FRT sites; distinguished from phage att sites by flanking context |
| `catalog_RBS.fasta.gz` | Ribosome binding sites | ENGINEERED | Synthetic designer RBS sequences with no natural chromosomal counterpart |
| `catalog_CDS.fasta.gz` | Resistance markers and tags | RECRUITED | Derived from natural resistance determinants; suspicious in context but not individually diagnostic |
| `catalog_promoter.fasta.gz` | Promoters | ENGINEERED / RECRUITED | Fully synthetic promoters (T7, CMV, etc.) are ENGINEERED; natural sigma-factor promoters in vectors are RECRUITED |
| `catalog_terminator.fasta.gz` | Transcription terminators | ENGINEERED / RECRUITED | Synthetic terminators (BGH, T7Te) are ENGINEERED; natural terminators (rrnB) are RECRUITED |
| `catalog_primer_bind.fasta.gz` | Sequencing primer binding sites | ENGINEERED | Designed for laboratory use; no meaningful natural occurrence |
| `catalog_protein_bind.fasta.gz` | Protein binding sites (operators) | RECRUITED | Derived from natural operator sequences |
| `catalog_regulatory.fasta.gz` | Synthetic regulatory elements | ENGINEERED | IRES, synthetic Kozak, splice signals used in expression vectors |
| `catalog_enhancer.fasta.gz` | Enhancer elements | RECRUITED | Derived from viral or cellular enhancers; not specific to synthetic constructs alone |
| `catalog_LTR.fasta.gz` | Retroviral LTR sequences | ENGINEERED | Appear in bacterial/non-viral assemblies only through deliberate cloning |
| `catalog_mobile_element.fasta.gz` | Mobile element remnants | RECRUITED | Transposon and IS borders occur in natural chromosomes and plasmids |
| `catalog_misc_feature.fasta.gz` | Miscellaneous synthetic features | RECRUITED | Mixed; assessed in combination with other evidence |

## Scoring Principle

One **ENGINEERED** hit with sufficient identity and coverage is treated as direct evidence
of synthetic origin.

For **RECRUITED** hits, the scoring model weighs **functional class diversity**: a contig
that carries hits from multiple independent functional classes (replication, resistance,
expression control) is much more likely to be synthetic than a contig with many hits all
from the same functional class. Hit count alone is not a reliable signal for RECRUITED
evidence.

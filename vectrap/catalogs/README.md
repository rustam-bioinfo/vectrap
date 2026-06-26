# VecTrap Sequence Catalogs

This directory holds the VecTrap sequence catalog files and search indexes.

Catalog files are **not bundled with the repository** due to their size. They must be downloaded separately from Zenodo before running the pipeline.

---

## Download

```bash
vectrap-build-db --download
```

This command downloads all catalog FASTA files from Zenodo, verifies their checksums, builds the minimap2 index, the k-mer automaton, and the metadata pickle in this directory.

## Zenodo Deposit

DOI: [10.5281/zenodo.20844271](https://doi.org/10.5281/zenodo.20844271)

---

## Files After Build

| File | Description |
| :--- | :--- |
| `catalog_*.fasta.gz` | Compressed FASTA files, one per feature-type category |
| `catalog_manifest.tsv` | One row per catalog entry: `catalog_id`, `feature_type`, `label`, `tier`, `confidence`, `reasoning`, `length`, `md5` |
| `combined_long.fasta` | All sequences >= 50 bp, input for minimap2 indexing |
| `combined_long.mmi` | minimap2 index for long-feature scanning (sequences >= 50 bp) |
| `short_index.pkl` | Aho-Corasick k-mer index for short-feature scanning (sequences < 50 bp) |
| `catalog_metadata.pkl` | Dict `catalog_id → {feature_type, label, tier, confidence, reasoning}` — loaded at runtime to annotate hits |

`catalog_metadata.pkl` is the runtime lookup table. It is compiled from `catalog_manifest.tsv` during `vectrap-build-db` and is the source of truth for the `tier`, `confidence`, and `reasoning` fields that appear in `*_hits.tsv` output.

---

## catalog_manifest.tsv columns

| Column | Description |
| :--- | :--- |
| `catalog_id` | Unique identifier, matches the FASTA sequence header |
| `feature_type` | Functional category (e.g. `promoter`, `LTR`, `CDS`) |
| `label` | Human-readable element name (e.g. `CMV enhancer`, `HIV-1 3' LTR`) |
| `tier` | `ENGINEERED`, `CONTEXT_DEPENDENT`, or `WEAK` |
| `confidence` | `HIGH`, `MEDIUM`, or `LOW` |
| `reasoning` | One-sentence rationale for the tier assignment |
| `length` | Sequence length in bp |
| `md5` | MD5 checksum of the raw sequence |

---

## Evidence tiers

All tier assignments are made relative to the **bacterial genome context** — the defining question is whether a given sequence could plausibly appear in a wild-type bacterial chromosome or a natural bacterial plasmid without deliberate cloning.

Every catalog entry is assigned to one of three evidence tiers. Every entry also has a **confidence** level (`HIGH`, `MEDIUM`, `LOW`) that reflects how unambiguous that assignment is.

### ENGINEERED

Sequences with **no natural occurrence in any wild-type bacterium**. They originate from eukaryotic viruses, mammalian cell biology, phage RNA polymerase systems, or are fully synthetic constructs. There is no plausible route for these sequences to appear in a bacterial assembly except deliberate cloning. A single HIGH-confidence ENGINEERED hit is sufficient to flag a contig.

- **Eukaryotic and viral promoters**: CMV, EF1α, PGK, CAG, CaMV 35S. Active only in eukaryotic cells; have no function in bacteria. Their presence in a bacterial assembly is unambiguous.
- **Phage RNA polymerase promoters**: T7, SP6, T3. Require T7/SP6/T3 RNA polymerase, which is never endogenous in wild-type bacteria. Used exclusively in engineered expression systems.
- **Eukaryotic and viral replication origins**: SV40 ori, EBV oriP. Never functional in bacteria.
- **Filamentous phage origins used as cloning tools**: f1 ori, M13 ori. Used in phagemid systems; do not occur in bacterial chromosomes or natural conjugative plasmids.
- **Eukaryotic transcription terminators and polyadenylation signals**: BGH polyA, SV40 polyA. Eukaryotic-specific; have no role in bacterial transcription.
- **Mammalian selectable markers**: Puromycin resistance (*pac*), Blasticidin resistance (*bsr*), Hygromycin resistance (*hph*) in mammalian vector context.
- **Fluorescent and reporter proteins**: EGFP, mCherry, mTurquoise, mVenus, mRuby, tdTomato, Luciferase, and all GFP-family variants. Never encoded in wild-type bacterial genomes.
- **CRISPR and genome-editing elements**: Cas9, Cas12, Cpf1, sgRNA scaffolds. Assembled into configurations that do not occur naturally.
- **Recombinase systems**: Cre recombinase, Flp recombinase, PhiC31 integrase. Derived from phage but used exclusively as engineered tools.
- **Site-specific recombination sites**: attB, attP, attL, attR (Lambda and PhiC31), loxP, loxN, FRT, rox, and designer variants, distinguished from naturally occurring phage att sites by flanking sequence context stored in the catalog.
- **Synthetic regulatory elements**: IRES sequences (EMCV, HCV), WPRE, P2A/T2A/E2A self-cleaving peptides, synthetic Kozak sequences, synthetic splice donor/acceptor sites. Eukaryotic-specific or fully synthetic.
- **Sequencing primer binding sites**: M13 forward/reverse, T7, SP6, T3, pUC, and similar universal sequencing primers. Designed for laboratory use; no meaningful natural occurrence in bacteria.
- **Retroviral and lentiviral elements**: LTR sequences (HIV-1, MMLV, HTLV-1, MSCV), packaging signals (HIV-1 PSI, MMLV PSI), Rev response element (RRE), cPPT/CTS. Viral in origin; never present in bacterial genomes outside deliberate cloning.
- **Viral enhancers**: CMV immediate-early enhancer, SV40 early enhancer, baculovirus hr5 enhancer. Viral regulatory elements; no natural occurrence in bacteria.

### CONTEXT_DEPENDENT

Sequences that **exist naturally in bacteria or natural bacterial plasmids**. Their presence is not automatically diagnostic — it depends on the host organism and broader genomic context. For example, a ColE1 origin in an *E. coli* assembly may reflect a natural endogenous plasmid; the same origin in a *Salmonella* or *Bacillus* assembly is anomalous. Antibiotic resistance genes such as KanR and AmpR are endemic in environmental and clinical populations via horizontal gene transfer and cannot individually prove synthetic contamination.

Their evidential weight comes from **functional class co-occurrence**: a contig simultaneously carrying a replication origin, a resistance marker, and a promoter from this tier is almost certainly synthetic, whereas a contig with a single hit from one functional class is ambiguous.

- **Bacterial replication origins**: ColE1, pMB1, pUC ori, p15A, pSC101, RK2, pBBR1, pSa, and similar. Derived from natural *E. coli* and broad-host-range plasmids.
- **Transfer origins (oriT)**: IncP, IncQ, IncW, and similar conjugative transfer sequences. Found in natural conjugative plasmids.
- **Bacterial antibiotic resistance markers**: *aphA* (KanR), *bla* (AmpR), *cat* (CmR), *tetA* (TetR), *aac(3)* (GentR), *aadA* (SpecR), and similar. Originally derived from transposons; widely disseminated via horizontal gene transfer.
- **T7 expression system elements**: T7 terminator, T7 tag. The T7 promoter is ENGINEERED; these associated elements are CONTEXT_DEPENDENT because the terminator occurs in some natural phage genomes and the tag is a protein epitope sequence.
- **Inducible expression system components**: *lacI*, *lacO* arrays, tac promoter, trc promoter, IPTG-responsive elements. Derived from the natural *E. coli* lac operon but assembled into configurations not found in wild-type chromosomes.
- **Agrobacterium T-DNA border repeats**: RB and LB T-DNA repeats. Occur naturally on Ti plasmids in *Agrobacterium*; only diagnostic outside that native context.
- **Mobile element remnants**: Transposon end sequences (Tn3, Tn5, Tn10, Tn903) and IS-element borders incorporated into vector backbones during historical cloning. Also occur in natural chromosomes and plasmids.
- **Protein-binding operator sequences**: *lacO*, *tetO* arrays and similar, when cataloged in their natural-length form rather than as synthetic multimers.
- **pBR322 basis-of-mobility (bom)**: Derived from natural bacterial plasmid transfer sequences; suggestive mainly in a vector-like context.

### WEAK

Sequences **widespread in wild-type bacterial genomes** with low specificity for synthetic origin. A WEAK hit contributes to a verdict only when co-occurring with ENGINEERED or CONTEXT_DEPENDENT evidence. WEAK hits alone are never sufficient for a `CONTAMINATION` or `SUSPECTED` verdict.

- **Common transcription terminators**: rrnB T1/T2, tonB terminator, lambda t0, fd terminator. Used in expression vectors precisely because they are efficient in many bacteria — and therefore common in natural genomes.
- **Natural sigma-factor-dependent promoters**: lac, ara, trp, phoA in their unmodified natural form. Widespread in Enterobacteriaceae and related organisms.
- **Short regulatory motifs**: Human RNA polymerase II pause signals and similar short elements that are common enough to appear by chance in unrelated sequences.

### Confidence levels

| Confidence | Meaning |
| :--- | :--- |
| `HIGH` | Presence of this element is unambiguous evidence for its tier. A viral LTR or EGFP in a bacterial assembly is unambiguously ENGINEERED. |
| `MEDIUM` | The element could have an alternative explanation in some genomic contexts. |
| `LOW` | The element is short, degenerate, or common enough to appear by chance; weak evidence even within its tier. |

---

## Catalog categories

| Catalog file | Feature type | Default tier(s) | Rationale |
| :--- | :--- | :--- | :--- |
| `catalog_enhancer.fasta.gz` | `enhancer` | ENGINEERED | Viral and mammalian enhancers (CMV, SV40, hr5) have no natural occurrence in bacteria |
| `catalog_LTR.fasta.gz` | `LTR` | ENGINEERED | Retroviral LTR sequences appear in bacterial assemblies only through deliberate cloning |
| `catalog_primer_bind.fasta.gz` | `primer_bind` | ENGINEERED | Universal sequencing primers are designed for laboratory use; no meaningful natural occurrence |
| `catalog_regulatory.fasta.gz` | `regulatory` | ENGINEERED | IRES, WPRE, P2A/T2A, synthetic Kozak, splice signals — eukaryotic-specific or fully synthetic |
| `catalog_misc_recomb.fasta.gz` | `misc_recomb` | ENGINEERED | Designer att/lox/FRT sites; no natural occurrence in wild-type bacteria |
| `catalog_rep_origin.fasta.gz` | `rep_origin` | CONTEXT_DEPENDENT | Derived from natural plasmid replicons; natural counterparts exist in bacteria; tier reflects host-organism dependency |
| `catalog_oriT.fasta.gz` | `oriT` | CONTEXT_DEPENDENT | Found in natural conjugative plasmids; not diagnostic alone |
| `catalog_protein_bind.fasta.gz` | `protein_bind` | CONTEXT_DEPENDENT | Derived from natural operator sequences; individually ambiguous |
| `catalog_mobile_element.fasta.gz` | `mobile_element` | CONTEXT_DEPENDENT | Transposon and IS borders occur in natural chromosomes and plasmids |
| `catalog_CDS.fasta.gz` | `CDS` | ENGINEERED / CONTEXT_DEPENDENT | Fluorescent proteins, Cas9, recombinases: ENGINEERED; Tn-derived resistance markers: CONTEXT_DEPENDENT |
| `catalog_promoter.fasta.gz` | `promoter` | ENGINEERED / CONTEXT_DEPENDENT / WEAK | T7/CMV/CAG: ENGINEERED; tac/trc/lacI-regulated: CONTEXT_DEPENDENT; unmodified lac/ara/trp: WEAK |
| `catalog_terminator.fasta.gz` | `terminator` | ENGINEERED / WEAK | BGH/SV40 polyA: ENGINEERED; rrnB/tonB/lambda t0: WEAK |
| `catalog_RBS.fasta.gz` | `RBS` | CONTEXT_DEPENDENT / WEAK | Synthetic designer RBS: CONTEXT_DEPENDENT; minimal natural-sequence variants: WEAK |
| `catalog_misc_feature.fasta.gz` | `misc_feature` | ENGINEERED / CONTEXT_DEPENDENT / WEAK | Mixed; tier assigned per entry in `catalog_manifest.tsv`; assessed in combination with other evidence |

---

## Scoring principles

**ENGINEERED hits** — one HIGH-confidence ENGINEERED hit with sufficient identity and coverage is treated as direct evidence of synthetic origin. The sample-level threshold for a `CONTAMINATION` verdict is configurable (`--min-engineered-contamination`, default 3 hits). A single ENGINEERED hit triggers `SUSPECTED`.

**CONTEXT_DEPENDENT hits** — the scoring model weighs **functional class co-occurrence**. A contig that carries hits from multiple independent functional classes (replication, resistance, expression control) is much more likely to be synthetic than a contig with hits from only one functional class. Functional class diversity is tracked via `unique_feature_types` in the contig summary.

**WEAK hits** — never scored in isolation. They contribute only when ENGINEERED or CONTEXT_DEPENDENT evidence is already present on the same contig.

**Overlap merging** — when multiple hits overlap on a contig, the `covered_bp` field reports the total non-overlapping coverage (computed by interval merging). This prevents double-counting adjacent or duplicate catalog entries.

**Per-entry tier and confidence** — the tier and confidence in `*_hits.tsv` are per-entry values from `catalog_manifest.tsv`, not category-level defaults. Within mixed-tier categories (e.g. `misc_feature`, `CDS`), individual entries carry their own independently assigned tier and reasoning. The category column in the table above shows the most common tier(s) in that catalog file, but the manifest is authoritative.

---

## Adding entries to the catalog

To add new entries:

1. Add the sequence(s) to the appropriate `catalog_*.fasta.gz` file.
2. Add a corresponding row to `catalog_manifest.tsv` with all required columns (`catalog_id`, `feature_type`, `label`, `tier`, `confidence`, `reasoning`, `length`, `md5`).
3. Re-run `vectrap-build-db` (without `--download`) to rebuild `combined_long.mmi`, `short_index.pkl`, and `catalog_metadata.pkl` from the updated files.

Sequences < 50 bp will go into `short_index.pkl`; sequences >= 50 bp will go into `combined_long.mmi`. Both thresholds are applied automatically by `vectrap-build-db`.

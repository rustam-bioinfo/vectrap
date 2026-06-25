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

VecTrap assigns every catalog category to one of three evidence tiers: **ENGINEERED**,
**CONTEXT_DEPENDENT**, or **WEAK**. All tier assignments are made relative to the
**bacterial genome context** — the defining question is whether a given sequence could
plausibly appear in a wild-type bacterial chromosome or a natural bacterial plasmid
without deliberate cloning.

### ENGINEERED

Sequences in this tier have **no natural occurrence in any wild-type bacterium**. They
originate from eukaryotic viruses, mammalian cell biology, phage RNA polymerase systems,
or are fully synthetic constructs. There is no plausible route for these sequences to
appear in a bacterial assembly except deliberate cloning. A single high-confidence
ENGINEERED hit is sufficient to classify a contig as synthetic.

This tier includes:

- **Eukaryotic and viral promoters**: CMV, EF1α, PGK, CAG, CaMV 35S. These are active
  only in eukaryotic cells and have no function in bacteria. Their presence in a bacterial
  assembly is unambiguous.
- **Phage RNA polymerase promoters**: T7, SP6, T3. These require T7/SP6/T3 RNA
  polymerase, which is never endogenous in wild-type bacteria. Their use is exclusively
  in engineered expression systems.
- **Eukaryotic and viral replication origins**: SV40 ori, EBV oriP. Never functional in
  bacteria.
- **Filamentous phage origins used as cloning tools**: f1 ori, M13 ori. Phage-derived
  elements used in phagemid systems; they do not occur in bacterial chromosomes or natural
  conjugative plasmids.
- **Eukaryotic transcription terminators and polyadenylation signals**: BGH polyA, SV40
  polyA. Eukaryotic-specific; have no role in bacterial transcription.
- **Mammalian selectable markers**: Puromycin resistance (pac), Blasticidin resistance
  (bsr), Hygromycin resistance (hph) when present in a mammalian vector context.
- **Fluorescent and reporter proteins**: EGFP, mCherry, mTurquoise, mVenus, mRuby,
  tdTomato, Luciferase, and all GFP-family variants. Never encoded in wild-type bacterial
  genomes.
- **CRISPR and genome-editing elements**: Cas9, Cas12, Cpf1, sgRNA scaffolds. Engineered
  systems assembled into configurations that do not occur naturally.
- **Recombinase systems**: Cre recombinase, Flp recombinase, PhiC31 integrase. Derived
  from phage but used exclusively as engineered tools.
- **Site-specific recombination sites**: attB, attP, attL, attR (Lambda and PhiC31),
  loxP, loxN, FRT, rox, and designer variants, distinguished from naturally occurring
  phage att sites by their flanking sequence context stored in the catalog.
- **Synthetic regulatory elements**: IRES sequences, WPRE, P2A/T2A/E2A self-cleaving
  peptides, synthetic Kozak sequences, synthetic splice donor/acceptor sites.
  Eukaryotic-specific or fully synthetic.
- **Sequencing primer binding sites**: M13 forward/reverse, T7, SP6, T3, pUC, and similar
  universal sequencing primers. Designed for laboratory use; no meaningful natural
  occurrence.
- **Retroviral and lentiviral elements**: LTR sequences, lentiviral packaging signals
  (PSI), HIV-derived elements. Viral in origin; never present in bacterial genomes outside
  deliberate cloning.

### CONTEXT_DEPENDENT

Sequences in this tier **exist naturally in bacteria or natural bacterial plasmids**. Their
presence in an assembly is not automatically diagnostic — it depends on the host organism
being assembled and the broader genomic context. For example, a ColE1 origin found in an
*E. coli* assembly may reflect a natural endogenous plasmid; the same origin in a
*Salmonella* or *Bacillus* assembly is anomalous. Antibiotic resistance genes such as KanR
and AmpR are endemic in environmental and clinical bacterial populations via horizontal
gene transfer and cannot individually prove synthetic contamination.

Their evidential weight comes from **functional class co-occurrence**: a contig
simultaneously carrying a replication origin, a resistance marker, and a promoter from
this tier is almost certainly synthetic, whereas a contig with a single hit from one
functional class is ambiguous.

This tier includes:

- **Bacterial replication origins**: ColE1, pMB1, pUC ori, p15A, pSC101, RK2, pBBR1,
  pSa, and similar. Derived from natural *E. coli* and broad-host-range plasmids.
- **Transfer origins (oriT)**: IncP, IncQ, IncW, and similar conjugative transfer
  sequences. Found in natural conjugative plasmids.
- **Bacterial antibiotic resistance markers**: aphA (KanR), bla (AmpR), cat (CmR),
  tetA (TetR), aac(3) (GentR), aadA (SpecR), and similar. Originally derived from
  transposons; widely disseminated in natural populations via horizontal gene transfer.
- **T7 expression system elements**: T7 terminator, T7 tag. The T7 promoter itself is
  ENGINEERED; these associated elements are CONTEXT_DEPENDENT because the terminator
  occurs in some natural phage genomes and the tag is a protein epitope.
- **Inducible expression system components**: lacI, lacO arrays, tac promoter, trc
  promoter, IPTG-responsive elements. Derived from the natural *E. coli* lac operon but
  assembled into configurations not found in wild-type chromosomes.
- **Mobile element remnants**: Transposon end sequences (Tn3, Tn5, Tn10, Tn903) and
  IS-element borders incorporated into vector backbones during historical cloning from
  natural sources. These also occur in natural chromosomes and plasmids.
- **Protein-binding operator sequences**: lacO, tetO arrays, and similar, when cataloged
  in their natural-length form rather than as synthetic multimers.

### WEAK

Sequences in this tier are **widespread in wild-type bacterial genomes** and carry low
specificity for synthetic origin. A WEAK hit contributes to a verdict only when it
co-occurs with ENGINEERED or CONTEXT_DEPENDENT evidence. WEAK hits alone are never
sufficient for a `VECTOR` or `SUSPECTED` verdict.

This tier includes:

- **Common transcription terminators**: rrnB T1/T2, tonB terminator, lambda t0, fd
  terminator. These are standard terminators used in expression vectors precisely because
  they are efficient in many bacteria — and therefore also common in natural genomes.
- **Natural sigma-factor-dependent promoters**: lac, ara, trp, phoA promoters in their
  unmodified natural form. Widespread in Enterobacteriaceae and related organisms.
- **Ribosomal RNA and tRNA elements**: rRNA genes, tRNA genes. Present in all bacteria.
- **Miscellaneous features with natural counterparts**: catalog entries that share the
  property of having a recognisable wild-type bacterial ancestor but insufficient
  specificity for synthetic origin on their own.

## Catalog Categories

| Catalog file | Feature type | Evidence tier | Rationale |
| :--- | :--- | :--- | :--- |
| `catalog_rep_origin.fasta.gz` | Replication origins | CONTEXT_DEPENDENT | Derived from natural plasmid/phage replicons; natural counterparts exist in bacteria; tier reflects host-organism dependency |
| `catalog_oriT.fasta.gz` | Transfer origins | CONTEXT_DEPENDENT | Found in natural conjugative plasmids; not diagnostic alone |
| `catalog_misc_recomb.fasta.gz` | Recombination sites | ENGINEERED | Designed att/lox/FRT sites; no natural occurrence in wild-type bacteria |
| `catalog_RBS.fasta.gz` | Ribosome binding sites | CONTEXT_DEPENDENT / WEAK | Synthetic designer RBS: CONTEXT_DEPENDENT; minimal natural-sequence variants: WEAK |
| `catalog_CDS.fasta.gz` | Resistance markers, reporter genes, epitope tags | ENGINEERED / CONTEXT_DEPENDENT | Fluorescent proteins, Cas9, recombinases: ENGINEERED; Tn-derived resistance genes: CONTEXT_DEPENDENT |
| `catalog_promoter.fasta.gz` | Promoters | ENGINEERED / CONTEXT_DEPENDENT / WEAK | T7/CMV/CAG and phage RNAP promoters: ENGINEERED; tac/trc/lacI-regulated: CONTEXT_DEPENDENT; unmodified lac/ara/trp: WEAK |
| `catalog_terminator.fasta.gz` | Transcription terminators | ENGINEERED / WEAK | BGH/SV40 polyA: ENGINEERED; rrnB/tonB/lambda t0: WEAK |
| `catalog_primer_bind.fasta.gz` | Sequencing primer binding sites | ENGINEERED | Designed for laboratory use; no meaningful natural occurrence in bacteria |
| `catalog_protein_bind.fasta.gz` | Protein binding sites (operators) | CONTEXT_DEPENDENT | Derived from natural operator sequences; individually ambiguous |
| `catalog_regulatory.fasta.gz` | Synthetic regulatory elements | ENGINEERED | IRES, WPRE, P2A/T2A, synthetic Kozak, splice signals — eukaryotic-specific or fully synthetic |
| `catalog_enhancer.fasta.gz` | Enhancer elements | ENGINEERED | Viral and mammalian enhancers (CMV, SV40) have no natural occurrence in bacteria |
| `catalog_LTR.fasta.gz` | Retroviral LTR sequences | ENGINEERED | Appear in bacterial assemblies only through deliberate cloning |
| `catalog_mobile_element.fasta.gz` | Mobile element remnants | CONTEXT_DEPENDENT | Transposon and IS borders occur in natural chromosomes and plasmids |
| `catalog_misc_feature.fasta.gz` | Miscellaneous synthetic features | CONTEXT_DEPENDENT / WEAK | Mixed; assessed in combination with other evidence |

## Scoring Principle

One **ENGINEERED** hit with sufficient identity and coverage is treated as direct evidence
of synthetic origin and is sufficient for a `VECTOR` verdict on its own.

For **CONTEXT_DEPENDENT** hits, the scoring model weighs **functional class co-occurrence**:
a contig that carries hits from multiple independent functional classes (replication,
resistance, expression control) is much more likely to be synthetic than a contig with
hits from only one functional class. Hit count alone is not a reliable signal for
CONTEXT_DEPENDENT evidence.

**WEAK** hits are never scored in isolation. They contribute only when ENGINEERED or
CONTEXT_DEPENDENT evidence is already present on the same contig.

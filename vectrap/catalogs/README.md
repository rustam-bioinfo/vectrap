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

## Catalog Categories

| Catalog file | Feature type | Evidence tier |
| :--- | :--- | :--- |
| `catalog_rep_origin.fasta.gz` | Replication origins | PRIMARY |
| `catalog_oriT.fasta.gz` | Transfer origins | PRIMARY |
| `catalog_misc_recomb.fasta.gz` | Recombination sites | PRIMARY |
| `catalog_RBS.fasta.gz` | Ribosome binding sites | PRIMARY |
| `catalog_CDS.fasta.gz` | Coding sequences (resistance markers, tags) | SUPPORTIVE |
| `catalog_promoter.fasta.gz` | Engineered promoters | SUPPORTIVE |
| `catalog_terminator.fasta.gz` | Transcription terminators | SUPPORTIVE |
| `catalog_primer_bind.fasta.gz` | Sequencing primer binding sites | SUPPORTIVE |
| `catalog_protein_bind.fasta.gz` | Protein binding sites (operators) | SUPPORTIVE |
| `catalog_regulatory.fasta.gz` | Regulatory elements | SUPPORTIVE |
| `catalog_enhancer.fasta.gz` | Enhancer elements | SUPPORTIVE |
| `catalog_LTR.fasta.gz` | Retroviral LTR sequences | SUPPORTIVE |
| `catalog_mobile_element.fasta.gz` | Mobile element remnants | SUPPORTIVE |
| `catalog_misc_feature.fasta.gz` | Miscellaneous synthetic features | SUPPORTIVE |

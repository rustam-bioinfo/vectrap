# VecTrap Development Roadmap

This file tracks all pending tasks for the complete rewrite of VecTrap into a catalog-driven homology-based detection pipeline.

---

## 1. Catalog Deposition

- [ ] Upload all `catalog_*.fasta` files to Zenodo
- [ ] Obtain Zenodo DOI for the catalog dataset
- [ ] Record the mapping file `vectrap_catalog_mapping.tsv` locally (do not ship with tool, internal provenance only)
- [ ] Add Zenodo DOI badge and download link to `README.md`
- [ ] Decide on catalog versioning scheme (e.g. `v1.0`) for future updates

---

## 2. Repository Structure

- [ ] Create `vectrap/modules/` new source files (see section 3)
- [ ] Create `db/build_db.py` -- one-time index preparation script
- [ ] Create `vectrap.py` -- unified CLI entry point
- [ ] Remove all `.gitkeep` placeholders once real files are added
- [ ] Add `vectrap/catalogs/README.md` explaining catalog format and how to obtain catalog files from Zenodo
- [ ] Update `requirements.txt` with all real dependencies once scanning strategy is finalized
- [ ] Add `.gitignore` entries for BLAST db index files (`*.nin`, `*.nhr`, `*.nsq`, `*.nsi`, `*.nsd`, `*.not`, `*.ntf`, `*.nto`) and minimap2 indexes (`*.mmi`)

---

## 3. Core Module Development

### 3.1 Scanning strategy decision
- [ ] Decide between BLAST and minimap2 for long features (>= 50 bp)
- [ ] Decide between exact k-mer hashing and short-read aligner for short features (< 50 bp: `primer_bind`, `RBS`)

### 3.2 `vectrap/modules/utils.py`
- [ ] `read_fasta(path)` -- plain and gzipped FASTA reader
- [ ] `rev_comp(seq)` -- reverse complement
- [ ] `open_text(path)` -- transparent gz/plain opener

### 3.3 `vectrap/modules/homology_scanner.py`
- [ ] Implement long-feature homology scan (BLAST or minimap2)
- [ ] Implement short-feature exact k-mer hash scan
- [ ] Bidirectional strand scanning with 0-based coordinate output
- [ ] Per-catalog configurable identity and coverage thresholds
- [ ] Return unified `HomologyHit` dataclass per hit

### 3.4 `vectrap/modules/scorer.py`
- [ ] Define `CATALOG_TIERS` dict (PRIMARY vs SUPPORTIVE per feature type)
- [ ] Define `SUPPORTIVE_WEIGHTS` dict (per feature type)
- [ ] Implement per-contig evidence aggregation
- [ ] Implement classification logic: VECTOR / SUSPECTED / CLEAN
- [ ] Define and document `VECTOR_THRESHOLD` default value
- [ ] Validate weights empirically using known positive and negative control sequences

### 3.5 `vectrap.py`
- [ ] Argument parser: `-i`, `-o`, `-c`, `--min-identity`, `--min-coverage`
- [ ] Orchestrate scanner and scorer
- [ ] Write `*.all_hits.tsv` per-hit output
- [ ] Write `*.verdicts.tsv` per-contig classification output
- [ ] Print run summary to stdout

### 3.6 `db/build_db.py`
- [ ] `--download` mode: fetch catalog FASTA files from Zenodo DOI
- [ ] `--catalog-dir` mode: use locally provided FASTA files
- [ ] Build BLAST or minimap2 indexes depending on chosen strategy
- [ ] Validate that all expected catalog files are present after download

---

## 4. Testing

- [ ] Collect positive control sequences: known vector contaminants from NCBI with confirmed annotation
- [ ] Collect negative control sequences: real bacterial chromosomal contigs with no vector elements
- [ ] Write `tests/test_scanner.py` -- unit tests for hit detection on positive controls
- [ ] Write `tests/test_scorer.py` -- unit tests for classification logic
- [ ] Write `tests/test_cli.py` -- integration test for full pipeline run
- [ ] Define acceptable sensitivity and specificity thresholds
- [ ] Benchmark runtime on a representative assembly dataset

---

## 5. Documentation

- [ ] Finalize `README.md` with complete installation and usage instructions once CLI is stable
- [ ] Document all CLI arguments with examples
- [ ] Add `CHANGELOG.md`
- [ ] Add `CONTRIBUTING.md`
- [ ] Write methods section text describing the catalog construction procedure (for manuscript)

---

## 6. Publication Preparation

- [ ] Finalize manuscript draft
- [ ] Add VecTrap citation placeholder to `README.md` once preprint is posted
- [ ] Add Zenodo catalog DOI to manuscript methods
- [ ] Add tool version to all output file headers for reproducibility

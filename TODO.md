# VecTrap Development Roadmap

This file tracks all pending tasks for the complete rewrite of VecTrap into a catalog-driven homology-based detection pipeline.

---

## 1. Catalog Deposition

- [x] Upload all `catalog_*.fasta` files to Zenodo
- [x] Obtain Zenodo DOI for the catalog dataset -- `10.5281/zenodo.20844271`
- [x] Record the mapping file `vectrap_catalog_mapping.tsv` locally (do not ship with tool, internal provenance only)
- [x] Add Zenodo DOI badge and download link to `README.md`
- [x] Decide on catalog versioning scheme -- `v1.0`

---

## 2. Repository Structure

- [x] Restructure for PyPI packaging with `pyproject.toml`
- [x] Create `vectrap/cli/` with `build_db.py` and `run.py` entry points
- [x] Create `vectrap/modules/__init__.py` and `vectrap/__init__.py`
- [x] Add `vectrap/catalogs/README.md` explaining catalog format and Zenodo download
- [x] Remove old `db/` top-level directory
- [ ] Create `vectrap/modules/homology_scanner.py`
- [ ] Create `vectrap/modules/scorer.py`
- [ ] Wire scanner and scorer into `vectrap/cli/run.py`
- [ ] Update `requirements.txt` with all dependencies once scanning strategy is finalized
- [ ] Add `.gitignore` entries for minimap2 indexes (`*.mmi`) and pickle files (`*.pkl`)

---

## 3. Core Module Development

### 3.1 Scanning strategy
- [x] Use minimap2 for long features (>= 50 bp)
- [x] Use exact k-mer hashing for short features (< 50 bp)

### 3.2 `vectrap/modules/utils.py`
- [x] `read_fasta(path)` -- plain and gzipped FASTA reader
- [x] `rev_comp(seq)` -- reverse complement
- [x] `open_text(path)` -- transparent gz/plain opener
- [x] `write_fasta(records, path)` -- plain and gzipped FASTA writer

### 3.3 `vectrap/cli/build_db.py`
- [x] `--download` mode: fetch catalog FASTA files from Zenodo DOI `10.5281/zenodo.20844271`
- [x] `--catalog-dir` mode: use locally provided FASTA files
- [x] Length-based routing: long to minimap2 index, short to k-mer hash pickle
- [x] Validate downloaded files against `catalog_manifest.tsv` checksums

### 3.4 `vectrap/modules/homology_scanner.py`
- [ ] `HomologyHit` dataclass (contig, start, end, strand, identity, coverage, catalog_id)
- [ ] minimap2 PAF parser
- [ ] minimap2 subprocess runner with identity and coverage filters
- [ ] Exact k-mer hash scanner for short sequences
- [ ] Bidirectional strand handling with 0-based coordinate output
- [ ] Single `scan(query_fasta, catalog_dir, min_identity, min_coverage)` interface

### 3.5 `vectrap/modules/scorer.py`
- [ ] Define `CATALOG_TIERS` dict (PRIMARY vs SUPPORTIVE per feature type)
- [ ] Define `SUPPORTIVE_WEIGHTS` dict (per feature type)
- [ ] Implement per-contig evidence aggregation
- [ ] Implement classification logic: VECTOR / SUSPECTED / CLEAN
- [ ] Define and document `VECTOR_THRESHOLD` default value
- [ ] Validate weights empirically using known positive and negative control sequences

### 3.6 `vectrap/cli/run.py`
- [x] Argument parser: `-i`, `-o`, `-c`, `--min-identity`, `--min-coverage`, `--version`
- [ ] Wire scanner and scorer
- [ ] Write `*.all_hits.tsv` per-hit output
- [ ] Write `*.verdicts.tsv` per-contig classification output
- [ ] Print run summary to stdout

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
- [ ] Publish to PyPI

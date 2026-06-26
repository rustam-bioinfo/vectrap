# VecTrap Development Roadmap

---

## 1. Catalog & Metadata

- [x] Upload all `catalog_*.fasta` files to Zenodo
- [x] Obtain Zenodo DOI — `10.5281/zenodo.20844271`
- [x] Add `tier`, `confidence`, `reasoning`, `label` columns to `catalog_manifest.tsv`
- [x] Record internal provenance mapping (`vectrap_catalog_mapping.tsv`) — not shipped with tool
- [x] Add Zenodo DOI badge and download link to `README.md`
- [x] Decide on catalog versioning scheme — `v1.0`

---

## 2. Repository & Packaging

- [x] Restructure for PyPI packaging with `pyproject.toml`
- [x] Create `vectrap/cli/` with `build_db.py` and `run.py` entry points
- [x] Create `vectrap/modules/__init__.py` and `vectrap/__init__.py`
- [x] Add `vectrap/catalogs/README.md` explaining catalog format and Zenodo download
- [x] Set all dependencies in `pyproject.toml` (`mappy`, `pyahocorasick`, `numpy`, `biopython`)
- [x] Delete `requirements.txt` — redundant with `pyproject.toml`
- [ ] Add `.gitignore` entries for `*.mmi` and `*.pkl` (currently missing from `.gitignore`)

---

## 3. Core Modules

### 3.1 `vectrap/modules/utils.py` — COMPLETE
- [x] `read_fasta`, `rev_comp`, `open_text`, `write_fasta`

### 3.2 `vectrap/cli/build_db.py` — COMPLETE
- [x] `--download` and `--catalog-dir` modes
- [x] Length-based routing: long → `combined_long.mmi`, short → `short_index.pkl`
- [x] MD5 checksum validation against `catalog_manifest.tsv`
- [x] Build `catalog_metadata.pkl` from `catalog_manifest.tsv`

### 3.3 `vectrap/modules/homology_scanner.py` — COMPLETE
- [x] `HomologyHit` dataclass with all fields (`tier`, `confidence`, `reasoning`, `label`, `feature_type`, `source`)
- [x] mappy aligner with identity and coverage filters
- [x] Aho-Corasick k-mer scanner for short sequences
- [x] Forward-strand 0-based coordinate output for both scanners
- [x] Metadata enrichment from `catalog_metadata.pkl`
- [x] Indexes accepted as pre-loaded objects (not reloaded per sample)

### 3.4 `vectrap/modules/scorer.py` — COMPLETE
- [x] `ContigSummary` dataclass
- [x] Interval merging for `covered_bp` / `covered_fraction`
- [x] `evidence_summary` string generation
- [x] `summarize(hits, contig_lengths)` public interface

### 3.5 `vectrap/modules/reporter.py` — COMPLETE
- [x] `build_sample_row()` with `CONTAMINATION` / `SUSPECTED` / `CLEAN` verdict logic
- [x] `write_summary()` writing `vectrap_summary.tsv`
- [x] `ERROR` row on per-sample exception

### 3.6 `vectrap/cli/run.py` — COMPLETE
- [x] All CLI arguments (`-i`, `-o`, `-c`, `--min-identity`, `--min-coverage`, `--min-len`, `--best-n`, `--threads`, `--min-engineered-contamination`, `--min-context-suspected`, `-v`, `--version`)
- [x] Single-file and directory batch modes
- [x] Indexes loaded once before the sample loop
- [x] `{stem}_hits.tsv` and `{stem}_contig_verdicts.tsv` per sample
- [x] `vectrap_summary.tsv` across all samples
- [x] Per-sample error isolation with partial file cleanup

---

## 4. Testing

- [x] `tests/test_utils.py` — `read_fasta`, `rev_comp`, `open_text`, `write_fasta`
- [x] `tests/test_build_db.py` — index building, manifest validation, routing logic
- [x] `tests/test_scanner.py` — hit detection, coordinate correctness, metadata enrichment
- [ ] `tests/test_scorer.py` — unit tests for `summarize()`: covered_bp merging, evidence_summary format, sort order, edge cases (no hits, unannotated-only)
- [ ] `tests/test_reporter.py` — verdict thresholds (`CONTAMINATION` / `SUSPECTED` / `CLEAN`), `ERROR` row format, `vectrap_summary.tsv` column completeness
- [ ] `tests/test_cli.py` — integration test: build fake catalog → run pipeline on synthetic FASTA → assert output TSVs exist and contain expected verdict
- [ ] Collect positive control sequences: known vector-contaminated contigs with confirmed annotation from NCBI
- [ ] Collect negative control sequences: real bacterial chromosomal contigs with no vector elements
- [ ] Define and document acceptable sensitivity / specificity thresholds
- [ ] Benchmark runtime on a representative dataset (e.g. 100 assemblies, varied contig counts)

---

## 5. Documentation

- [x] `README.md` — installation, quickstart, CLI reference, output file schemas
- [x] `CONTEXT.md` — up-to-date coding-session onboarding document
- [ ] `CHANGELOG.md` — create; add entry for v0.1.0
- [ ] `CONTRIBUTING.md` — create; describe how to add new catalog entries, run tests, submit PRs
- [ ] Write methods section text for manuscript describing catalog construction and scoring logic

---

## 6. Publication Preparation

- [ ] Finalize manuscript draft
- [ ] Add tool version to all output file headers for reproducibility (currently absent from TSV headers)
- [ ] Add VecTrap citation placeholder to `README.md` once preprint is posted
- [ ] Add Zenodo catalog DOI to manuscript methods
- [ ] Publish to PyPI (package is installable locally; PyPI upload not yet done)

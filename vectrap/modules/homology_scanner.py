"""Homology-based vector sequence scanner.

This module provides the core scanning engine for VecTrap. It exposes a single
public function, ``scan()``, which runs two complementary detection strategies
against the query assembly:

1. **mappy aligner** -- for long catalog sequences (>= min_len bp).
2. **Aho-Corasick k-mer scanner** -- for short catalog sequences (< min_len bp).

Performance notes
-----------------
- The mappy index (``.mmi``) is loaded once via ``load_mappy_index()`` and
  reused across all samples in a batch run.
- The Aho-Corasick automaton is built once via ``build_kmer_automaton()`` and
  reused across all samples.
- Each FASTA file is read in a single pass that feeds both the contig-length
  dict and the scanners, avoiding redundant I/O.
"""

from __future__ import annotations

import pickle
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import ahocorasick
import mappy

from vectrap.modules.utils import rev_comp

_DEFAULT_MIN_LEN: int = 50
_DEFAULT_BEST_N: int = 10


def _ts() -> str:
    return datetime.now().strftime("[%H:%M:%S]")


def _elapsed(t0: float) -> str:
    secs = time.time() - t0
    if secs < 60:
        return f"{secs:.1f}s"
    m, s = divmod(int(secs), 60)
    return f"{m}m{s:02d}s"


# ---------------------------------------------------------------------------
# Public data structure
# ---------------------------------------------------------------------------

@dataclass
class HomologyHit:
    """A single homology match between a catalog sequence and a query contig."""

    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float
    catalog_id: str
    source: str = field(default="mappy")
    feature_type: str = field(default="")
    label: str = field(default="")
    tier: str = field(default="")
    confidence: str = field(default="")
    reasoning: str = field(default="")

    @property
    def length(self) -> int:
        return self.end - self.start

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"HomologyHit({self.contig!r}, {self.start}-{self.end} "
            f"{self.strand}, id={self.identity:.3f}, cov={self.coverage:.3f}, "
            f"cat={self.catalog_id!r}, tier={self.tier!r}, src={self.source!r})"
        )


# ---------------------------------------------------------------------------
# Catalog metadata loader
# ---------------------------------------------------------------------------

def _load_metadata(catalog_dir: Path) -> Dict[str, dict]:
    pkl_path = catalog_dir / "catalog_metadata.pkl"
    if not pkl_path.exists():
        print(
            f"[vectrap] WARNING: catalog_metadata.pkl not found in {catalog_dir}.\n"
            "  Run 'vectrap-build-db' to build it. Hits will lack tier/confidence/reasoning.",
            file=sys.stderr,
        )
        return {}
    with open(pkl_path, "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------------------
# One-time index loaders (call once per batch, reuse across samples)
# ---------------------------------------------------------------------------

def load_mappy_index(
    index_path: Path,
    best_n: int = _DEFAULT_BEST_N,
    threads: int = 4,
) -> mappy.Aligner:
    """Load the mappy minimap2 index once and return the ``Aligner`` object.

    The returned aligner should be reused for all samples in a batch run.

    Parameters
    ----------
    index_path : Path
        Path to the ``.mmi`` index file built by ``vectrap-build-db``.
    best_n : int
        Maximum alignments reported per query sequence.
    threads : int
        Number of aligner threads.
    """
    aligner = mappy.Aligner(str(index_path), preset="map-ont", n_threads=threads, best_n=best_n)
    if not aligner:
        raise RuntimeError(
            f"mappy could not load index: {index_path}\n"
            "Run 'vectrap-build-db' to rebuild the catalog indexes."
        )
    return aligner


def build_kmer_automaton(pkl_path: Path) -> Optional[ahocorasick.Automaton]:
    """Build the Aho-Corasick automaton from the short-sequence index once.

    Returns ``None`` if the index is empty.
    The returned automaton should be reused for all samples in a batch run.

    Parameters
    ----------
    pkl_path : Path
        Path to ``short_index.pkl`` built by ``vectrap-build-db``.
    """
    with open(pkl_path, "rb") as fh:
        short_index: dict = pickle.load(fh)

    if not short_index:
        return None

    A = ahocorasick.Automaton()
    for kmer, catalog_ids in short_index.items():
        klen = len(kmer)
        fwd_payload = [(cid, "+", klen) for cid in catalog_ids]
        if kmer in A:
            A.get(kmer).extend(fwd_payload)
        else:
            A.add_word(kmer, list(fwd_payload))

        rc = rev_comp(kmer)
        if rc == kmer:
            continue
        rc_payload = [(cid, "-", klen) for cid in catalog_ids]
        if rc in A:
            A.get(rc).extend(rc_payload)
        else:
            A.add_word(rc, list(rc_payload))

    A.make_automaton()
    return A


# ---------------------------------------------------------------------------
# Per-sample scanners (accept pre-built index objects)
# ---------------------------------------------------------------------------

def _identity_from_hit(hit: "mappy.Alignment") -> float:
    if hit.blen > 0:
        return hit.mlen / hit.blen
    return 0.0


def _scan_with_mappy(
    contigs: list[tuple[str, str]],
    aligner: mappy.Aligner,
    min_identity: float,
    min_coverage: float,
) -> List[HomologyHit]:
    """Align *contigs* against a pre-loaded mappy aligner."""
    hits: List[HomologyHit] = []
    for contig_name, seq in contigs:
        for hit in aligner.map(seq, cs=True):
            identity = _identity_from_hit(hit)
            if identity < min_identity:
                continue
            coverage = (hit.r_en - hit.r_st) / hit.ctg_len if hit.ctg_len > 0 else 0.0
            if coverage < min_coverage:
                continue
            hits.append(HomologyHit(
                contig=contig_name,
                start=hit.q_st,
                end=hit.q_en,
                strand="+" if hit.strand == 1 else "-",
                identity=identity,
                coverage=coverage,
                catalog_id=hit.ctg,
                source="mappy",
            ))
    return hits


def _scan_with_automaton(
    contigs: list[tuple[str, str]],
    automaton: ahocorasick.Automaton,
) -> List[HomologyHit]:
    """Scan *contigs* using a pre-built Aho-Corasick automaton."""
    hits: List[HomologyHit] = []
    for contig_name, seq in contigs:
        for end_idx, payload in automaton.iter(seq):
            for catalog_id, strand, klen in payload:
                hits.append(HomologyHit(
                    contig=contig_name,
                    start=end_idx - klen + 1,
                    end=end_idx + 1,
                    strand=strand,
                    identity=1.0,
                    coverage=1.0,
                    catalog_id=catalog_id,
                    source="kmer",
                ))
    return hits


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def scan(
    query_fasta: str | Path,
    catalog_dir: str | Path,
    min_identity: float = 0.90,
    min_coverage: float = 0.80,
    min_len: int = _DEFAULT_MIN_LEN,
    best_n: int = _DEFAULT_BEST_N,
    threads: int = 4,
    verbose: bool = False,
    # Optional pre-built objects — pass these in batch runs to avoid
    # re-loading the index and rebuilding the automaton for every sample.
    aligner: Optional[mappy.Aligner] = None,
    automaton: Optional[ahocorasick.Automaton] = None,
    metadata: Optional[Dict[str, dict]] = None,
) -> tuple[List[HomologyHit], Dict[str, int]]:
    """Scan *query_fasta* against the VecTrap catalog.

    Parameters
    ----------
    query_fasta : str or Path
        Input assembly FASTA file (plain or gzipped).
    catalog_dir : str or Path
        Directory containing the indexes built by ``vectrap-build-db``.
    min_identity : float
        Minimum per-base sequence identity for mappy hits (default: 0.90).
    min_coverage : float
        Minimum fraction of the catalog sequence covered by a hit (default: 0.80).
    min_len : int
        Sequences >= this length use mappy; shorter use exact k-mer matching.
        Must match the value used during index build (default: 50).
    best_n : int
        Maximum number of mappy alignments reported per query sequence (default: 10).
    threads : int
        Number of threads passed to mappy (default: 4).
    verbose : bool
        Print timestamped progress messages to stderr (default: False).
    aligner : mappy.Aligner, optional
        Pre-loaded mappy aligner. If omitted, loaded from *catalog_dir*.
        Pass a shared instance in batch runs to avoid redundant index loading.
    automaton : ahocorasick.Automaton, optional
        Pre-built Aho-Corasick automaton. If omitted, built from *catalog_dir*.
        Pass a shared instance in batch runs to avoid redundant automaton builds.
    metadata : dict, optional
        Pre-loaded catalog metadata dict. If omitted, loaded from *catalog_dir*.

    Returns
    -------
    tuple[list[HomologyHit], dict[str, int]]
        ``(hits, contig_lengths)`` — all enriched hits plus contig name->length
        mapping, both derived from a single FASTA pass.
    """
    query_fasta = Path(query_fasta)
    catalog_dir = Path(catalog_dir)

    mmi_path = catalog_dir / "combined_long.mmi"
    pkl_path = catalog_dir / "short_index.pkl"

    for p in (mmi_path, pkl_path):
        if not p.exists():
            raise FileNotFoundError(
                f"Catalog index not found: {p}\n"
                "Run 'vectrap-build-db' to build the catalog indexes."
            )

    def _log(msg: str) -> None:
        if verbose:
            print(f"{_ts()} {msg}", file=sys.stderr, flush=True)

    t_total = time.time()

    _log(f"input            : {query_fasta}")
    _log(f"catalog          : {catalog_dir}")
    _log(f"min identity     : {min_identity}")
    _log(f"min coverage     : {min_coverage}")
    _log(f"min len          : {min_len}")
    _log(f"best n           : {best_n}")
    _log(f"threads          : {threads}")

    # --- single-pass FASTA read ---
    _log("\u2500" * 48)
    _log("stage 1/3  reading FASTA")
    t0 = time.time()
    import gzip
    contigs: list[tuple[str, str]] = []
    contig_lengths: Dict[str, int] = {}

    def _open(p: Path):
        if str(p).endswith(".gz"):
            return gzip.open(p, "rt")
        return open(p, "r")

    with _open(query_fasta) as fh:
        name = None
        buf: list[str] = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(buf)
                    contigs.append((name, seq))
                    contig_lengths[name] = len(seq)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seq = "".join(buf)
            contigs.append((name, seq))
            contig_lengths[name] = len(seq)

    _log(f"    contigs read   : {len(contigs):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    # --- mappy aligner ---
    _log("\u2500" * 48)
    _log("stage 2/3  mappy aligner (long sequences)")
    t0 = time.time()
    if aligner is None:
        _log("    loading mappy index ...")
        aligner = load_mappy_index(mmi_path, best_n=best_n, threads=threads)
    mappy_hits = _scan_with_mappy(contigs, aligner, min_identity, min_coverage)
    _log(f"    hits           : {len(mappy_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    # --- k-mer scanner ---
    _log("\u2500" * 48)
    _log("stage 3/3  k-mer scanner (short sequences)")
    t0 = time.time()
    kmer_hits: List[HomologyHit] = []
    if automaton is None:
        _log("    building automaton ...")
        automaton = build_kmer_automaton(pkl_path)
    if automaton is not None:
        kmer_hits = _scan_with_automaton(contigs, automaton)
    _log(f"    hits           : {len(kmer_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    all_hits = mappy_hits + kmer_hits

    # --- enrich with metadata ---
    if metadata is None:
        metadata = _load_metadata(catalog_dir)
    unannotated = 0
    for hit in all_hits:
        m = metadata.get(hit.catalog_id)
        if m is None:
            unannotated += 1
            continue
        hit.feature_type = m.get("feature_type", "")
        hit.label        = m.get("label", "")
        hit.tier         = m.get("tier", "")
        hit.confidence   = m.get("confidence", "")
        hit.reasoning    = m.get("reasoning", "")
    if unannotated and verbose:
        _log(f"    WARNING: {unannotated:,} hit(s) had no metadata entry")

    _log("\u2500" * 48)
    _log(f"total hits        : {len(all_hits):,}")
    _log(f"total elapsed     : {_elapsed(t_total)}")

    return all_hits, contig_lengths

"""Homology-based vector sequence scanner.

This module provides the core scanning engine for VecTrap. It exposes a single
public function, ``scan()``, which runs two complementary detection strategies
against the query assembly:

1. **mappy aligner** -- for long catalog sequences (>= min_len bp).
2. **Aho-Corasick k-mer scanner** -- for short catalog sequences (< min_len bp).

All hits from both strategies are returned as a flat list of ``HomologyHit``
objects, enriched with tier/confidence/reasoning from the catalog metadata.
"""

from __future__ import annotations

import pickle
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import ahocorasick
import mappy

from vectrap.modules.utils import read_fasta, rev_comp

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
# mappy aligner (long sequences)
# ---------------------------------------------------------------------------

def _identity_from_hit(hit: "mappy.Alignment") -> float:
    if hit.blen > 0:
        return hit.mlen / hit.blen
    return 0.0


def _run_mappy(
    query_fasta: Path,
    index_path: Path,
    min_identity: float,
    min_coverage: float,
    best_n: int = _DEFAULT_BEST_N,
    threads: int = 4,
    log=None,
) -> List[HomologyHit]:
    if log is None:
        log = lambda _: None

    aligner = mappy.Aligner(str(index_path), preset="map-ont", n_threads=threads, best_n=best_n)
    if not aligner:
        raise RuntimeError(
            f"mappy could not load index: {index_path}\n"
            "Run 'vectrap-build-db' to rebuild the catalog indexes."
        )

    hits: List[HomologyHit] = []
    n_contigs = 0
    for header, seq in read_fasta(query_fasta):
        contig_name = header.split()[0]
        n_contigs += 1
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

    log(f"    contigs scanned : {n_contigs:,}")
    return hits


# ---------------------------------------------------------------------------
# Aho-Corasick k-mer scanner (short sequences)
# ---------------------------------------------------------------------------

def _build_automaton_with_klen(short_index: dict) -> ahocorasick.Automaton:
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


def _kmer_scan_sequence(
    contig_name: str,
    contig_seq: str,
    automaton: ahocorasick.Automaton,
) -> List[HomologyHit]:
    hits: List[HomologyHit] = []
    for end_idx, payload in automaton.iter(contig_seq):
        for catalog_id, strand, klen in payload:
            start = end_idx - klen + 1
            end   = end_idx + 1
            hits.append(HomologyHit(
                contig=contig_name,
                start=start,
                end=end,
                strand=strand,
                identity=1.0,
                coverage=1.0,
                catalog_id=catalog_id,
                source="kmer",
            ))
    return hits


def _run_kmer_scanner(
    query_fasta: Path,
    pkl_path: Path,
    log=None,
) -> List[HomologyHit]:
    if log is None:
        log = lambda _: None

    with open(pkl_path, "rb") as fh:
        short_index: dict = pickle.load(fh)

    if not short_index:
        return []

    log(f"    building automaton ({len(short_index):,} patterns) ...")
    t_ac = time.time()
    automaton = _build_automaton_with_klen(short_index)
    log(f"    automaton ready  : {time.time() - t_ac:.2f}s")

    hits: List[HomologyHit] = []
    n_contigs = 0
    for header, seq in read_fasta(query_fasta):
        contig_name = header.split()[0]
        n_contigs += 1
        hits.extend(_kmer_scan_sequence(contig_name, seq, automaton))

    log(f"    contigs scanned : {n_contigs:,}")
    log(f"    k-mers in index : {len(short_index):,}")
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
) -> List[HomologyHit]:
    """Scan *query_fasta* against the VecTrap catalog and return all hits.

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
        Sequences >= this length use the mappy aligner; shorter sequences use
        exact k-mer matching. Must match the value used during index build (default: 50).
    best_n : int
        Maximum number of mappy alignments reported per query sequence (default: 10).
    threads : int
        Number of threads passed to mappy (default: 4).
    verbose : bool
        Print timestamped progress messages to stderr (default: False).
    """
    query_fasta = Path(query_fasta)
    catalog_dir = Path(catalog_dir)

    mmi_path = catalog_dir / "combined_long.mmi"
    pkl_path = catalog_dir / "short_index.pkl"

    if not mmi_path.exists():
        raise FileNotFoundError(
            f"minimap2 index not found: {mmi_path}\n"
            "Run 'vectrap-build-db' to build the catalog indexes."
        )
    if not pkl_path.exists():
        raise FileNotFoundError(
            f"k-mer index not found: {pkl_path}\n"
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

    # --- mappy aligner ---
    _log("\u2500" * 48)
    _log("stage 1/3  mappy aligner (long sequences)")
    t0 = time.time()
    mappy_hits = _run_mappy(
        query_fasta=query_fasta,
        index_path=mmi_path,
        min_identity=min_identity,
        min_coverage=min_coverage,
        best_n=best_n,
        threads=threads,
        log=_log,
    )
    _log(f"    hits           : {len(mappy_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    # --- k-mer scanner ---
    _log("\u2500" * 48)
    _log("stage 2/3  k-mer scanner (short sequences)")
    t0 = time.time()
    kmer_hits = _run_kmer_scanner(
        query_fasta=query_fasta,
        pkl_path=pkl_path,
        log=_log,
    )
    _log(f"    hits           : {len(kmer_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    all_hits = mappy_hits + kmer_hits

    # --- enrich hits with catalog metadata ---
    _log("\u2500" * 48)
    _log("stage 3/3  enriching hits with catalog metadata")
    t0 = time.time()
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
    if unannotated:
        _log(f"    WARNING: {unannotated:,} hit(s) had no metadata entry")
    _log(f"    elapsed        : {_elapsed(t0)}")

    _log("\u2500" * 48)
    _log(f"total hits        : {len(all_hits):,}")
    _log(f"total elapsed     : {_elapsed(t_total)}")

    return all_hits

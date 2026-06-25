"""Homology-based vector sequence scanner.

This module provides the core scanning engine for VecTrap. It exposes a single
public function, ``scan()``, which runs two complementary detection strategies
against the query assembly:

1. **mappy aligner** -- for long catalog sequences (>= MIN_LEN bp).
   Uses ``mappy.Aligner`` (the official Python bindings for minimap2) against
   the pre-built ``combined_long.mmi`` index from ``vectrap-build-db``.  No
   external binary is required -- mappy is a pure pip dependency.

2. **Aho-Corasick k-mer scanner** -- for short catalog sequences (< MIN_LEN bp).
   Builds an ``ahocorasick.Automaton`` from the short-sequence index, then
   scans every contig in a single O(contig_len) pass per strand.  This
   replaces the previous O(kmers * contig_len) str.find loop and is
   typically 100-1000x faster for the catalog sizes used in VecTrap.

All hits from both strategies are returned as a flat list of ``HomologyHit``
objects for downstream scoring.
"""

from __future__ import annotations

import pickle
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List

import ahocorasick
import mappy

from vectrap.modules.utils import read_fasta, rev_comp

# Sequences shorter than this threshold are stored in the k-mer hash index;
# sequences >= this threshold are indexed by mappy.  Must match the value
# used in vectrap/cli/build_db.py when the indexes were built.
MIN_LEN: int = 50


def _ts() -> str:
    """Return a formatted timestamp string for log messages: [HH:MM:SS]"""
    return datetime.now().strftime("[%H:%M:%S]")


def _elapsed(t0: float) -> str:
    """Return a human-readable elapsed time string since *t0* (from time.time())."""
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
    """A single homology match between a catalog sequence and a query contig.

    All coordinates are 0-based and refer to the forward strand of the query
    contig, regardless of the strand on which the match was found.

    Attributes
    ----------
    contig : str
        FASTA header of the query contig (first whitespace-delimited token).
    start : int
        0-based start coordinate on the forward strand.
    end : int
        0-based end coordinate on the forward strand (exclusive).
    strand : str
        ``'+'`` if the match is on the forward strand, ``'-'`` if on the
        reverse complement strand.
    identity : float
        Sequence identity as a fraction in [0, 1].  ``1.0`` for exact k-mer
        hits; derived from mappy alignment for long-sequence hits.
    coverage : float
        Fraction of the *catalog* sequence covered by the hit, in [0, 1].
        ``1.0`` for exact k-mer hits.
    catalog_id : str
        Identifier of the matching catalog entry.
    source : str
        ``'mappy'`` or ``'kmer'`` -- indicates which scanner produced the hit.
    """

    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float
    catalog_id: str
    source: str = field(default="mappy")

    @property
    def length(self) -> int:
        """Length of the hit on the query contig in bp."""
        return self.end - self.start

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"HomologyHit({self.contig!r}, {self.start}-{self.end} "
            f"{self.strand}, id={self.identity:.3f}, cov={self.coverage:.3f}, "
            f"cat={self.catalog_id!r}, src={self.source!r})"
        )


# ---------------------------------------------------------------------------
# mappy aligner (long sequences)
# ---------------------------------------------------------------------------

def _identity_from_hit(hit: "mappy.Alignment") -> float:
    """Derive per-base sequence identity from a mappy Alignment object.

    mappy exposes ``mlen`` (number of matching bases) and ``blen`` (alignment
    block length including gaps).  Identity = mlen / blen when blen > 0.
    """
    if hit.blen > 0:
        return hit.mlen / hit.blen
    return 0.0


def _run_mappy(
    query_fasta: Path,
    index_path: Path,
    min_identity: float,
    min_coverage: float,
    threads: int = 4,
    log=None,
) -> List[HomologyHit]:
    """Align *query_fasta* against *index_path* using mappy and return hits."""
    if log is None:
        log = lambda _: None

    aligner = mappy.Aligner(str(index_path), preset="map-ont", n_threads=threads, best_n=10)
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

def _build_automaton(short_index: dict) -> ahocorasick.Automaton:
    """Build an Aho-Corasick automaton from *short_index*.

    Each pattern is stored with its list of catalog IDs as the payload.
    Both the forward k-mer and (where non-palindromic) its reverse complement
    are added so a single forward-strand scan catches both orientations.

    The payload for each pattern is a list of (catalog_id, strand) tuples
    where strand is '+' for the original k-mer and '-' for the RC.
    """
    A = ahocorasick.Automaton()
    for kmer, catalog_ids in short_index.items():
        # Forward pattern
        fwd_payload = [(cid, "+") for cid in catalog_ids]
        if kmer in A:
            A.get(kmer).extend(fwd_payload)
        else:
            A.add_word(kmer, list(fwd_payload))

        # Reverse-complement pattern (skip palindromes)
        rc = rev_comp(kmer)
        if rc == kmer:
            continue
        rc_payload = [(cid, "-") for cid in catalog_ids]
        if rc in A:
            A.get(rc).extend(rc_payload)
        else:
            A.add_word(rc, list(rc_payload))

    A.make_automaton()
    return A


def _ac_scan_sequence(
    contig_name: str,
    contig_seq: str,
    automaton: ahocorasick.Automaton,
) -> List[HomologyHit]:
    """Scan *contig_seq* with the Aho-Corasick automaton in a single O(n) pass.

    The automaton contains both forward k-mers and their reverse complements,
    so one pass over the forward strand is sufficient.

    For a match ending at position ``end_idx`` (0-based, inclusive as returned
    by ahocorasick):
      - Forward hit  (strand='+'):  start = end_idx - klen + 1,  end = end_idx + 1
      - RC hit       (strand='-'):  the RC pattern matched at [start, end) on the
        forward strand, meaning the *original* k-mer is on the minus strand at the
        same coordinates.  We report the forward-strand coordinates of the RC match
        directly (start=end_idx-klen+1, end=end_idx+1) with strand='-'.
    """
    hits: List[HomologyHit] = []
    seq_len = len(contig_seq)

    for end_idx, payload in automaton.iter(contig_seq):
        # payload is a list of (catalog_id, strand) tuples
        # ahocorasick returns end_idx as the index of the LAST character of the match
        for catalog_id, strand in payload:
            # infer klen from the matched pattern length
            # ahocorasick doesn't expose pattern length directly but we can get it
            # from the payload key length -- instead, use end_idx and a sentinel
            # approach: we stored klen in the payload during build.
            # WORKAROUND: we don't have klen here; rebuild with klen in payload.
            pass

    # The above approach requires klen in the payload.  Rebuild is done in
    # _build_automaton_with_klen below -- this function is superseded.
    return hits


def _build_automaton_with_klen(short_index: dict) -> ahocorasick.Automaton:
    """Build an Aho-Corasick automaton where each payload includes the k-mer length.

    Payload per pattern: list of (catalog_id, strand, klen) tuples.
    """
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
    """Scan *contig_seq* with the Aho-Corasick automaton in a single O(n) pass.

    Parameters
    ----------
    contig_name : str
        FASTA header of the contig (first token).
    contig_seq : str
        Nucleotide sequence (upper-case).
    automaton : ahocorasick.Automaton
        Built by ``_build_automaton_with_klen``.  Payload per entry is a list
        of ``(catalog_id, strand, klen)`` tuples.

    Returns
    -------
    list[HomologyHit]
        All exact matches with ``identity=1.0``, ``coverage=1.0``.
    """
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
    """Load the short-sequence k-mer index, build an Aho-Corasick automaton,
    and scan all contigs in *query_fasta* in a single pass per contig.

    Parameters
    ----------
    query_fasta : Path
        Input assembly FASTA (plain or gzipped).
    pkl_path : Path
        Path to the ``short_index.pkl`` file built by ``vectrap-build-db``.
    log : callable or None
        Optional logging function.

    Returns
    -------
    list[HomologyHit]
        All exact k-mer matches.
    """
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
    threads: int = 4,
    verbose: bool = False,
) -> List[HomologyHit]:
    """Scan *query_fasta* against the VecTrap catalog and return all hits.

    This function runs both the mappy aligner (long sequences) and the exact
    k-mer hash scanner (short sequences) and returns their combined output as
    a flat list of ``HomologyHit`` objects.

    Parameters
    ----------
    query_fasta : str or Path
        Input assembly FASTA file (plain or gzipped).
    catalog_dir : str or Path
        Directory containing the indexes built by ``vectrap-build-db``:
        ``combined_long.mmi`` and ``short_index.pkl``.
    min_identity : float, optional
        Minimum per-base sequence identity for mappy hits (default: 0.90).
        Exact k-mer hits always have identity 1.0 regardless of this setting.
    min_coverage : float, optional
        Minimum fraction of the catalog sequence that must be covered by a hit
        (default: 0.80). Exact k-mer hits always have coverage 1.0.
    threads : int, optional
        Number of threads passed to mappy (default: 4).
    verbose : bool, optional
        Print timestamped progress messages to stderr (default: False).

    Returns
    -------
    list[HomologyHit]
        All hits passing the identity and coverage thresholds, from both
        scanners combined.

    Raises
    ------
    FileNotFoundError
        If the catalog index files are missing.
    RuntimeError
        If the mappy index cannot be loaded.
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
    _log(f"threads          : {threads}")

    # --- mappy aligner ---
    _log("─" * 48)
    _log("stage 1/2  mappy aligner (long sequences)")
    t0 = time.time()
    mappy_hits = _run_mappy(
        query_fasta=query_fasta,
        index_path=mmi_path,
        min_identity=min_identity,
        min_coverage=min_coverage,
        threads=threads,
        log=_log,
    )
    _log(f"    hits           : {len(mappy_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    # --- k-mer scanner ---
    _log("─" * 48)
    _log("stage 2/2  k-mer scanner (short sequences)")
    t0 = time.time()
    kmer_hits = _run_kmer_scanner(
        query_fasta=query_fasta,
        pkl_path=pkl_path,
        log=_log,
    )
    _log(f"    hits           : {len(kmer_hits):,}")
    _log(f"    elapsed        : {_elapsed(t0)}")

    # --- summary ---
    all_hits = mappy_hits + kmer_hits
    _log("─" * 48)
    _log(f"total hits        : {len(all_hits):,}")
    _log(f"total elapsed     : {_elapsed(t_total)}")

    return all_hits

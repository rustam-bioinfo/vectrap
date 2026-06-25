"""Homology-based vector sequence scanner.

This module provides the core scanning engine for VecTrap. It exposes a single
public function, ``scan()``, which runs two complementary detection strategies
against the query assembly:

1. **mappy aligner** -- for long catalog sequences (>= MIN_LEN bp).
   Uses ``mappy.Aligner`` (the official Python bindings for minimap2) against
   the pre-built ``combined_long.mmi`` index from ``vectrap-build-db``.  No
   external binary is required -- mappy is a pure pip dependency.

2. **Exact k-mer hash scanner** -- for short catalog sequences (< MIN_LEN bp).
   Loads the ``short_index.pkl`` dictionary built by ``vectrap-build-db`` and
   searches every contig for exact occurrences of each short sequence and its
   reverse complement. Returns all non-overlapping matches.

All hits from both strategies are returned as a flat list of ``HomologyHit``
objects for downstream scoring.
"""

from __future__ import annotations

import pickle
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import mappy

from vectrap.modules.utils import read_fasta, rev_comp

# Sequences shorter than this threshold are stored in the k-mer hash index;
# sequences >= this threshold are indexed by mappy.  Must match the value
# used in vectrap/cli/build_db.py when the indexes were built.
MIN_LEN: int = 50


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
) -> List[HomologyHit]:
    """Align *query_fasta* against *index_path* using mappy and return hits.

    Parameters
    ----------
    query_fasta : Path
        Input assembly FASTA (plain or gzipped).
    index_path : Path
        Pre-built minimap2 ``.mmi`` index (compatible with mappy).
    min_identity : float
        Minimum per-base sequence identity threshold (e.g. 0.90).
    min_coverage : float
        Minimum catalog sequence coverage threshold (e.g. 0.80).
    threads : int
        Number of threads for mappy (default: 4).

    Returns
    -------
    list[HomologyHit]
        Filtered hits passing both thresholds.

    Raises
    ------
    RuntimeError
        If the mappy index cannot be loaded.
    """
    aligner = mappy.Aligner(str(index_path), preset="map-ont", n_threads=threads, best_n=10)
    if not aligner:
        raise RuntimeError(
            f"mappy could not load index: {index_path}\n"
            "Run 'vectrap-build-db' to rebuild the catalog indexes."
        )

    hits: List[HomologyHit] = []
    for header, seq in read_fasta(query_fasta):
        contig_name = header.split()[0]
        for hit in aligner.map(seq, cs=True):
            # hit.ctg       -- catalog sequence name (target)
            # hit.ctg_len   -- catalog sequence length
            # hit.r_st/r_en -- target (catalog) start/end
            # hit.q_st/q_en -- query (contig) start/end, always forward-strand
            # hit.strand    -- +1 or -1

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


# ---------------------------------------------------------------------------
# Exact k-mer hash scanner (short sequences)
# ---------------------------------------------------------------------------

def _kmer_scan_sequence(
    contig_name: str,
    contig_seq: str,
    short_index: dict[str, list[str]],
) -> List[HomologyHit]:
    """Search *contig_seq* for all exact occurrences of every k-mer in *short_index*.

    Strategy
    --------
    For each catalog k-mer we perform two passes over the contig:

    **Forward pass** -- search ``contig_seq`` directly for ``kmer``.
    Matches report ``strand='+'``.

    **Reverse-complement pass** -- the RC of the contig is
    ``rc_seq = rev_comp(contig_seq)``.  A k-mer that sits on the minus strand
    of the original contig appears as ``kmer`` itself inside ``rc_seq``
    (because ``rev_comp(rev_comp(kmer)) == kmer``).  So we search ``rc_seq``
    for the *original* ``kmer`` (not ``rev_comp(kmer)``), then map the
    match position back to forward-strand coordinates::

        fwd_start = seq_len - (rc_idx + klen)
        fwd_end   = seq_len - rc_idx

    Palindromes (``rev_comp(kmer) == kmer``) are already fully captured by
    the forward pass, so the RC pass is skipped for them to avoid double
    reporting.

    Parameters
    ----------
    contig_name : str
        FASTA header of the contig being searched.
    contig_seq : str
        Nucleotide sequence of the contig (upper-case).
    short_index : dict[str, list[str]]
        Mapping of ``{sequence: [catalog_id, ...]}`` as produced by
        ``build_db.py``.

    Returns
    -------
    list[HomologyHit]
        Exact matches with ``identity=1.0`` and ``coverage=1.0``.
    """
    hits: List[HomologyHit] = []
    seq_len = len(contig_seq)
    rc_seq  = rev_comp(contig_seq)

    for kmer, catalog_ids in short_index.items():
        klen = len(kmer)
        if klen > seq_len:
            continue

        # Forward pass
        pos = 0
        while True:
            idx = contig_seq.find(kmer, pos)
            if idx == -1:
                break
            for cat_id in catalog_ids:
                hits.append(HomologyHit(
                    contig=contig_name,
                    start=idx,
                    end=idx + klen,
                    strand="+",
                    identity=1.0,
                    coverage=1.0,
                    catalog_id=cat_id,
                    source="kmer",
                ))
            pos = idx + 1

        # Reverse-complement pass -- skip palindromes
        if rev_comp(kmer) == kmer:
            continue

        pos = 0
        while True:
            idx = rc_seq.find(kmer, pos)
            if idx == -1:
                break
            fwd_start = seq_len - (idx + klen)
            fwd_end   = seq_len - idx
            for cat_id in catalog_ids:
                hits.append(HomologyHit(
                    contig=contig_name,
                    start=fwd_start,
                    end=fwd_end,
                    strand="-",
                    identity=1.0,
                    coverage=1.0,
                    catalog_id=cat_id,
                    source="kmer",
                ))
            pos = idx + 1

    return hits


def _run_kmer_scanner(
    query_fasta: Path,
    pkl_path: Path,
) -> List[HomologyHit]:
    """Load the short-sequence k-mer index and scan all contigs in *query_fasta*.

    Parameters
    ----------
    query_fasta : Path
        Input assembly FASTA (plain or gzipped).
    pkl_path : Path
        Path to the ``short_index.pkl`` file built by ``vectrap-build-db``.

    Returns
    -------
    list[HomologyHit]
        All exact k-mer matches.
    """
    with open(pkl_path, "rb") as fh:
        short_index: dict[str, list[str]] = pickle.load(fh)

    if not short_index:
        return []

    hits: List[HomologyHit] = []
    for header, seq in read_fasta(query_fasta):
        contig_name = header.split()[0]
        hits.extend(_kmer_scan_sequence(contig_name, seq, short_index))
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
        Print progress messages to stderr (default: False).

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
            print(msg, file=sys.stderr)

    _log("[vectrap] Running mappy aligner (long sequences) ...")
    mappy_hits = _run_mappy(
        query_fasta=query_fasta,
        index_path=mmi_path,
        min_identity=min_identity,
        min_coverage=min_coverage,
        threads=threads,
    )
    _log(f"[vectrap] mappy: {len(mappy_hits):,} hits after filtering.")

    _log("[vectrap] Running k-mer scanner (short sequences) ...")
    kmer_hits = _run_kmer_scanner(
        query_fasta=query_fasta,
        pkl_path=pkl_path,
    )
    _log(f"[vectrap] k-mer: {len(kmer_hits):,} hits.")

    all_hits = mappy_hits + kmer_hits
    _log(f"[vectrap] Total hits: {len(all_hits):,}.")
    return all_hits

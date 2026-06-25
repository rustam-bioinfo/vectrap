"""Homology-based vector sequence scanner.

This module provides the core scanning engine for VecTrap. It exposes a single
public function, ``scan()``, which runs two complementary detection strategies
against the query assembly:

1. **minimap2 PAF scanner** -- for long catalog sequences (>= MIN_LEN bp).
   Uses the pre-built ``combined_long.mmi`` index from ``vectrap-build-db``.
   Each alignment is filtered by identity and query coverage, then coordinates
   are normalised to 0-based forward-strand positions.

2. **Exact k-mer hash scanner** -- for short catalog sequences (< MIN_LEN bp).
   Loads the ``short_index.pkl`` dictionary built by ``vectrap-build-db`` and
   searches every contig for exact occurrences of each short sequence and its
   reverse complement. Returns all non-overlapping matches.

All hits from both strategies are returned as a flat list of ``HomologyHit``
objects for downstream scoring.
"""

from __future__ import annotations

import pickle
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

from vectrap.modules.utils import read_fasta, rev_comp

# Sequences shorter than this threshold are stored in the k-mer hash index;
# sequences >= this threshold are indexed by minimap2.  Must match the value
# used in vectrap/cli/build_db.py when the indexes were built.
MIN_LEN: int = 50

# Minimum minimap2 mapping quality to accept a hit.
_MIN_MAPQ: int = 0


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
        hits; derived from PAF divergence tag (``de:f:``) for minimap2 hits.
    coverage : float
        Fraction of the *catalog* sequence covered by the hit, in [0, 1].
        ``1.0`` for exact k-mer hits.
    catalog_id : str
        Identifier of the matching catalog entry.
    source : str
        ``'minimap2'`` or ``'kmer'`` -- indicates which scanner produced the
        hit.
    """

    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float
    catalog_id: str
    source: str = field(default="minimap2")

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

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
# minimap2 PAF scanner
# ---------------------------------------------------------------------------

_PAF_QNAME = 0
_PAF_QLEN  = 1
_PAF_QSTART = 2
_PAF_QEND   = 3
_PAF_STRAND = 4
_PAF_TNAME  = 5
_PAF_TLEN   = 6
_PAF_TSTART = 7
_PAF_TEND   = 8
_PAF_NMATCH = 9
_PAF_ALEN   = 10
_PAF_MAPQ   = 11


def _parse_paf_tags(fields: list[str]) -> dict[str, str]:
    """Parse optional PAF tags (col 12 onward) into a ``{tag: value}`` dict."""
    tags: dict[str, str] = {}
    for f in fields[12:]:
        parts = f.split(":")
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def _identity_from_tags(tags: dict[str, str], n_match: int, aln_len: int) -> float:
    """Derive per-base sequence identity from PAF optional tags.

    Preference order:
    1. ``de:f:`` tag (minimap2 gap-compressed divergence) -- most accurate.
    2. ``dv:f:`` tag (sequence divergence).
    3. Fallback: ``NM / aln_len`` approximation.
    """
    if "de" in tags:
        try:
            return max(0.0, 1.0 - float(tags["de"]))
        except ValueError:
            pass
    if "dv" in tags:
        try:
            return max(0.0, 1.0 - float(tags["dv"]))
        except ValueError:
            pass
    if aln_len > 0:
        return n_match / aln_len
    return 0.0


def _parse_paf_line(line: str, min_identity: float, min_coverage: float) -> HomologyHit | None:
    """Parse a single PAF line and return a HomologyHit or None if filtered."""
    if not line or line.startswith("#"):
        return None

    parts = line.rstrip().split("\t")
    if len(parts) < 12:
        return None

    try:
        q_name  = parts[_PAF_QNAME].split()[0]
        q_start = int(parts[_PAF_QSTART])
        q_end   = int(parts[_PAF_QEND])
        strand  = parts[_PAF_STRAND]            # '+' or '-'
        t_name  = parts[_PAF_TNAME].split()[0]
        t_len   = int(parts[_PAF_TLEN])
        t_start = int(parts[_PAF_TSTART])
        t_end   = int(parts[_PAF_TEND])
        n_match = int(parts[_PAF_NMATCH])
        aln_len = int(parts[_PAF_ALEN])
        mapq    = int(parts[_PAF_MAPQ])
    except (ValueError, IndexError):
        return None

    if mapq < _MIN_MAPQ:
        return None

    tags = _parse_paf_tags(parts)
    identity = _identity_from_tags(tags, n_match, aln_len)
    if identity < min_identity:
        return None

    # Coverage is measured on the *catalog* (target) sequence.
    covered = t_end - t_start
    coverage = covered / t_len if t_len > 0 else 0.0
    if coverage < min_coverage:
        return None

    # Coordinates in PAF are always on the query forward strand when the
    # strand is '+'.  When strand is '-' they are still query-forward, but
    # the *alignment* is to the reverse complement of the query region.
    # We therefore use q_start/q_end as-is for both strands.
    return HomologyHit(
        contig=q_name,
        start=q_start,
        end=q_end,
        strand=strand,
        identity=identity,
        coverage=coverage,
        catalog_id=t_name,
        source="minimap2",
    )


def _run_minimap2(
    query_fasta: Path,
    index_path: Path,
    min_identity: float,
    min_coverage: float,
    threads: int = 4,
) -> List[HomologyHit]:
    """Align *query_fasta* against *index_path* using minimap2 and return hits.

    Parameters
    ----------
    query_fasta : Path
        Input assembly FASTA (plain or gzipped).
    index_path : Path
        Pre-built minimap2 ``.mmi`` index.
    min_identity : float
        Minimum per-base sequence identity threshold (e.g. 0.90).
    min_coverage : float
        Minimum catalog sequence coverage threshold (e.g. 0.80).
    threads : int
        Number of minimap2 threads.

    Returns
    -------
    list[HomologyHit]
        Filtered hits passing both thresholds.

    Raises
    ------
    FileNotFoundError
        If ``minimap2`` is not available in PATH.
    RuntimeError
        If minimap2 exits with a non-zero return code.
    """
    cmd = [
        "minimap2",
        "-c",             # output CIGAR strings in PAF
        "--cs=short",     # short cs tag for divergence estimation
        "-t", str(threads),
        str(index_path),
        str(query_fasta),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            "minimap2 not found. Install it and ensure it is in your PATH. "
            "See https://github.com/lh3/minimap2"
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"minimap2 exited with code {result.returncode}.\n"
            f"stderr:\n{result.stderr}"
        )

    hits: List[HomologyHit] = []
    for line in result.stdout.splitlines():
        hit = _parse_paf_line(line, min_identity, min_coverage)
        if hit is not None:
            hits.append(hit)

    return hits


# ---------------------------------------------------------------------------
# Exact k-mer hash scanner
# ---------------------------------------------------------------------------

def _kmer_scan_sequence(
    contig_name: str,
    contig_seq: str,
    short_index: dict[str, list[str]],
) -> List[HomologyHit]:
    """Search *contig_seq* for all exact occurrences of every k-mer in *short_index*.

    Both the forward sequence and a reverse complement pass are performed.
    Overlapping matches for the same k-mer are all reported (non-deduplicated);
    deduplication of overlapping hits from different catalog entries is left to
    the scorer.

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
    rc_seq = rev_comp(contig_seq)

    for kmer, catalog_ids in short_index.items():
        klen = len(kmer)
        if klen > seq_len:
            continue

        # Forward strand search
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

        # Reverse complement search -- convert rc coordinates back to
        # forward strand: fwd_start = seq_len - rc_end
        rc_kmer = rev_comp(kmer)
        # Avoid double-counting palindromes
        if rc_kmer == kmer:
            continue

        pos = 0
        while True:
            idx = rc_seq.find(rc_kmer, pos)
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

    This function runs both the minimap2 PAF scanner (long sequences) and the
    exact k-mer hash scanner (short sequences) and returns their combined
    output as a flat list of ``HomologyHit`` objects.

    Parameters
    ----------
    query_fasta : str or Path
        Input assembly FASTA file (plain or gzipped).
    catalog_dir : str or Path
        Directory containing the indexes built by ``vectrap-build-db``:
        ``combined_long.mmi`` and ``short_index.pkl``.
    min_identity : float, optional
        Minimum per-base sequence identity for minimap2 hits (default: 0.90).
        Exact k-mer hits always have identity 1.0 regardless of this setting.
    min_coverage : float, optional
        Minimum fraction of the catalog sequence that must be covered by a hit
        (default: 0.80). Exact k-mer hits always have coverage 1.0.
    threads : int, optional
        Number of threads passed to minimap2 (default: 4).
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
        If the catalog index files are missing or ``minimap2`` is not in PATH.
    RuntimeError
        If minimap2 exits with a non-zero return code.
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

    _log("[vectrap] Running minimap2 scanner (long sequences) ...")
    mm2_hits = _run_minimap2(
        query_fasta=query_fasta,
        index_path=mmi_path,
        min_identity=min_identity,
        min_coverage=min_coverage,
        threads=threads,
    )
    _log(f"[vectrap] minimap2: {len(mm2_hits):,} hits after filtering.")

    _log("[vectrap] Running k-mer scanner (short sequences) ...")
    kmer_hits = _run_kmer_scanner(
        query_fasta=query_fasta,
        pkl_path=pkl_path,
    )
    _log(f"[vectrap] k-mer: {len(kmer_hits):,} hits.")

    all_hits = mm2_hits + kmer_hits
    _log(f"[vectrap] Total hits: {len(all_hits):,}.")
    return all_hits

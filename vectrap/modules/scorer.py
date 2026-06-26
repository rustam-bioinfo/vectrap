"""Per-contig hit summarisation for VecTrap.

This module aggregates enriched ``HomologyHit`` objects (produced by
``homology_scanner.scan()``) into per-contig statistics without applying any
classification logic.  The summary is written to ``verdicts.tsv`` by
``vectrap.cli.run``.

Public API
----------
    summarize(hits) -> List[ContigSummary]
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import List

from vectrap.modules.homology_scanner import HomologyHit


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------

@dataclass
class ContigSummary:
    """Hit statistics for a single query contig.

    Attributes
    ----------
    contig : str
        Contig name (FASTA header first token).
    total_hits : int
        Total number of hits (mappy + kmer, before overlap merging).
    engineered_hits : int
        Hits with tier == 'ENGINEERED'.
    context_dependent_hits : int
        Hits with tier == 'CONTEXT_DEPENDENT'.
    weak_hits : int
        Hits with tier == 'WEAK'.
    unannotated_hits : int
        Hits with no tier annotation (empty string).
    unique_feature_types : int
        Number of distinct feature_type values among all hits.
    unique_labels : int
        Number of distinct label values among all hits.
    covered_bp : int
        Total non-overlapping base pairs on the contig covered by hits.
    top_label : str
        Label of the hit with the highest identity*coverage score.
    top_tier : str
        Tier of the top hit.
    top_confidence : str
        Confidence of the top hit.
    evidence_summary : str
        Compact semicolon-separated list of unique labels with their tier.
    """

    contig: str
    total_hits: int = 0
    engineered_hits: int = 0
    context_dependent_hits: int = 0
    weak_hits: int = 0
    unannotated_hits: int = 0
    unique_feature_types: int = 0
    unique_labels: int = 0
    covered_bp: int = 0
    top_label: str = ""
    top_tier: str = ""
    top_confidence: str = ""
    evidence_summary: str = ""


# ---------------------------------------------------------------------------
# Overlap merge helper
# ---------------------------------------------------------------------------

def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent intervals and return sorted non-overlapping list.

    Parameters
    ----------
    intervals : list of (start, end) tuples
        0-based half-open intervals.

    Returns
    -------
    list of (start, end) tuples
        Merged, sorted intervals.
    """
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _covered_bp(hits: list[HomologyHit]) -> int:
    """Return the total number of non-overlapping bases covered by *hits*."""
    intervals = [(h.start, h.end) for h in hits]
    return sum(e - s for s, e in _merge_intervals(intervals))


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def summarize(hits: List[HomologyHit]) -> List[ContigSummary]:
    """Aggregate enriched hits into per-contig statistics.

    Parameters
    ----------
    hits : list[HomologyHit]
        All hits returned by ``homology_scanner.scan()``, already enriched
        with tier/confidence/reasoning.

    Returns
    -------
    list[ContigSummary]
        One entry per contig that had at least one hit, sorted by
        (engineered_hits DESC, context_dependent_hits DESC, contig ASC).
    """
    # Group hits by contig
    by_contig: dict[str, list[HomologyHit]] = defaultdict(list)
    for hit in hits:
        by_contig[hit.contig].append(hit)

    summaries: list[ContigSummary] = []

    for contig, chits in by_contig.items():
        s = ContigSummary(contig=contig)
        s.total_hits = len(chits)

        for h in chits:
            t = (h.tier or "").upper()
            if t == "ENGINEERED":
                s.engineered_hits += 1
            elif t == "CONTEXT_DEPENDENT":
                s.context_dependent_hits += 1
            elif t == "WEAK":
                s.weak_hits += 1
            else:
                s.unannotated_hits += 1

        s.unique_feature_types = len({h.feature_type for h in chits if h.feature_type})
        s.unique_labels        = len({h.label for h in chits if h.label})
        s.covered_bp           = _covered_bp(chits)

        # Top hit by identity * coverage
        top = max(chits, key=lambda h: h.identity * h.coverage)
        s.top_label      = top.label
        s.top_tier       = top.tier
        s.top_confidence = top.confidence

        # Evidence summary: group labels by tier, show up to 5 per tier
        tier_order = ["ENGINEERED", "CONTEXT_DEPENDENT", "WEAK", ""]
        tier_labels: dict[str, list[str]] = defaultdict(list)
        seen: set[tuple[str, str]] = set()
        for h in chits:
            key = (h.tier or "", h.label or h.catalog_id)
            if key not in seen:
                seen.add(key)
                tier_labels[h.tier or ""].append(h.label or h.catalog_id)

        parts = []
        tier_abbrev = {
            "ENGINEERED": "ENG",
            "CONTEXT_DEPENDENT": "CD",
            "WEAK": "WEAK",
            "": "?",
        }
        for tier in tier_order:
            labels = tier_labels.get(tier)
            if not labels:
                continue
            abbr = tier_abbrev.get(tier, tier)
            shown = labels[:5]
            suffix = f"+{len(labels) - 5} more" if len(labels) > 5 else ""
            label_str = ", ".join(shown)
            if suffix:
                label_str += f", {suffix}"
            parts.append(f"{abbr}[{len(labels)}]: {label_str}")
        s.evidence_summary = "; ".join(parts)

        summaries.append(s)

    # Sort: most engineered first, then most context-dependent, then by name
    summaries.sort(key=lambda s: (-s.engineered_hits, -s.context_dependent_hits, s.contig))
    return summaries

"""Sample-level summary report for VecTrap.

Aggregates per-contig statistics produced by ``scorer.summarize()`` into a
single row per genome and writes ``vectrap_summary.tsv`` to the output
directory.

Public API
----------
    build_sample_row(sample, hits, summaries, contig_lengths,
                     min_engineered_contamination, min_context_suspected)
        -> dict

    write_summary(rows, output_dir) -> Path
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List

from vectrap.modules.homology_scanner import HomologyHit
from vectrap.modules.scorer import ContigSummary, _covered_bp


_SUMMARY_FIELDS = [
    "sample",
    "total_contigs",
    "total_bp",
    "total_hits",
    "contigs_with_hits",
    "engineered_hits",
    "context_dependent_hits",
    "weak_hits",
    "sample_verdict",
    "top_labels",
    "unique_feature_types",
    "unique_labels",
    "covered_bp",
    "covered_fraction",
    "top_contig",
    "top_contig_length",
    "top_contig_covered_fraction",
]


def _sample_verdict(
    engineered_hits: int,
    context_dependent_hits: int,
    min_engineered_contamination: int,
    min_context_suspected: int,
) -> str:
    if engineered_hits >= min_engineered_contamination:
        return "CONTAMINATION"
    if engineered_hits > 0 or context_dependent_hits >= min_context_suspected:
        return "SUSPECTED"
    return "CLEAN"


def build_sample_row(
    sample: str,
    hits: List[HomologyHit],
    summaries: List[ContigSummary],
    contig_lengths: Dict[str, int],
    min_engineered_contamination: int = 3,
    min_context_suspected: int = 1,
) -> dict:
    """Build one summary row for *sample*.

    Parameters
    ----------
    sample : str
        Sample stem name.
    hits : list[HomologyHit]
        All enriched hits for this sample.
    summaries : list[ContigSummary]
        Per-contig summaries produced by ``scorer.summarize()``.
    contig_lengths : dict[str, int]
        Mapping contig name -> length in bp.
    min_engineered_contamination : int
        Engineered-hit count threshold for CONTAMINATION verdict (default: 3).
    min_context_suspected : int
        Context-dependent-hit count threshold for SUSPECTED verdict (default: 1).
    """
    total_contigs = len(contig_lengths)
    total_bp      = sum(contig_lengths.values())
    total_hits    = len(hits)
    contigs_with  = len(summaries)

    eng  = sum(s.engineered_hits          for s in summaries)
    cd   = sum(s.context_dependent_hits   for s in summaries)
    weak = sum(s.weak_hits                for s in summaries)

    verdict = _sample_verdict(eng, cd, min_engineered_contamination, min_context_suspected)

    # top 5 labels by hit count
    label_counts = Counter(h.label for h in hits if h.label)
    top_labels = ", ".join(lbl for lbl, _ in label_counts.most_common(5))

    unique_ft     = len({h.feature_type for h in hits if h.feature_type})
    unique_labels = len({h.label        for h in hits if h.label})

    covered = _covered_bp(hits)
    covered_frac = covered / total_bp if total_bp > 0 else 0.0

    # contig with highest covered_fraction
    top_s = max(summaries, key=lambda s: s.covered_fraction) if summaries else None
    top_contig              = top_s.contig            if top_s else ""
    top_contig_length       = top_s.contig_length     if top_s else 0
    top_contig_covered_frac = top_s.covered_fraction  if top_s else 0.0

    return {
        "sample":                       sample,
        "total_contigs":                total_contigs,
        "total_bp":                     total_bp,
        "total_hits":                   total_hits,
        "contigs_with_hits":            contigs_with,
        "engineered_hits":              eng,
        "context_dependent_hits":       cd,
        "weak_hits":                    weak,
        "sample_verdict":               verdict,
        "top_labels":                   top_labels,
        "unique_feature_types":         unique_ft,
        "unique_labels":                unique_labels,
        "covered_bp":                   covered,
        "covered_fraction":             f"{covered_frac:.6f}",
        "top_contig":                   top_contig,
        "top_contig_length":            top_contig_length,
        "top_contig_covered_fraction":  f"{top_contig_covered_frac:.6f}",
    }


def write_summary(rows: list[dict], output_dir: Path) -> Path:
    """Write *rows* to ``vectrap_summary.tsv`` in *output_dir*.

    Parameters
    ----------
    rows : list[dict]
        One dict per sample, as returned by ``build_sample_row()``.
    output_dir : Path
        Directory where the file will be written.

    Returns
    -------
    Path
        Absolute path to the written file.
    """
    out_path = output_dir / "vectrap_summary.tsv"
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return out_path

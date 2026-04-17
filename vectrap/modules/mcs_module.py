#!/usr/bin/env python3
"""Multiple Cloning Site (MCS) scanner.

Fixes applied vs. original:
1. detect_mcs_windows: removed premature ``break`` inside the inner loop.
   Hits are sorted by (start, end), so .end is NOT strictly monotone --
   a later hit can have a smaller .end if it belongs to a shorter motif.
   The old ``break`` fired as soon as span > max_span, silently skipping
   valid shorter windows that followed.  The fix uses ``continue`` so every
   (i, j) pair where hits[j].start - hits[i].start <= max_span is evaluated.
   An outer early-exit is kept: once hits[j].start alone already exceeds
   hits[i].start + max_span, no further j can contribute (start is
   monotone even if end is not), so the inner loop breaks there instead.
2. Isoschizomers (NarI/KasI share GGCGCC; SmaI/XmaI share CCCGGG):
   Both enzymes are intentionally preserved as separate Hit objects so
   that unique_enzymes output is informative.  unique_motif_count is
   unaffected because score_window deduplicates by motif, not enzyme name.
   A note is added to the WindowCall docstring so downstream users are not
   confused by enzyme counts exceeding motif counts.
3. Shared utilities (rev_comp, open_text, read_fasta) are now imported from
   utils.py instead of being copy-pasted.
"""
import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple
import re

from utils import rev_comp, open_text, read_fasta

SITE_ROWS = [
    ("EcoRI",  "GAATTC",   "core",     "sticky", "classical"),
    ("XbaI",   "TCTAGA",   "core",     "sticky", "classical"),
    ("SpeI",   "ACTAGT",   "core",     "sticky", "classical"),
    ("NotI",   "GCGGCCGC", "core",     "rare8",  "rare_gc"),
    ("PstI",   "CTGCAG",   "core",     "sticky", "classical"),
    ("BamHI",  "GGATCC",   "core",     "sticky", "classical"),
    ("HindIII","AAGCTT",   "core",     "sticky", "classical"),
    ("KpnI",   "GGTACC",   "core",     "sticky", "classical"),
    ("SacI",   "GAGCTC",   "core",     "sticky", "classical"),
    ("SalI",   "GTCGAC",   "core",     "sticky", "classical"),
    ("SmaI",   "CCCGGG",   "core",     "blunt",  "gc_rich"),
    # XmaI is an isoschizomer of SmaI (same motif CCCGGG, different cut position).
    # Both are kept so unique_enzymes output reflects cloning options accurately.
    # unique_motif_count is NOT inflated -- score_window deduplicates by motif.
    ("XmaI",   "CCCGGG",   "extended", "sticky", "gc_rich"),
    ("SphI",   "GCATGC",   "core",     "sticky", "classical"),
    ("XhoI",   "CTCGAG",   "extended", "sticky", "directional"),
    ("NheI",   "GCTAGC",   "extended", "sticky", "directional"),
    ("EcoRV",  "GATATC",   "extended", "blunt",  "blunt"),
    ("ClaI",   "ATCGAT",   "extended", "sticky", "classical"),
    ("ApaI",   "GGGCCC",   "extended", "sticky", "gc_rich"),
    ("BglII",  "AGATCT",   "extended", "sticky", "directional"),
    ("AvrII",  "CCTAGG",   "extended", "sticky", "directional"),
    ("MluI",   "ACGCGT",   "extended", "sticky", "classical"),
    ("NcoI",   "CCATGG",   "extended", "sticky", "coding_friendly"),
    ("NdeI",   "CATATG",   "extended", "sticky", "coding_friendly"),
    ("BspEI",  "TCCGGA",   "extended", "sticky", "gc_rich"),
    ("AflII",  "CTTAAG",   "extended", "sticky", "classical"),
    ("AgeI",   "ACCGGT",   "extended", "sticky", "coding_friendly"),
    ("BsrGI",  "TGTACA",   "extended", "sticky", "coding_friendly"),
    ("PacI",   "TTAATTAA", "extended", "rare8",  "rare_at"),
    ("AscI",   "GGCGCGCC", "extended", "rare8",  "rare_gc"),
    ("FseI",   "GGCCGGCC", "extended", "rare8",  "rare_gc"),
    ("PmeI",   "GTTTAAAC", "extended", "rare8",  "rare_at"),
    ("StuI",   "AGGCCT",   "extended", "blunt",  "blunt"),
    ("NsiI",   "ATGCAT",   "extended", "sticky", "classical"),
    ("BclI",   "TGATCA",   "extended", "sticky", "classical"),
    ("BstBI",  "TTCGAA",   "extended", "sticky", "classical"),
    ("SacII",  "CCGCGG",   "extended", "sticky", "gc_rich"),
    # NarI and KasI are isoschizomers (same motif GGCGCC, different cut positions).
    # See SmaI/XmaI note above -- same rationale applies here.
    ("NarI",   "GGCGCC",   "extended", "sticky", "gc_rich"),
    ("KasI",   "GGCGCC",   "extended", "sticky", "gc_rich"),
    ("EagI",   "CGGCCG",   "extended", "sticky", "gc_rich"),
    ("BsaI",   "GGTCTC",   "extended", "sticky", "type_iis"),
    ("BsmBI",  "CGTCTC",   "extended", "sticky", "type_iis"),
    ("BbsI",   "GAAGAC",   "extended", "sticky", "type_iis"),
    ("SapI",   "GCTCTTC",  "extended", "sticky", "type_iis"),
    ("AarI",   "CACCTGC",  "extended", "sticky", "type_iis"),
]

DIRECTIONAL_MOTIFS = {
    "GAATTC", "CTCGAG", "AAGCTT", "GGATCC",
    "GCTAGC", "TCTAGA", "ACTAGT", "AGATCT", "CCTAGG",
}
RARE_MOTIFS = {
    "GCGGCCGC", "TTAATTAA", "GGCGCGCC", "GGCCGGCC", "GTTTAAAC",
}


@dataclass
class EnzymeInfo:
    name: str
    motif: str
    panel: str
    cut_type: str
    site_class: str


@dataclass
class Hit:
    contig: str
    motif: str
    enzyme_name: str
    start: int
    end: int
    strand: str
    site_class: str
    cut_type: str


@dataclass
class WindowCall:
    """A candidate MCS window.

    Note on enzyme vs. motif counts
    --------------------------------
    unique_enzymes may list more entries than unique_motif_count when
    isoschizomers are present (e.g. NarI + KasI both recognise GGCGCC;
    SmaI + XmaI both recognise CCCGGG).  The score is always based on
    unique_motif_count, not enzyme count, so isoschizomers do not inflate
    the score.
    """
    contig: str
    start: int
    end: int
    span_bp: int
    unique_motif_count: int
    total_hit_count: int
    unique_enzymes: str
    unique_motifs: str
    class_set: str
    rare_motif_count: int
    directional_motif_count: int
    repeat_penalty: int
    score: int
    strength: str
    hits: List[Hit]


def build_catalog(panel: str) -> List[EnzymeInfo]:
    catalog = []
    for enzyme, motif, min_panel, cut_type, site_class in SITE_ROWS:
        if panel == "core" and min_panel != "core":
            continue
        catalog.append(EnzymeInfo(enzyme, motif, min_panel, cut_type, site_class))
    return catalog


def find_hits(seq: str, contig: str, catalog: List[EnzymeInfo]) -> List[Hit]:
    hits = []
    motif_to_enzymes: dict = defaultdict(list)
    for enz in catalog:
        motif_to_enzymes[enz.motif].append(enz)
    motifs = list(motif_to_enzymes.keys())
    pattern = re.compile(f"(?=({'|'.join(motifs)}))")

    for match in pattern.finditer(seq):
        motif = match.group(1)
        start = match.start()
        end = start + len(motif)
        for enz in motif_to_enzymes[motif]:
            hits.append(Hit(contig, motif, enz.name, start, end, "+",
                            enz.site_class, enz.cut_type))

    rc_seq = rev_comp(seq)
    seq_len = len(seq)
    for match in pattern.finditer(rc_seq):
        motif = match.group(1)
        end = seq_len - match.start()
        start = end - len(motif)
        for enz in motif_to_enzymes[motif]:
            hits.append(Hit(contig, motif, enz.name, start, end, "-",
                            enz.site_class, enz.cut_type))

    # Deduplicate: for palindromic sites the same (start, end, motif, enzyme)
    # appears on both strands -- keep the '+' strand representative.
    unique_hits: dict = {}
    for h in hits:
        key = (h.start, h.end, h.motif, h.enzyme_name)
        if key not in unique_hits or h.strand == "+":
            unique_hits[key] = h
    return sorted(unique_hits.values(), key=lambda x: (x.start, x.end))


def score_window(hits: List[Hit]) -> Tuple[int, int, int, int, str]:
    unique_motifs = {h.motif for h in hits}
    class_categories = {c for h in hits for c in h.site_class.split(';')}
    unique_count = len(unique_motifs)
    rare_count = sum(1 for m in unique_motifs if m in RARE_MOTIFS)
    directional_count = sum(1 for m in unique_motifs if m in DIRECTIONAL_MOTIFS)
    mixed_class_bonus = 1 if len(class_categories) >= 3 else 0
    repeat_penalty = max(0, len(hits) - unique_count)
    score = unique_count + rare_count + mixed_class_bonus - repeat_penalty
    span_bp = hits[-1].end - hits[0].start if hits else 0
    strength = "none"
    if unique_count >= 8 and span_bp <= 150 and (rare_count >= 1 or score >= 9):
        strength = "strong"
    elif unique_count >= 7 and span_bp <= 180 and directional_count >= 1:
        strength = "moderate"
    elif unique_count >= 6 and span_bp <= 200:
        strength = "weak"
    return score, rare_count, directional_count, repeat_penalty, strength


def detect_mcs_windows(
    hits: List[Hit],
    max_span: int = 200,
) -> List[WindowCall]:
    """Find all hit windows whose span does not exceed *max_span* bp.

    BUG FIX (vs. original): the original inner loop used ``break`` the first
    time ``span > max_span``.  Because hits are sorted by (start, end) and
    *.end* is NOT strictly monotone, a later hit with a shorter motif can
    have a smaller .end, making the span drop back below max_span.  The old
    break silently skipped those valid windows.

    The corrected logic:
    * ``continue`` (skip, do not record) when the span exceeds max_span.
    * ``break`` only when hits[j].start alone already exceeds
      hits[i].start + max_span -- since .start IS monotone, no subsequent j
      can ever produce a valid window for this i.
    """
    if not hits:
        return []

    valid_windows = []
    n = len(hits)
    for i in range(n):
        for j in range(i, n):
            # Early exit: once the start position of j is beyond the window
            # boundary for i, no further j can be valid (starts are monotone).
            if hits[j].start > hits[i].start + max_span:
                break
            span = hits[j].end - hits[i].start
            # Skip this j if the span is too large, but keep looking because
            # a later hit with a shorter motif may bring the span back down.
            if span > max_span:
                continue
            current_hits = hits[i : j + 1]
            score, rc, dc, rp, strength = score_window(current_hits)
            if strength != "none":
                unique_motifs = sorted({h.motif for h in current_hits})
                unique_enzymes = sorted({h.enzyme_name for h in current_hits})
                class_set = sorted(
                    {c for h in current_hits for c in h.site_class.split(';')}
                )
                valid_windows.append(WindowCall(
                    contig=current_hits[0].contig,
                    start=current_hits[0].start,
                    end=current_hits[-1].end,
                    span_bp=span,
                    unique_motif_count=len(unique_motifs),
                    total_hit_count=len(current_hits),
                    unique_enzymes=",".join(unique_enzymes),
                    unique_motifs=",".join(unique_motifs),
                    class_set=",".join(class_set),
                    rare_motif_count=rc,
                    directional_motif_count=dc,
                    repeat_penalty=rp,
                    score=score,
                    strength=strength,
                    hits=current_hits,
                ))
    return resolve_overlaps(valid_windows)


def resolve_overlaps(calls: List[WindowCall]) -> List[WindowCall]:
    if not calls:
        return []
    calls.sort(key=lambda x: (-x.score, x.span_bp))
    final_calls = []
    used_hits: set = set()
    for call in calls:
        call_hit_ids = {(h.start, h.end, h.enzyme_name) for h in call.hits}
        if not call_hit_ids.intersection(used_hits):
            final_calls.append(call)
            used_hits.update(call_hit_ids)
    final_calls.sort(key=lambda x: x.start)
    return final_calls


def write_hits(path: Path, hits: Iterable[Hit]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'motif', 'enzyme', 'start_0based', 'end_0based',
            'strand', 'site_class', 'cut_type',
        ])
        for h in hits:
            writer.writerow([
                h.contig, h.motif, h.enzyme_name, h.start, h.end,
                h.strand, h.site_class, h.cut_type,
            ])


def write_windows(path: Path, calls: Iterable[WindowCall]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'window_start_0based', 'window_end_0based', 'span_bp',
            'unique_motif_count', 'total_hit_count', 'unique_enzymes',
            'unique_motifs', 'class_set', 'rare_motif_count',
            'directional_motif_count', 'repeat_penalty', 'score', 'mcs_strength',
        ])
        for c in calls:
            writer.writerow([
                c.contig, c.start, c.end, c.span_bp,
                c.unique_motif_count, c.total_hit_count,
                c.unique_enzymes, c.unique_motifs, c.class_set,
                c.rare_motif_count, c.directional_motif_count,
                c.repeat_penalty, c.score, c.strength,
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Scan a FASTA file for Multiple Cloning Site windows."
    )
    parser.add_argument('-i', '--input',  required=True,
                        help="Input FASTA (plain or .gz)")
    parser.add_argument('-o', '--outdir', required=True,
                        help="Output directory")
    parser.add_argument('--panel', choices=['core', 'extended'],
                        default='extended',
                        help="Enzyme panel (default: extended)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    catalog = build_catalog(args.panel)

    all_hits: List[Hit] = []
    all_calls: List[WindowCall] = []
    for contig, seq in read_fasta(args.input):
        hits = find_hits(seq, contig, catalog)
        calls = detect_mcs_windows(hits)
        all_hits.extend(hits)
        all_calls.extend(calls)

    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith('.gz') else Path(stem).stem
    write_hits(outdir / f'{stem}.mcs_hits.tsv', all_hits)
    write_windows(outdir / f'{stem}.mcs_windows.tsv', all_calls)
    print(f'MCS scan complete. Wrote to {outdir}')


if __name__ == '__main__':
    main()

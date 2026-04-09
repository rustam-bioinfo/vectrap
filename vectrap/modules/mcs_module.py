#!/usr/bin/env python3
import argparse
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple
import re

SITE_ROWS = [
    ("EcoRI", "GAATTC", "core", "sticky", "classical"),
    ("XbaI", "TCTAGA", "core", "sticky", "classical"),
    ("SpeI", "ACTAGT", "core", "sticky", "classical"),
    ("NotI", "GCGGCCGC", "core", "rare8", "rare_gc"),
    ("PstI", "CTGCAG", "core", "sticky", "classical"),
    ("BamHI", "GGATCC", "core", "sticky", "classical"),
    ("HindIII", "AAGCTT", "core", "sticky", "classical"),
    ("KpnI", "GGTACC", "core", "sticky", "classical"),
    ("SacI", "GAGCTC", "core", "sticky", "classical"),
    ("SalI", "GTCGAC", "core", "sticky", "classical"),
    ("SmaI", "CCCGGG", "core", "blunt", "gc_rich"),
    ("XmaI", "CCCGGG", "extended", "sticky", "gc_rich"),
    ("SphI", "GCATGC", "core", "sticky", "classical"),
    ("XhoI", "CTCGAG", "extended", "sticky", "directional"),
    ("NheI", "GCTAGC", "extended", "sticky", "directional"),
    ("EcoRV", "GATATC", "extended", "blunt", "blunt"),
    ("ClaI", "ATCGAT", "extended", "sticky", "classical"),
    ("ApaI", "GGGCCC", "extended", "sticky", "gc_rich"),
    ("BglII", "AGATCT", "extended", "sticky", "directional"),
    ("AvrII", "CCTAGG", "extended", "sticky", "directional"),
    ("MluI", "ACGCGT", "extended", "sticky", "classical"),
    ("NcoI", "CCATGG", "extended", "sticky", "coding_friendly"),
    ("NdeI", "CATATG", "extended", "sticky", "coding_friendly"),
    ("BspEI", "TCCGGA", "extended", "sticky", "gc_rich"),
    ("AflII", "CTTAAG", "extended", "sticky", "classical"),
    ("AgeI", "ACCGGT", "extended", "sticky", "coding_friendly"),
    ("BsrGI", "TGTACA", "extended", "sticky", "coding_friendly"),
    ("PacI", "TTAATTAA", "extended", "rare8", "rare_at"),
    ("AscI", "GGCGCGCC", "extended", "rare8", "rare_gc"),
    ("FseI", "GGCCGGCC", "extended", "rare8", "rare_gc"),
    ("PmeI", "GTTTAAAC", "extended", "rare8", "rare_at"),
    ("StuI", "AGGCCT", "extended", "blunt", "blunt"),
    ("NsiI", "ATGCAT", "extended", "sticky", "classical"),
    ("BclI", "TGATCA", "extended", "sticky", "classical"),
    ("BstBI", "TTCGAA", "extended", "sticky", "classical"),
    ("SacII", "CCGCGG", "extended", "sticky", "gc_rich"),
    ("NarI", "GGCGCC", "extended", "sticky", "gc_rich"),
    ("KasI", "GGCGCC", "extended", "sticky", "gc_rich"),
    ("EagI", "CGGCCG", "extended", "sticky", "gc_rich"),
    ("BsaI", "GGTCTC", "extended", "sticky", "type_iis"),
    ("BsmBI", "CGTCTC", "extended", "sticky", "type_iis"),
    ("BbsI", "GAAGAC", "extended", "sticky", "type_iis"),
    ("SapI", "GCTCTTC", "extended", "sticky", "type_iis"),
    ("AarI", "CACCTGC", "extended", "sticky", "type_iis"),
]

DIRECTIONAL_MOTIFS = {"GAATTC", "CTCGAG", "AAGCTT", "GGATCC", "GCTAGC", "TCTAGA", "ACTAGT", "AGATCT", "CCTAGG"}
RARE_MOTIFS = {"GCGGCCGC", "TTAATTAA", "GGCGCGCC", "GGCCGGCC", "GTTTAAAC"}

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

def rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def build_catalog(panel: str) -> List[EnzymeInfo]:
    catalog = []
    for enzyme, motif, min_panel, cut_type, site_class in SITE_ROWS:
        if panel == "core" and min_panel != "core":
            continue
        catalog.append(EnzymeInfo(enzyme, motif, min_panel, cut_type, site_class))
    return catalog

def find_hits(seq: str, contig: str, catalog: List[EnzymeInfo]) -> List[Hit]:
    hits = []
    motif_to_enzymes = defaultdict(list)
    for enz in catalog:
        motif_to_enzymes[enz.motif].append(enz)
    motifs = list(motif_to_enzymes.keys())
    pattern = re.compile(f"(?=({'|'.join(motifs)}))")
    for match in pattern.finditer(seq):
        motif = match.group(1)
        start = match.start()
        end = start + len(motif)
        for enz in motif_to_enzymes[motif]:
            hits.append(Hit(contig, motif, enz.name, start, end, "+", enz.site_class, enz.cut_type))
    rc_seq = rev_comp(seq)
    seq_len = len(seq)
    for match in pattern.finditer(rc_seq):
        motif = match.group(1)
        end = seq_len - match.start()
        start = end - len(motif)
        for enz in motif_to_enzymes[motif]:
            hits.append(Hit(contig, motif, enz.name, start, end, "-", enz.site_class, enz.cut_type))
    unique_hits = {}
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

def detect_mcs_windows(hits: List[Hit], max_span: int = 200) -> List[WindowCall]:
    if not hits:
        return []
    valid_windows = []
    n = len(hits)
    for i in range(n):
        for j in range(i, n):
            span = hits[j].end - hits[i].start
            if span > max_span:
                break
            current_hits = hits[i:j+1]
            score, rc, dc, rp, strength = score_window(current_hits)
            if strength != "none":
                unique_motifs = sorted({h.motif for h in current_hits})
                unique_enzymes = sorted({h.enzyme_name for h in current_hits})
                class_set = sorted({c for h in current_hits for c in h.site_class.split(';')})
                valid_windows.append(WindowCall(
                    contig=current_hits[0].contig, start=current_hits[0].start,
                    end=current_hits[-1].end, span_bp=span,
                    unique_motif_count=len(unique_motifs), total_hit_count=len(current_hits),
                    unique_enzymes=",".join(unique_enzymes), unique_motifs=",".join(unique_motifs),
                    class_set=",".join(class_set), rare_motif_count=rc,
                    directional_motif_count=dc, repeat_penalty=rp, score=score,
                    strength=strength, hits=current_hits
                ))
    return resolve_overlaps(valid_windows)

def resolve_overlaps(calls: List[WindowCall]) -> List[WindowCall]:
    if not calls:
        return []
    calls.sort(key=lambda x: (-x.score, x.span_bp))
    final_calls = []
    used_hits = set()
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
        writer.writerow(['contig', 'motif', 'enzyme', 'start_0based', 'end_0based', 'strand', 'site_class', 'cut_type'])
        for h in hits:
            writer.writerow([h.contig, h.motif, h.enzyme_name, h.start, h.end, h.strand, h.site_class, h.cut_type])

def write_windows(path: Path, calls: Iterable[WindowCall]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'window_start_0based', 'window_end_0based', 'span_bp', 'unique_motif_count',
            'total_hit_count', 'unique_enzymes', 'unique_motifs', 'class_set', 'rare_motif_count',
            'directional_motif_count', 'repeat_penalty', 'score', 'mcs_strength'
        ])
        for c in calls:
            writer.writerow([
                c.contig, c.start, c.end, c.span_bp, c.unique_motif_count, c.total_hit_count,
                c.unique_enzymes, c.unique_motifs, c.class_set, c.rare_motif_count,
                c.directional_motif_count, c.repeat_penalty, c.score, c.strength
            ])

def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        name = None
        chunks = []
        for line in handle:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None: yield name, "".join(chunks).upper()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('--panel', choices=['core', 'extended'], default='extended')
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    catalog = build_catalog(args.panel)
    all_hits = []
    all_calls = []
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

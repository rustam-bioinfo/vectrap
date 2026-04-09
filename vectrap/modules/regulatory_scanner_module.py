#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple

REGULATORY_CATALOG = {
    "tac_promoter": (r'(?=((TTGACA)(.{16})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco', 'PRIMARY'),
    "trc_promoter": (r'(?=((TTGACA)(.{17})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco', 'PRIMARY'),
    "T5_lacO_promoter": (r'(?=((TTGCTT)(.{17})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco', 'PRIMARY'),
    "lacUV5_promoter": (r'(?=((TTTACA)(.{18})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco', 'PRIMARY'),
    "T7_promoter_core": (r'(?=((TAATACGACTCACTATA)))', 'monolithic', 'SUPPORTIVE'),
    "T3_promoter_core": (r'(?=((AATTAACCCTCACTAAA)))', 'monolithic', 'SUPPORTIVE'),
    "SP6_promoter_core": (r'(?=((ATTTAGGTGACACTATAG)))', 'monolithic', 'SUPPORTIVE'),
    "lacWT_promoter": (r'(?=((TTTACA)(.{18})(TATGTT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco', 'SUPPORTIVE'),
    "lambda_pL_promoter": (r'(?=((TTGACA)(.{17})(GATACT)))', 'bipartite', 'SUPPORTIVE'),
    "T7_terminator": (r'(?=((GCTAGTTATTGCTCAGCGG)))', 'monolithic', 'SUPPORTIVE'),
    "VF2_primer": (r'(?=((TGCCACCTGACGTCTAAGAA)))', 'monolithic', 'PRIMARY'),
    "VR_primer": (r'(?=((ATTACCGCCTTTGAGTGAGC)))', 'monolithic', 'PRIMARY'),
    "M13_fwd_21": (r'(?=((TGTAAAACGACGGCCAGT)))', 'monolithic', 'PRIMARY'),
    "M13_fwd_40": (r'(?=((GTTTTCCCAGTCACGAC)))', 'monolithic', 'PRIMARY'),
    "M13_rev_27": (r'(?=((CAGGAAACAGCTATGAC)))', 'monolithic', 'PRIMARY'),
    "M13_pUC_rev_48": (r'(?=((AGCGGATAACAATTTCACACAGG)))', 'monolithic', 'PRIMARY'),
    "RFC10_prefix_NC": (r'(?=((GAATTCGCGGCCGCTTCTAGAG)))', 'monolithic', 'PRIMARY'),
    "RFC10_prefix_CDS": (r'(?=((GAATTCGCGGCCGCTTCTAG)))', 'monolithic', 'PRIMARY'),
    "RFC10_suffix": (r'(?=((TACTAGTAGCGGCCGCTGCAG)))', 'monolithic', 'PRIMARY'),
    "RFC10_scar": (r'(?=((TACTAGAG)))', 'monolithic', 'CONTEXT_DEPENDENT'),
}

COMPILED_CATALOG = {
    name: (re.compile(pattern, re.IGNORECASE), arch, ev_class)
    for name, (pattern, arch, ev_class) in REGULATORY_CATALOG.items()
}

@dataclass
class RegulatoryHit:
    contig: str
    start_0based: int
    end_0based: int
    strand: str
    marker_name: str
    evidence_class: str
    architecture: str
    matched_sequence: str
    minus_35: str = "N/A"
    minus_10: str = "N/A"
    spacer_len: str = "N/A"
    lac_operator: str = "N/A"

def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        name = None
        chunks: List[str] = []
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

def rev_comp(seq: str) -> str:
    # N->N is IUPAC-correct; lowercase handled defensively in case called outside read_fasta
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def scan_regulatory_elements(seq: str, contig: str) -> List[RegulatoryHit]:
    raw_hits = []
    seq_len = len(seq)
    rc_seq = rev_comp(seq)

    for name, (pattern, arch, ev_class) in COMPILED_CATALOG.items():
        for strand, target_seq in [("+", seq), ("-", rc_seq)]:
            for match in pattern.finditer(target_seq):
                full_seq = match.group(1)
                m35, m10, laco, sp_len = "N/A", "N/A", "N/A", "N/A"
                if arch in ['bipartite', 'bipartite_laco']:
                    m35 = match.group(2)
                    sp_len = str(len(match.group(3)))
                    m10 = match.group(4)
                if arch == 'bipartite_laco':
                    laco = match.group(6)
                if strand == "+":
                    start = match.start()
                    end = start + len(full_seq)
                else:
                    start = seq_len - match.start() - len(full_seq)
                    end = seq_len - match.start()
                raw_hits.append(RegulatoryHit(
                    contig=contig, start_0based=start, end_0based=end, strand=strand,
                    marker_name=name, evidence_class=ev_class, architecture=arch,
                    matched_sequence=full_seq, minus_35=m35, minus_10=m10,
                    spacer_len=sp_len, lac_operator=laco
                ))

    valid_biobrick_anchors = {"RFC10_prefix_NC", "RFC10_prefix_CDS", "RFC10_suffix"}
    found_markers = {h.marker_name for h in raw_hits}
    has_biobrick_context = bool(valid_biobrick_anchors.intersection(found_markers))

    final_hits = []
    for h in raw_hits:
        if h.marker_name == "RFC10_scar":
            if has_biobrick_context:
                h.evidence_class = "PRIMARY"
                final_hits.append(h)
        else:
            final_hits.append(h)
    return sorted(final_hits, key=lambda x: x.start_0based)

def write_hits(path: Path, hits: Iterable[RegulatoryHit]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'start_0based', 'end_0based', 'strand', 'marker_name',
            'evidence_class', 'architecture', 'minus_35_seq', 'spacer_len',
            'minus_10_seq', 'lac_operator_seq', 'matched_sequence'
        ])
        for h in hits:
            writer.writerow([
                h.contig, h.start_0based, h.end_0based, h.strand, h.marker_name,
                h.evidence_class, h.architecture, h.minus_35, h.spacer_len,
                h.minus_10, h.lac_operator, h.matched_sequence
            ])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outdir', required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_hits = []
    for contig, seq in read_fasta(args.input):
        all_hits.extend(scan_regulatory_elements(seq, contig))
    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith('.gz') else Path(stem).stem
    outfile = outdir / f'{stem}.regulatory_markers.tsv'
    write_hits(outfile, all_hits)
    print(f'Found {len(all_hits)} regulatory markers. Wrote to {outfile}')

if __name__ == '__main__':
    main()

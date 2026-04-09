#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple

PROMOTER_CATALOG = {
    "tac": (r'(?=((TTGACA)(.{16})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco'),
    "trc": (r'(?=((TTGACA)(.{17})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco'),
    "lacUV5": (r'(?=((TTTACA)(.{18})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco'),
    "lacWT": (r'(?=((TTTACA)(.{18})(TATGTT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco'),
    "T5_lacO": (r'(?=((TTGCTT)(.{17})(TATAAT)(.{4,12})(TGAGCGGATAACAA)))', 'bipartite_laco'),
    "lambda_pL": (r'(?=((TTGACA)(.{17})(GATACT)))', 'bipartite'),
    "T7": (r'(?=((TAATACGACTCACTATAGGG)))', 'monolithic'),
    "T3": (r'(?=((AATTAACCCTCACTAAAGGG)))', 'monolithic'),
    "SP6": (r'(?=((ATTTAGGTGACACTATAGAA)))', 'monolithic'),
}

COMPILED_CATALOG = {
    name: (re.compile(pattern, re.IGNORECASE), arch)
    for name, (pattern, arch) in PROMOTER_CATALOG.items()
}

@dataclass
class PromoterHit:
    contig: str
    start: int
    end: int
    strand: str
    promoter_name: str
    minus_35: str
    minus_10: str
    spacer_len: str
    lac_operator: str
    matched_seq: str

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
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def find_promoters(seq: str, contig: str) -> List[PromoterHit]:
    hits = []
    seq_len = len(seq)
    rc_seq = rev_comp(seq)
    for name, (pattern, arch) in COMPILED_CATALOG.items():
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
                    rc_start = match.start()
                    rc_end = rc_start + len(full_seq)
                    end = seq_len - rc_start
                    start = seq_len - rc_end
                hits.append(PromoterHit(
                    contig=contig, start=start, end=end, strand=strand,
                    promoter_name=name, minus_35=m35, minus_10=m10,
                    spacer_len=sp_len, lac_operator=laco, matched_seq=full_seq
                ))
    return sorted(hits, key=lambda x: x.start)

def write_hits(path: Path, hits: Iterable[PromoterHit]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'start_0based', 'end_0based', 'strand', 'promoter_type',
            'minus_35_seq', 'spacer_len', 'minus_10_seq', 'lac_operator_seq', 'full_match'
        ])
        for h in hits:
            writer.writerow([
                h.contig, h.start, h.end, h.strand, h.promoter_name,
                h.minus_35, h.spacer_len, h.minus_10, h.lac_operator, h.matched_seq
            ])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outdir', required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_hits: List[PromoterHit] = []
    for contig, seq in read_fasta(args.input):
        all_hits.extend(find_promoters(seq, contig))
    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith('.gz') else Path(stem).stem
    outfile = outdir / f'{stem}.promoter_hits.tsv'
    write_hits(outfile, all_hits)
    print(f'Found {len(all_hits)} engineered promoters. Wrote to {outfile}')

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

# --- Configuration: Origin Sequence Anchors ---
# We use highly conserved structural regions (RNA transcripts or DNA iterons) 
# instead of whole-origin alignment to guarantee detection even in heavily mutated vectors.
ORIGIN_CATALOG = {
# 1. The ColE1 / pMB1 Family
"pUC_ori_HighCopy": r"GTTCTCAGTAATTC",
"pMB1_pBR322_ori_MediumCopy": r"GTTCTGAGTAATTC",

# 2. Compatible Cloning Origins
"p15A_ori_LowCopy": r"CTGACGAGTGGAAGGAGGACG",
# Increased required iteron tandem repeats to 3 for absolute specificity
"pSC101_ori_LowCopy": r"(?:TGACAGGT){3,}",

# 3. Broad-Host-Range Origins
"RK2_oriV_IncP_BroadHost": r"TGACGCACCC(?=.{0,20}TGACGCACCC)",
# Increased required iteron tandem repeats to 3
"RSF1010_IncQ_BroadHost": r"(?:CCTTGCA){3,}",
"pBBR1_BroadHost": r"GATCGCCCGC(?=.{0,30}GATCGCCCGC)",

# 4. Phage / Single-Stranded Origins
"f1_M13_Phagemid_ori": r"TGCAAACTATTAACTGGCGAACTACT",

# 5. Eukaryotic Shuttle Origins (STRICT anchors)
# The 11bp ARS core [AT]TTTAT[AG]TTT[AT] followed tightly by the 2-micron specific poly-A tract.
# This prevents random 11bp AT-rich genomic noise from triggering a hit.
"Yeast_2u_Shuttle": r"[AT]TTTAT[AG]TTT[AT].{5,15}AAAATAAAA",

# SV40 Core Origin (27bp perfect inverted repeat)
"Mammalian_SV40_Shuttle": r"GCCTCGGCCTCTGCATAAATAAAAAAA"
}

# Pre-compile the regex dictionary for maximum execution speed
COMPILED_ORIGINS = {
name: re.compile(pattern, re.IGNORECASE) for name, pattern in ORIGIN_CATALOG.items()
}

@dataclass
class OriginHit:
    contig: str
    start_0based: int
    end_0based: int
    strand: str
    origin_type: str
    matched_sequence: str

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
    return seq.translate(str.maketrans('ATCGNacgtn', 'TAGCNtgcan'))[::-1]

def scan_origins(seq: str, contig: str) -> List[OriginHit]:
    hits = []
    seq_len = len(seq)
    rc_seq = rev_comp(seq)

    for name, pattern in COMPILED_ORIGINS.items():
        for match in pattern.finditer(seq):
            hits.append(OriginHit(
                contig=contig,
                start_0based=match.start(),
                end_0based=match.end(),
                strand="+",
                origin_type=name,
                matched_sequence=match.group(0)
            ))
        for match in pattern.finditer(rc_seq):
            fwd_start = seq_len - match.end()
            fwd_end = seq_len - match.start()
            hits.append(OriginHit(
                contig=contig,
                start_0based=fwd_start,
                end_0based=fwd_end,
                strand="-",
                origin_type=name,
                matched_sequence=match.group(0)
            ))
    return sorted(hits, key=lambda x: x.start_0based)

def write_hits(path: Path, hits: Iterable[OriginHit]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['contig', 'start_0based', 'end_0based', 'strand', 'origin_type', 'matched_sequence'])
        for h in hits:
            writer.writerow([h.contig, h.start_0based, h.end_0based, h.strand, h.origin_type, h.matched_sequence])

def main():
    parser = argparse.ArgumentParser(description='Ultra-fast nucleotide scanner for synthetic Origins of Replication.')
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outdir', required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_hits: List[OriginHit] = []
    for contig, seq in read_fasta(args.input):
        all_hits.extend(scan_origins(seq, contig))
    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith('.gz') else Path(stem).stem
    outfile = outdir / f'{stem}.origin_hits.tsv'
    write_hits(outfile, all_hits)
    print(f'Found {len(all_hits)} engineered origins. Wrote to {outfile}')

if __name__ == '__main__':
    main()

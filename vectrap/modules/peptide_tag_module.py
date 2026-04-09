#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

MOTIF_CATALOG = {
    "His_tag": r"H{6,}",
    "FLAG_tag": r"DYKDDDDK",
    "Strep_tag_II": r"WSHPQFEK",
    "TEV_cleavage": r"ENLYFQ[GS]",
    "Thrombin_cleavage": r"LVPRGS",
    "Flexible_linker_GGS": r"(?:G{3,4}S){2,}",
    "Rigid_linker_EAAAK": r"(?:EAAAK){2,}",
    "Myc_tag": r"EQKLISEEDL",
    "HA_tag": r"YPYDVPDYA",
    "V5_tag": r"GKPIPNPLLGLDST",
}

COMPILED_MOTIFS = {name: re.compile(pattern) for name, pattern in MOTIF_CATALOG.items()}

@dataclass
class FastORF:
    contig: str
    start: int
    end: int
    strand: str
    aa_seq: str

@dataclass
class ProteinMotifHit:
    contig: str
    orf_start: int
    orf_end: int
    strand: str
    tag_start: int
    tag_end: int
    motif_type: str
    matched_aa: str
    dist_n_term: int
    dist_c_term: int
    orf_len_aa: int

def translate_dna(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i+3], 'X') for i in range(0, len(seq) - 2, 3))

def rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans('ATCGNacgtn', 'TAGCNtgcan'))[::-1]

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

def find_orfs(seq: str, contig: str, min_aa: int = 30) -> List[FastORF]:
    orfs = []
    seq_len = len(seq)
    for frame in range(3):
        aa_seq = translate_dna(seq[frame:])
        peptides = aa_seq.split('*')
        pos_aa = 0
        for pep in peptides:
            if len(pep) >= min_aa:
                nt_start = frame + (pos_aa * 3)
                nt_end = nt_start + (len(pep) * 3)
                orfs.append(FastORF(contig, nt_start, nt_end, '+', pep))
            pos_aa += len(pep) + 1
    rev_seq = rev_comp(seq)
    for frame in range(3):
        aa_seq = translate_dna(rev_seq[frame:])
        peptides = aa_seq.split('*')
        pos_aa = 0
        for pep in peptides:
            if len(pep) >= min_aa:
                rev_nt_start = frame + (pos_aa * 3)
                rev_nt_end = rev_nt_start + (len(pep) * 3)
                fwd_start = seq_len - rev_nt_end
                fwd_end = seq_len - rev_nt_start
                orfs.append(FastORF(contig, fwd_start, fwd_end, '-', pep))
            pos_aa += len(pep) + 1
    return orfs

def scan_peptide_motifs(orf: FastORF) -> List[ProteinMotifHit]:
    hits = []
    orf_aa_len = len(orf.aa_seq)
    for motif_name, pattern in COMPILED_MOTIFS.items():
        for match in pattern.finditer(orf.aa_seq):
            start_aa = match.start()
            end_aa = match.end()
            if orf.strand == '+':
                tag_start = orf.start + (start_aa * 3)
                tag_end = orf.start + (end_aa * 3)
            else:
                tag_end = orf.end - (start_aa * 3)
                tag_start = orf.end - (end_aa * 3)
            hits.append(ProteinMotifHit(
                contig=orf.contig, orf_start=orf.start, orf_end=orf.end, strand=orf.strand,
                tag_start=tag_start, tag_end=tag_end, motif_type=motif_name,
                matched_aa=match.group(0), dist_n_term=start_aa,
                dist_c_term=orf_aa_len - end_aa, orf_len_aa=orf_aa_len
            ))
    return hits

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('--min-aa', type=int, default=30)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_hits: List[ProteinMotifHit] = []
    for contig, seq in read_fasta(args.input):
        orfs = find_orfs(seq, contig, args.min_aa)
        for orf in orfs:
            all_hits.extend(scan_peptide_motifs(orf))
    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith('.gz') else Path(stem).stem
    outfile = outdir / f'{stem}.peptide_tags.tsv'
    with open(outfile, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'contig', 'orf_start_0based', 'orf_end_0based', 'strand',
            'tag_start_0based', 'tag_end_0based', 'motif_type', 'matched_aa',
            'dist_n_term_aa', 'dist_c_term_aa', 'orf_len_aa'
        ])
        for h in all_hits:
            writer.writerow([
                h.contig, h.orf_start, h.orf_end, h.strand, h.tag_start, h.tag_end,
                h.motif_type, h.matched_aa, h.dist_n_term, h.dist_c_term, h.orf_len_aa
            ])
    print(f'Found {len(all_hits)} peptide motif hits. Wrote to {outfile}')

if __name__ == '__main__':
    main()

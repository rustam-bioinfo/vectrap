#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Tuple

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
    "His_tag":              r"H{6,}",
    "FLAG_tag":             r"DYKDDDDK",
    "Strep_tag_II":         r"WSHPQFEK",
    "TEV_cleavage":         r"ENLYFQ[GS]",
    "Thrombin_cleavage":    r"LVPRGS",
    "Flexible_linker_GGS":  r"(?:G{3,4}S){2,}",
    "Rigid_linker_EAAAK":   r"(?:EAAAK){2,}",
    "Myc_tag":              r"EQKLISEEDL",
    "HA_tag":               r"YPYDVPDYA",
    "V5_tag":               r"GKPIPNPLLGLDST",
}

COMPILED_MOTIFS = {name: re.compile(pattern) for name, pattern in MOTIF_CATALOG.items()}

# --------------------------------------------------------------------------- #
# Data structures                                                              #
# --------------------------------------------------------------------------- #

@dataclass
class ProteinMotifHit:
    contig:          str
    tag_start_0based: int   # absolute DNA coord, 0-based inclusive
    tag_end_0based:   int   # absolute DNA coord, 0-based exclusive
    strand:          str
    motif_type:      str
    matched_aa:      str
    orf_start_0based: int   # absolute DNA coord, 0-based inclusive
    orf_end_0based:   int   # absolute DNA coord, 0-based exclusive
    dist_to_start_aa: int   # AA distance from ORF start (M or frame start) to tag
    dist_to_stop_aa:  int   # AA distance from tag end to downstream stop
    has_start_codon:  bool

# --------------------------------------------------------------------------- #
# Sequence utilities                                                           #
# --------------------------------------------------------------------------- #

_COMP_TABLE = str.maketrans("ATCGNacgtn", "TAGCNtgcan")

def rev_comp(seq: str) -> str:
    return seq.translate(_COMP_TABLE)[::-1]

def translate_frame(seq: str, frame: int) -> str:
    """Translate a DNA sequence starting at the given frame offset."""
    s = seq[frame:]
    return "".join(CODON_TABLE.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3))

def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        name = None
        chunks: List[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks).upper()

# --------------------------------------------------------------------------- #
# Seed-and-extend core                                                         #
# --------------------------------------------------------------------------- #

def _aa_to_dna_fwd(frame: int, aa_pos: int) -> int:
    """Convert an amino acid index in a forward frame to an absolute 0-based
    DNA coordinate (start of the corresponding codon)."""
    return frame + aa_pos * 3

def _aa_to_dna_rev(seq_len: int, frame: int, aa_pos: int) -> int:
    """Convert an amino acid index in a reverse-complement frame to an absolute
    0-based DNA coordinate on the *forward* strand.

    The reverse-complement sequence has length seq_len.  After slicing off
    `frame` leading nucleotides the translated region covers nucleotides
    [frame, seq_len).  Amino acid aa_pos maps to rev-comp position
    frame + aa_pos*3 (start of codon).  Converting back to the forward
    strand: fwd_pos = seq_len - 1 - rev_pos for the last nt of that codon,
    which equals seq_len - frame - aa_pos*3 - 3 (0-based start, exclusive
    end = seq_len - frame - aa_pos*3).

    Returns the 0-based *start* of the codon on the forward strand (i.e. the
    smaller coordinate), and the *end* (exclusive) is +3 away."""
    rev_codon_start = frame + aa_pos * 3
    fwd_end_exclusive = seq_len - rev_codon_start        # exclusive
    fwd_start = fwd_end_exclusive - 3                    # inclusive
    return fwd_start  # caller adds 3 for the exclusive end

def scan_sequence(contig: str, seq: str) -> List[ProteinMotifHit]:
    """Seed-and-extend: translate all 6 frames, seed from motif matches,
    extend to flanking stop codons, find nearest upstream start codon."""
    hits: List[ProteinMotifHit] = []
    seq_len = len(seq)
    rev = rev_comp(seq)

    for strand, dna in (("+", seq), ("-", rev)):
        for frame in range(3):
            aa_seq = translate_frame(dna, frame)
            aa_len = len(aa_seq)

            for motif_name, pattern in COMPILED_MOTIFS.items():
                for m in pattern.finditer(aa_seq):
                    tag_aa_start = m.start()   # inclusive
                    tag_aa_end   = m.end()     # exclusive

                    # ---- extend downstream to the next stop codon ----------
                    stop_pos = aa_seq.find("*", tag_aa_end)
                    if stop_pos == -1:
                        stop_pos = aa_len      # no stop; use sequence end
                    dist_to_stop_aa = stop_pos - tag_aa_end

                    # ---- extend upstream to the nearest preceding stop -----
                    upstream_stop = aa_seq.rfind("*", 0, tag_aa_start)
                    # upstream_stop == -1 means no stop; frame starts at 0
                    frame_start_aa = upstream_stop + 1  # aa after stop (or 0)

                    # ---- find nearest start codon (M) between frame_start
                    #      and the tag -----------------------------------------
                    segment = aa_seq[frame_start_aa:tag_aa_start]
                    m_offset = segment.rfind("M")  # last M before the tag
                    if m_offset != -1:
                        orf_start_aa = frame_start_aa + m_offset
                        has_start_codon = True
                    else:
                        orf_start_aa = frame_start_aa
                        has_start_codon = False

                    dist_to_start_aa = tag_aa_start - orf_start_aa

                    # ---- convert AA coords to absolute DNA coords -----------
                    if strand == "+":
                        tag_dna_start = _aa_to_dna_fwd(frame, tag_aa_start)
                        tag_dna_end   = _aa_to_dna_fwd(frame, tag_aa_end)
                        orf_dna_start = _aa_to_dna_fwd(frame, orf_start_aa)
                        orf_dna_end   = _aa_to_dna_fwd(frame, stop_pos) + (
                            3 if stop_pos < aa_len else 0
                        )
                    else:
                        # For the reverse strand the coordinates are computed
                        # on the rev-comp string and then mirrored back.
                        # _aa_to_dna_rev returns the 0-based fwd start of the
                        # codon; the exclusive end is that value + 3.
                        #
                        # The TAG interval on the fwd strand is:
                        #   [fwd_pos(tag_aa_end - 1),  fwd_pos(tag_aa_start) + 3)
                        # i.e. the last codon's fwd_start to the first codon's
                        # fwd_end.  Because on the reverse strand higher AA
                        # index -> lower fwd coordinate:
                        tag_dna_start = _aa_to_dna_rev(seq_len, frame, tag_aa_end - 1)
                        tag_dna_end   = _aa_to_dna_rev(seq_len, frame, tag_aa_start) + 3

                        orf_dna_start = _aa_to_dna_rev(seq_len, frame, stop_pos - 1) if stop_pos < aa_len else (
                            seq_len - frame - (aa_len * 3)
                        )
                        orf_dna_end   = _aa_to_dna_rev(seq_len, frame, orf_start_aa) + 3

                        # guard: ensure start <= end
                        if orf_dna_start > orf_dna_end:
                            orf_dna_start, orf_dna_end = orf_dna_end, orf_dna_start
                        if tag_dna_start > tag_dna_end:
                            tag_dna_start, tag_dna_end = tag_dna_end, tag_dna_start

                    # clamp to sequence bounds
                    tag_dna_start = max(0, tag_dna_start)
                    tag_dna_end   = min(seq_len, tag_dna_end)
                    orf_dna_start = max(0, orf_dna_start)
                    orf_dna_end   = min(seq_len, orf_dna_end)

                    hits.append(ProteinMotifHit(
                        contig=contig,
                        tag_start_0based=tag_dna_start,
                        tag_end_0based=tag_dna_end,
                        strand=strand,
                        motif_type=motif_name,
                        matched_aa=m.group(0),
                        orf_start_0based=orf_dna_start,
                        orf_end_0based=orf_dna_end,
                        dist_to_start_aa=dist_to_start_aa,
                        dist_to_stop_aa=dist_to_stop_aa,
                        has_start_codon=has_start_codon,
                    ))
    return hits

# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Detect peptide/affinity tags in nucleotide sequences "
                    "using a seed-and-extend approach across all 6 reading frames."
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="Input FASTA (plain or .gz)")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_hits: List[ProteinMotifHit] = []
    for contig, seq in read_fasta(args.input):
        all_hits.extend(scan_sequence(contig, seq))

    stem = Path(args.input).name
    stem = Path(stem[:-3]).stem if stem.endswith(".gz") else Path(stem).stem
    outfile = outdir / f"{stem}.peptide_tags.tsv"

    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "contig",
            "tag_start_0based",
            "tag_end_0based",
            "strand",
            "motif_type",
            "matched_aa",
            "orf_start_0based",
            "orf_end_0based",
            "dist_to_start_aa",
            "dist_to_stop_aa",
            "has_start_codon",
        ])
        for h in all_hits:
            writer.writerow([
                h.contig,
                h.tag_start_0based,
                h.tag_end_0based,
                h.strand,
                h.motif_type,
                h.matched_aa,
                h.orf_start_0based,
                h.orf_end_0based,
                h.dist_to_start_aa,
                h.dist_to_stop_aa,
                h.has_start_codon,
            ])

    print(f"Found {len(all_hits)} peptide motif hits. Wrote to {outfile}")

if __name__ == "__main__":
    main()

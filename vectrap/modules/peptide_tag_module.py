#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

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

# All sequences verified against primary literature.
# See vectrap/docs/peptide_tag_references.txt for full citations.
MOTIF_CATALOG = {
    # ---- affinity / purification -------------------------------------------
    # His-tag: 6-9 consecutive histidines (Hochuli et al. 1988, J Chromatogr)
    "His_tag":               r"H{6,9}",
    # His-tag 10x: 10 or more histidines; higher Ni-NTA affinity tier
    "His_tag_10x":           r"H{10,}",
    # FLAG: Hopp et al. 1988, Biotechnology 6:1204
    "FLAG_tag":              r"DYKDDDDK",
    # 3xFLAG: DYKDHDG-DYKDHDI-DYKDDDDK (Sigma-Aldrich / Brizzard et al.)
    "FLAG_tag_3x":           r"DYKDHDGDYKDHDIDYKDDDDK",
    # Strep-tag I: Schmidt & Skerra 1993, Protein Eng 6:109
    "Strep_tag_I":           r"WRHPQFGG",
    # Strep-tag II: Schmidt & Skerra 1994 / 2007, Nat Methods
    "Strep_tag_II":          r"WSHPQFEK",
    # Twin-Strep-tag: two Strep-II units joined by a linker (IBA Lifesciences)
    "Twin_Strep_tag":        r"WSHPQFEK.{1,10}WSHPQFEK",
    # SpyTag: Zakeri et al. 2012, PNAS; confirmed sequence AHIVMVDAYKPTK
    "SpyTag":                r"AHIVMVDAYKPTK",
    # ALFA-tag: Gotzke et al. 2019, Nat Commun 10:4403; sequence SRLEEELRRRLTE
    "ALFA_tag":              r"SRLEEELRRRLTE",
    # AviTag: biotin-acceptor peptide, BirA substrate; Beckett et al. 1999
    "AviTag":                r"GLNDIFEAQKIEWHE",
    # SBP-tag: streptavidin-binding peptide; Keefe et al. 2001, Protein Expr Purif
    "SBP_tag":               r"MDEKTTGWRGGHVVEGLAGELEQLRARLEHHPQGQREP",
    # Softag3: mild-elution affinity tag used in tandem-affinity and Y2H systems
    "Softag3":               r"TQDPSRVVGQLEQRPPR",
    # ---- epitope tags -------------------------------------------------------
    # Myc-tag: Evans et al. 1985, Mol Cell Biol 5:3610
    "Myc_tag":               r"EQKLISEEDL",
    # HA-tag: Wilson et al. 1984, Cell 37:767
    "HA_tag":                r"YPYDVPDYA",
    # V5-tag: Southern et al. 1991, J Gen Virol; sequence GKPIPNPLLGLDST
    "V5_tag":                r"GKPIPNPLLGLDST",
    # T7-tag: N-terminal tag from T7 gene 10 leader; Novagen
    "T7_tag":                r"MASMTGGQQMG",
    # E-tag: 13-AA synthetic tag; Pharmacia/GE Healthcare
    "E_tag":                 r"GAPVPYPDPLEPR",
    # VSV-G tag: from vesicular stomatitis virus glycoprotein; Bhatt et al.
    "VSV_G_tag":             r"YTDIEMNRLGK",
    # AU1-tag: from BPV-1 major capsid protein; sequence DTYRYI (Luo et al. 1993)
    "AU1_tag":               r"DTYRYI",
    # AU5-tag: from BPV-1 major capsid protein; sequence TDFYLK (Luo et al. 1993)
    "AU5_tag":               r"TDFYLK",
    # Protein C tag: human protein C epitope; Ca2+-dependent elution
    "Protein_C_tag":         r"EDQVDPRLIDGK",
    # ---- protease cleavage sites --------------------------------------------
    # TEV: Parks et al. 1994, Anal Biochem; consensus ENLYFQ[GS]
    "TEV_cleavage":          r"ENLYFQ[GS]",
    # Thrombin: LVPRGS (Sigma / GE Healthcare standard)
    "Thrombin_cleavage":     r"LVPRGS",
    # PreScission (HRV 3C): LEVLFQGP; Walker et al. 1994 / GE Healthcare
    "PreScission_cleavage":  r"LEVLFQGP",
    # Enterokinase: DDDDK; cleavage site immediately downstream of FLAG
    "Enterokinase_cleavage": r"DDDDK",
    # Factor Xa: IEGR; cleavage after R (Nagai & Thogersen 1987)
    "Factor_Xa_cleavage":    r"IEG[R]",
    # ---- linkers ------------------------------------------------------------
    # Flexible GGS linker: (GGGS)n or (GGGGS)n repeats; Chen et al. 2013
    "Flexible_linker_GGS":   r"(?:G{3,4}S){2,}",
    # Rigid EAAAK linker: (EAAAK)n alpha-helical repeats; Arai et al. 2001
    "Rigid_linker_EAAAK":    r"(?:EAAAK){2,}",
    # Rigid PAPAP linker: proline-alanine repeats; semi-rigid spacer
    "Rigid_linker_PAPAP":    r"(?:PAPAP){2,}",
    # Helical poly-alanine linker: (AAAA)n; used in synthetic constructs
    "Helical_linker_A4":     r"(?:AAAA){3,}",
}

COMPILED_MOTIFS = {name: re.compile(pattern) for name, pattern in MOTIF_CATALOG.items()}

# --------------------------------------------------------------------------- #
# Data structures                                                              #
# --------------------------------------------------------------------------- #

@dataclass
class ProteinMotifHit:
    contig:           str
    tag_start_0based: int   # absolute DNA coord, 0-based inclusive
    tag_end_0based:   int   # absolute DNA coord, 0-based exclusive
    strand:           str
    motif_type:       str
    matched_aa:       str
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
    """Translate *seq* starting at *frame* (0, 1, or 2).

    Uses idiomatic (len // 3) * 3 range so the last incomplete codon is
    always excluded cleanly rather than relying on the fragile len-2 sentinel.
    """
    s = seq[frame:]
    return "".join(CODON_TABLE.get(s[i:i + 3], "X") for i in range(0, (len(s) // 3) * 3, 3))

def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        name = None
        chunks: list[str] = []
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
# Coordinate helpers                                                           #
# --------------------------------------------------------------------------- #

def _aa_to_dna_fwd(frame: int, aa_pos: int) -> int:
    """AA index in a forward frame -> absolute 0-based DNA start of that codon."""
    return frame + aa_pos * 3

def _aa_to_dna_rev(seq_len: int, frame: int, aa_pos: int) -> int:
    """AA index in a reverse-complement frame -> absolute 0-based DNA start of
    the corresponding codon on the *forward* strand.

    On the rev-comp string the codon begins at rev_codon_start = frame + aa_pos*3.
    Mirroring back: the codon occupies forward positions
        [seq_len - rev_codon_start - 3,  seq_len - rev_codon_start)
    so this function returns the inclusive start; the exclusive end is +3.
    """
    rev_codon_start = frame + aa_pos * 3
    return seq_len - rev_codon_start - 3

# --------------------------------------------------------------------------- #
# Seed-and-extend core                                                         #
# --------------------------------------------------------------------------- #

def scan_sequence(contig: str, seq: str) -> list[ProteinMotifHit]:
    """Translate all 6 reading frames, seed from compiled motif matches, then
    extend each seed to flanking stop codons and find the nearest upstream M.

    Coordinate contract
    -------------------
    Forward strand  : orf_start <= tag_start < tag_end <= orf_end  (all fwd)
    Reverse strand  : orf_start <= tag_start < tag_end <= orf_end  (all fwd)
                      orf_start is the 5'-most nt of the stop codon (low coord)
                      orf_end   is the 3'-most nt+1 of the M codon  (high coord)
    """
    if len(seq) < 3:
        warnings.warn(
            f"Contig '{contig}' is shorter than 3 nt ({len(seq)} nt) "
            "and cannot be translated. Skipping.",
            RuntimeWarning,
            stacklevel=2,
        )
        return []

    hits: list[ProteinMotifHit] = []
    seen: set[tuple[str, str, int, str]] = set()

    seq_len = len(seq)
    rev = rev_comp(seq)

    for strand, dna in (("+", seq), ("-", rev)):
        for frame in range(3):
            aa_seq = translate_frame(dna, frame)
            aa_len = len(aa_seq)

            for motif_name, pattern in COMPILED_MOTIFS.items():
                for m in pattern.finditer(aa_seq):
                    tag_aa_start = m.start()
                    tag_aa_end   = m.end()

                    stop_pos = aa_seq.find("*", tag_aa_end)
                    if stop_pos == -1:
                        stop_pos = aa_len
                    dist_to_stop_aa = stop_pos - tag_aa_end

                    upstream_stop  = aa_seq.rfind("*", 0, tag_aa_start)
                    frame_start_aa = upstream_stop + 1

                    segment  = aa_seq[frame_start_aa:tag_aa_start]
                    m_offset = segment.rfind("M")
                    if m_offset != -1:
                        orf_start_aa    = frame_start_aa + m_offset
                        has_start_codon = True
                    else:
                        orf_start_aa    = frame_start_aa
                        has_start_codon = False

                    dist_to_start_aa = tag_aa_start - orf_start_aa

                    if strand == "+":
                        tag_dna_start = _aa_to_dna_fwd(frame, tag_aa_start)
                        tag_dna_end   = _aa_to_dna_fwd(frame, tag_aa_end)
                        orf_dna_start = _aa_to_dna_fwd(frame, orf_start_aa)
                        stop_dna      = _aa_to_dna_fwd(frame, stop_pos)
                        orf_dna_end   = stop_dna + 3 if stop_pos < aa_len else stop_dna

                    else:
                        tag_dna_start = _aa_to_dna_rev(seq_len, frame, tag_aa_end - 1)
                        tag_dna_end   = _aa_to_dna_rev(seq_len, frame, tag_aa_start) + 3

                        if stop_pos < aa_len:
                            orf_dna_start = _aa_to_dna_rev(seq_len, frame, stop_pos)
                        else:
                            orf_dna_start = seq_len - frame - aa_len * 3

                        orf_dna_end = _aa_to_dna_rev(seq_len, frame, orf_start_aa) + 3

                        assert orf_dna_start <= tag_dna_start, (
                            f"[{contig} {strand} frame {frame}] "
                            f"orf_dna_start {orf_dna_start} > tag_dna_start {tag_dna_start}"
                        )
                        assert tag_dna_start <= tag_dna_end, (
                            f"[{contig} {strand} frame {frame}] "
                            f"tag_dna_start {tag_dna_start} > tag_dna_end {tag_dna_end}"
                        )
                        assert tag_dna_end <= orf_dna_end, (
                            f"[{contig} {strand} frame {frame}] "
                            f"tag_dna_end {tag_dna_end} > orf_dna_end {orf_dna_end}"
                        )

                    tag_dna_start = max(0, tag_dna_start)
                    tag_dna_end   = min(seq_len, tag_dna_end)
                    orf_dna_start = max(0, orf_dna_start)
                    orf_dna_end   = min(seq_len, orf_dna_end)

                    dedup_key = (contig, strand, tag_dna_start, motif_name)
                    if dedup_key in seen:
                        continue
                    seen.add(dedup_key)

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

    all_hits: list[ProteinMotifHit] = []
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

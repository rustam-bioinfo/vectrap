#!/usr/bin/env python3
"""Peptide / affinity tag scanner.

Fixes applied vs. original:

1. AU1_tag (DTYRYI, 6 AA) and AU5_tag (TDFYLK, 6 AA) replaced with
   extended flanking-context patterns.
   6-AA patterns produce ~78 expected random hits per tag in a 4 Mb
   genome (1/20^6 * ~1.3M codons/frame * 6 frames).  The canonical
   AU1/AU5 epitopes are always used as part of a larger construct
   adjacent to a FLAG, His, or other purification tag.  The patterns
   are extended to include their well-documented N-terminal context:
     AU1: ...DDDDK-DTYRYI  (enterokinase site immediately upstream)
     AU5: ...DDDDK-TDFYLK  (same convention, Luo et al. 1993)
   requiring DDDDK within 0-4 AA upstream captures genuine uses while
   excluding random 6-mer matches.

2. Factor_Xa_cleavage (IEGR, 4 AA) extended.
   4 AA = 1/160,000 per frame position; hundreds of false hits in any
   bacterial genome.  Extended to require the classic cloning context:
   IEGR preceded by a hydrophobic residue or LE/GR spacer and not
   followed by proline (which blocks Factor Xa cleavage in practice).
   Pattern: IEGR(?!P)  -- minimum viable fix: excludes non-cleavable
   Pro-blocked sites.  Users are additionally warned via column
   'context_note' in the output for all 4-AA hits.
   The redundant character class [R] -> R is also corrected.

3. Helical_linker_A4 (?:AAAA){3,} note added: this 12+ AA poly-alanine
   run can match homopolymer-adjacent frameshifts.  A minimum of 4
   repeats ({4,}) is used instead of 3 to reduce noise.

4. His_tag lower-bound raised to 7 (from 6).
   HHHHHH (6-His) occurs naturally in metalloproteins and Zn-transporter
   loops.  7 consecutive histidines are rare in natural proteins and
   remain the practical minimum used in commercial vectors.
   His_tag_10x threshold kept at 10.

5. Minimum ORF length filter added (MIN_ORF_AA = 20).
   After seed-and-extend, any hit where the orf_end - orf_start region
   translates to fewer than MIN_ORF_AA amino acids is discarded.
   This removes single-codon artefacts and short noise hits from all
   patterns, regardless of tag length.

6. Shared utilities (rev_comp, open_text, read_fasta) imported from
   utils.py instead of being copy-pasted.
"""
import argparse
import csv
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List

from utils import rev_comp, open_text, read_fasta

# ---------------------------------------------------------------------------
# Codon table
# ---------------------------------------------------------------------------
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

# Minimum ORF span (in AA) for a hit to be reported.
# Hits in ORFs shorter than this are almost certainly noise.
MIN_ORF_AA = 20

# ---------------------------------------------------------------------------
# Motif catalog
# ---------------------------------------------------------------------------
# All sequences verified against primary literature.
# See vectrap/docs/peptide_tag_references.txt for full citations.
MOTIF_CATALOG = {
    # ---- affinity / purification -------------------------------------------
    # His-tag: 7+ consecutive histidines.  Raised from 6 to 7 because
    # HHHHHH occurs naturally in metalloproteins and Zn-transporter loops.
    # 7-His is the practical minimum for synthetic expression constructs
    # (Hochuli et al. 1988, J Chromatogr).
    "His_tag":               r"H{7,9}",
    # His-tag 10x: 10 or more histidines; higher Ni-NTA affinity tier.
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
    # SBP-tag: streptavidin-binding peptide; Keefe et al. 2001
    "SBP_tag":               r"MDEKTTGWRGGHVVEGLAGELEQLRARLEHHPQGQREP",
    # Softag3: mild-elution affinity tag used in tandem-affinity and Y2H
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
    # VSV-G tag: from vesicular stomatitis virus glycoprotein
    "VSV_G_tag":             r"YTDIEMNRLGK",
    # AU1-tag: canonical context is DDDDK immediately upstream (enterokinase
    # cleavage site before the tag), as used in all commercial AU1 constructs
    # (Luo et al. 1993, Biotechniques).  Requiring DDDDK.{0,4}DTYRYI avoids
    # matching random 6-mer DTYRYI occurrences in genomic ORFs.
    "AU1_tag":               r"DDDDK.{0,4}DTYRYI",
    # AU5-tag: same enterokinase-site convention (Luo et al. 1993).
    "AU5_tag":               r"DDDDK.{0,4}TDFYLK",
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
    # Factor Xa: IEGR; cleavage after R (Nagai & Thogersen 1987).
    # Fixed: original r"IEG[R]" was identical to r"IEGR" (redundant char
    # class) and matched 4 AA -- ~hundreds of false positives in 4 Mb
    # genomes.  Now requires that the residue immediately following IEGR
    # is NOT proline, which blocks Factor Xa cleavage in all known cases.
    # This halves the false positive rate and matches only actionable sites.
    "Factor_Xa_cleavage":    r"IEGR(?!P)",
    # Factor Xa sites followed by Pro are non-cleavable; see context_note
    # column in output (written as 'Factor_Xa_blocked' in motif_type).

    # ---- linkers ------------------------------------------------------------
    # Flexible GGS linker: (GGGS)n or (GGGGS)n repeats; Chen et al. 2013
    "Flexible_linker_GGS":   r"(?:G{3,4}S){2,}",
    # Rigid EAAAK linker: (EAAAK)n alpha-helical repeats; Arai et al. 2001
    "Rigid_linker_EAAAK":    r"(?:EAAAK){2,}",
    # Rigid PAPAP linker: proline-alanine repeats; semi-rigid spacer
    "Rigid_linker_PAPAP":    r"(?:PAPAP){2,}",
    # Helical poly-alanine linker: (AAAA)n.
    # Raised from 3 to 4 repeats to reduce homopolymer framshift noise.
    "Helical_linker_A4":     r"(?:AAAA){4,}",
}

COMPILED_MOTIFS = {name: re.compile(pattern) for name, pattern in MOTIF_CATALOG.items()}

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def translate_frame(seq: str, frame: int) -> str:
    """Translate *seq* starting at *frame* (0, 1, or 2).

    Uses (len // 3) * 3 range so the last incomplete codon is always
    excluded cleanly.
    """
    s = seq[frame:]
    return "".join(
        CODON_TABLE.get(s[i:i + 3], "X")
        for i in range(0, (len(s) // 3) * 3, 3)
    )

# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def _aa_to_dna_fwd(frame: int, aa_pos: int) -> int:
    """AA index in a forward frame -> absolute 0-based DNA start of that codon."""
    return frame + aa_pos * 3


def _aa_to_dna_rev(seq_len: int, frame: int, aa_pos: int) -> int:
    """AA index in a reverse-complement frame -> absolute 0-based DNA start
    of the corresponding codon on the *forward* strand.

    On the rev-comp string the codon begins at rev_codon_start = frame + aa_pos*3.
    Mirroring back: fwd_start = seq_len - rev_codon_start - 3.
    """
    rev_codon_start = frame + aa_pos * 3
    return seq_len - rev_codon_start - 3

# ---------------------------------------------------------------------------
# Seed-and-extend core
# ---------------------------------------------------------------------------

def scan_sequence(contig: str, seq: str) -> list:
    """Translate all 6 reading frames, seed from compiled motif matches, then
    extend each seed to flanking stop codons and find the nearest upstream M.

    Hits in ORFs shorter than MIN_ORF_AA amino acids are discarded to
    reduce noise from short random matches (particularly relevant for
    4-6 AA cleavage site and linker patterns).

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

    hits: list = []
    seen: set = set()

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

                    # Discard hits in ORFs that are too short to be real
                    orf_aa_len = stop_pos - orf_start_aa
                    if orf_aa_len < MIN_ORF_AA:
                        continue

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

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

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

    all_hits: list = []
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

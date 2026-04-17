#!/usr/bin/env python3
"""Shared utilities for all vectrap scanner modules.

Consolidating rev_comp, open_text, and read_fasta here eliminates the four
identical (but subtly inconsistent) copy-pasted variants that previously
existed across mcs_module, origin_scanner_module, peptide_tag_module, and
regulatory_scanner_module.

rev_comp translation table canonical form: ACGTNacgtn -> TGCANtgcan
"""
import gzip
from pathlib import Path
from typing import Iterator, Tuple

_COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA string.

    Handles A/T/C/G/N in both upper and lower case.  Input is expected to be
    uppercase (as produced by read_fasta) but lowercase is handled defensively.
    """
    return seq.translate(_COMP_TABLE)[::-1]


def open_text(path: str):
    """Open a plain-text or gzip-compressed text file for reading."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) tuples from a FASTA file.

    *name* is the first whitespace-delimited token after '>'.  The sequence is
    returned uppercased with all whitespace stripped.  Both plain and
    gzip-compressed files are supported.
    """
    with open_text(path) as handle:
        name = None
        chunks = []
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

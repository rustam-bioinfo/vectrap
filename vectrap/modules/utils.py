import gzip
import io
from pathlib import Path
from typing import Generator, Tuple


def open_text(path: str | Path, mode: str = "rt", encoding: str = "utf-8"):
    """Open a plain or gzip-compressed text file transparently."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode=mode, encoding=encoding)
    return open(p, mode=mode, encoding=encoding)


def read_fasta(path: str | Path) -> Generator[Tuple[str, str], None, None]:
    """Yield (header, sequence) tuples from a plain or gzipped FASTA file."""
    header = None
    chunks = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def write_fasta(records: list[Tuple[str, str]], path: str | Path, line_width: int = 60) -> None:
    """Write a list of (header, sequence) tuples to a FASTA file."""
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")

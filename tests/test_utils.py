"""Tests for vectrap.modules.utils."""

import gzip
import io
import textwrap
from pathlib import Path

import pytest

from vectrap.modules.utils import open_text, read_fasta, rev_comp, write_fasta


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_plain_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def _write_gz_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# open_text
# ---------------------------------------------------------------------------

class TestOpenText:
    def test_opens_plain_file(self, tmp_path):
        p = tmp_path / "test.txt"
        p.write_text("hello world")
        with open_text(p) as fh:
            assert fh.read() == "hello world"

    def test_opens_gz_file(self, tmp_path):
        p = tmp_path / "test.txt.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("hello gz")
        with open_text(p) as fh:
            assert fh.read() == "hello gz"

    def test_accepts_string_path(self, tmp_path):
        p = tmp_path / "test.txt"
        p.write_text("str path")
        with open_text(str(p)) as fh:
            assert fh.read() == "str path"


# ---------------------------------------------------------------------------
# read_fasta
# ---------------------------------------------------------------------------

class TestReadFasta:
    def test_single_record_plain(self, tmp_path):
        p = tmp_path / "single.fasta"
        _write_plain_fasta(p, [("seq1", "ACGT")])
        records = list(read_fasta(p))
        assert records == [("seq1", "ACGT")]

    def test_multiple_records_plain(self, tmp_path):
        p = tmp_path / "multi.fasta"
        _write_plain_fasta(p, [("a", "AAAA"), ("b", "CCCC"), ("c", "TTTT")])
        records = list(read_fasta(p))
        assert records == [("a", "AAAA"), ("b", "CCCC"), ("c", "TTTT")]

    def test_multiline_sequence(self, tmp_path):
        p = tmp_path / "multiline.fasta"
        p.write_text(">seq1\nACGT\nACGT\nACGT\n")
        records = list(read_fasta(p))
        assert records == [("seq1", "ACGTACGTACGT")]

    def test_gzipped_fasta(self, tmp_path):
        p = tmp_path / "test.fasta.gz"
        _write_gz_fasta(p, [("gz_seq", "GGGGCCCC")])
        records = list(read_fasta(p))
        assert records == [("gz_seq", "GGGGCCCC")]

    def test_uppercase_normalisation(self, tmp_path):
        p = tmp_path / "lower.fasta"
        p.write_text(">seq1\nacgtACGT\n")
        records = list(read_fasta(p))
        assert records == [("seq1", "ACGTACGT")]

    def test_empty_lines_ignored(self, tmp_path):
        p = tmp_path / "gaps.fasta"
        p.write_text(">seq1\n\nACGT\n\n")
        records = list(read_fasta(p))
        assert records == [("seq1", "ACGT")]

    def test_header_with_description(self, tmp_path):
        p = tmp_path / "desc.fasta"
        p.write_text(">seq1 this is a description\nACGT\n")
        records = list(read_fasta(p))
        # read_fasta yields the full header (without '>') including description
        assert records[0][0] == "seq1 this is a description"

    def test_empty_file(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text("")
        records = list(read_fasta(p))
        assert records == []

    def test_sequence_with_n(self, tmp_path):
        p = tmp_path / "n.fasta"
        p.write_text(">seq1\nACGTNNACGT\n")
        records = list(read_fasta(p))
        assert records == [("seq1", "ACGTNNACGT")]


# ---------------------------------------------------------------------------
# rev_comp
# ---------------------------------------------------------------------------

class TestRevComp:
    @pytest.mark.parametrize("seq, expected", [
        ("A", "T"),
        ("T", "A"),
        ("G", "C"),
        ("C", "G"),
        ("ACGT", "ACGT"),       # palindrome
        ("AACCGGTT", "AACCGGTT"),  # palindrome
        ("AAAA", "TTTT"),
        ("ATCG", "CGAT"),
        ("GCATGCAT", "ATGCATGC"),
        ("N", "N"),
        ("ACGTN", "NACGT"),
        ("", ""),
    ])
    def test_known_sequences(self, seq, expected):
        assert rev_comp(seq) == expected

    def test_double_rev_comp_is_identity(self):
        seq = "ATCGATCGATCG"
        assert rev_comp(rev_comp(seq)) == seq

    def test_lowercase_input(self):
        # rev_comp handles lower-case via translation table
        assert rev_comp("acgt") == "acgt"


# ---------------------------------------------------------------------------
# write_fasta
# ---------------------------------------------------------------------------

class TestWriteFasta:
    def test_plain_roundtrip(self, tmp_path):
        p = tmp_path / "out.fasta"
        records = [("seq1", "ACGT"), ("seq2", "GGGGCCCC")]
        write_fasta(records, p)
        assert list(read_fasta(p)) == records

    def test_gz_roundtrip(self, tmp_path):
        p = tmp_path / "out.fasta.gz"
        records = [("seq1", "ACGT"), ("seq2", "TTTTAAAA")]
        write_fasta(records, p)
        assert list(read_fasta(p)) == records

    def test_line_wrapping(self, tmp_path):
        p = tmp_path / "wrapped.fasta"
        long_seq = "A" * 200
        write_fasta([("long", long_seq)], p, line_width=60)
        text = p.read_text()
        lines = text.strip().split("\n")
        # first line is header
        assert lines[0] == ">long"
        # sequence lines should all be <= 60 chars except possibly last
        for line in lines[1:-1]:
            assert len(line) == 60
        # reconstructed sequence should match original
        seq = "".join(lines[1:])
        assert seq == long_seq

    def test_empty_records(self, tmp_path):
        p = tmp_path / "empty.fasta"
        write_fasta([], p)
        assert p.read_text() == ""

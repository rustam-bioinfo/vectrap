"""Tests for vectrap.modules.homology_scanner.

The minimap2-dependent paths are tested with mocked subprocess.run so the
tests pass without minimap2 being installed.  The k-mer scanner and PAF
parsing logic are tested directly against synthetic data.
"""

import pickle
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vectrap.modules.homology_scanner import (
    HomologyHit,
    MIN_LEN,
    _kmer_scan_sequence,
    _parse_paf_line,
    _run_kmer_scanner,
    scan,
)
from vectrap.modules.utils import rev_comp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_paf_line(
    q_name="contig1",
    q_len=1000,
    q_start=100,
    q_end=200,
    strand="+",
    t_name="catalog_entry",
    t_len=100,
    t_start=0,
    t_end=100,
    n_match=98,
    aln_len=100,
    mapq=60,
    tags="de:f:0.02",
) -> str:
    return "\t".join([
        q_name, str(q_len), str(q_start), str(q_end), strand,
        t_name, str(t_len), str(t_start), str(t_end),
        str(n_match), str(aln_len), str(mapq),
        tags,
    ])


def _make_catalog_dir(tmp_path: Path, short_index: dict) -> Path:
    catalog_dir = tmp_path / "catalogs"
    catalog_dir.mkdir()
    (catalog_dir / "combined_long.mmi").touch()
    pkl = catalog_dir / "short_index.pkl"
    with open(pkl, "wb") as fh:
        pickle.dump(short_index, fh)
    return catalog_dir


# ---------------------------------------------------------------------------
# HomologyHit
# ---------------------------------------------------------------------------

class TestHomologyHit:
    def test_length_property(self):
        h = HomologyHit(
            contig="c1", start=10, end=50,
            strand="+", identity=0.95, coverage=0.90,
            catalog_id="cat1",
        )
        assert h.length == 40

    def test_default_source_is_minimap2(self):
        h = HomologyHit(
            contig="c1", start=0, end=10,
            strand="+", identity=1.0, coverage=1.0,
            catalog_id="cat1",
        )
        assert h.source == "minimap2"

    def test_kmer_source(self):
        h = HomologyHit(
            contig="c1", start=0, end=10,
            strand="+", identity=1.0, coverage=1.0,
            catalog_id="cat1", source="kmer",
        )
        assert h.source == "kmer"


# ---------------------------------------------------------------------------
# PAF line parsing
# ---------------------------------------------------------------------------

class TestParsePafLine:
    def test_valid_line_returns_hit(self):
        line = _make_paf_line()
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is not None
        assert hit.contig == "contig1"
        assert hit.start == 100
        assert hit.end == 200
        assert hit.strand == "+"
        assert hit.catalog_id == "catalog_entry"
        assert hit.source == "minimap2"

    def test_identity_from_de_tag(self):
        line = _make_paf_line(tags="de:f:0.05")
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is not None
        assert abs(hit.identity - 0.95) < 1e-9

    def test_identity_from_dv_tag_fallback(self):
        line = _make_paf_line(tags="dv:f:0.05")
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is not None
        assert abs(hit.identity - 0.95) < 1e-9

    def test_identity_from_nm_fallback(self):
        line = _make_paf_line(n_match=95, aln_len=100, tags="")
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is not None
        assert abs(hit.identity - 0.95) < 1e-9

    def test_filtered_by_low_identity(self):
        line = _make_paf_line(tags="de:f:0.15")  # identity = 0.85
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is None

    def test_filtered_by_low_coverage(self):
        line = _make_paf_line(t_start=0, t_end=50, t_len=100, tags="de:f:0.02")
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is None

    def test_reverse_strand_hit(self):
        line = _make_paf_line(strand="-")
        hit = _parse_paf_line(line, min_identity=0.80, min_coverage=0.80)
        assert hit is not None
        assert hit.strand == "-"

    def test_empty_line_returns_none(self):
        assert _parse_paf_line("", 0.90, 0.80) is None

    def test_comment_line_returns_none(self):
        assert _parse_paf_line("# comment", 0.90, 0.80) is None

    def test_short_line_returns_none(self):
        assert _parse_paf_line("too\tshort", 0.90, 0.80) is None

    def test_coverage_calculation(self):
        line = _make_paf_line(t_start=0, t_end=90, t_len=100, tags="de:f:0.01")
        hit = _parse_paf_line(line, min_identity=0.90, min_coverage=0.80)
        assert hit is not None
        assert abs(hit.coverage - 0.90) < 1e-9

    def test_contig_name_stripped_at_whitespace(self):
        line = _make_paf_line(q_name="contig1 extra description")
        hit = _parse_paf_line(line, min_identity=0.80, min_coverage=0.80)
        assert hit is not None
        assert hit.contig == "contig1"


# ---------------------------------------------------------------------------
# k-mer scanner
# ---------------------------------------------------------------------------

class TestKmerScanSequence:
    def test_exact_forward_match(self):
        kmer = "ATCGATCG"
        contig = "NNNN" + kmer + "NNNN"
        hits = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        fwd = [h for h in hits if h.strand == "+"]
        assert len(fwd) == 1
        assert fwd[0].start == 4
        assert fwd[0].end == 4 + len(kmer)
        assert fwd[0].identity == 1.0
        assert fwd[0].coverage == 1.0
        assert fwd[0].source == "kmer"

    def test_exact_reverse_match(self):
        """rev_comp(kmer) placed in contig must produce a '-' hit at those coords."""
        kmer = "ATCGATCG"
        rc   = rev_comp(kmer)       # "CGATCGAT"
        # contig: NNNN + CGATCGAT + NNNN  (length=16)
        contig  = "NNNN" + rc + "NNNN"
        seq_len = len(contig)       # 16
        klen    = len(kmer)         # 8

        hits = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        rev_hits = [h for h in hits if h.strand == "-"]
        assert len(rev_hits) == 1

        # rc_seq = rev_comp(contig) = NNNN + ATCGATCG + NNNN
        # kmer found at rc_idx=4
        # fwd_start = 16 - (4 + 8) = 4
        # fwd_end   = 16 - 4       = 12
        assert rev_hits[0].start == 4
        assert rev_hits[0].end   == 12

    def test_multiple_occurrences_found(self):
        kmer = "ATCG"
        contig = kmer + "NNNN" + kmer
        hits = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        fwd = [h for h in hits if h.strand == "+"]
        assert len(fwd) == 2
        assert fwd[0].start == 0
        assert fwd[1].start == 8

    def test_palindrome_not_double_counted(self):
        palindrome = "ACGT"   # rev_comp("ACGT") == "ACGT"
        assert rev_comp(palindrome) == palindrome
        contig = palindrome * 3
        hits = _kmer_scan_sequence("c1", contig, {palindrome: ["cat1"]})
        strands = {h.strand for h in hits}
        assert "-" not in strands

    def test_kmer_longer_than_contig_skipped(self):
        kmer = "A" * 100
        hits = _kmer_scan_sequence("c1", "AAAA", {kmer: ["cat1"]})
        assert hits == []

    def test_multiple_catalog_ids_for_same_kmer(self):
        kmer = "ATCGATCG"
        hits = _kmer_scan_sequence("c1", kmer, {kmer: ["cat1", "cat2", "cat3"]})
        fwd = [h for h in hits if h.strand == "+"]
        assert len(fwd) == 3
        assert {h.catalog_id for h in fwd} == {"cat1", "cat2", "cat3"}

    def test_empty_index_returns_no_hits(self):
        assert _kmer_scan_sequence("c1", "ACGTACGT", {}) == []

    def test_no_match_returns_empty(self):
        assert _kmer_scan_sequence("c1", "AAAAAAA", {"CCCCCC": ["cat1"]}) == []

    def test_rc_hit_coordinates_roundtrip(self):
        """contig[fwd_start:fwd_end] must equal rev_comp(kmer)."""
        kmer   = "GGTTAACC"
        rc     = rev_comp(kmer)
        contig = "TTTT" + rc + "TTTT"
        hits   = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        rev_hits = [h for h in hits if h.strand == "-"]
        assert len(rev_hits) == 1
        h = rev_hits[0]
        assert contig[h.start:h.end] == rc


# ---------------------------------------------------------------------------
# _run_kmer_scanner (file-level)
# ---------------------------------------------------------------------------

class TestRunKmerScanner:
    def test_scans_all_contigs(self, tmp_path):
        kmer  = "ATCGATCG"
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\n" + kmer + "NNNN\n>c2\nNNNN" + kmer + "\n")
        pkl = tmp_path / "short_index.pkl"
        with open(pkl, "wb") as fh:
            pickle.dump({kmer: ["cat1"]}, fh)
        hits = _run_kmer_scanner(fasta, pkl)
        contigs = {h.contig for h in hits if h.strand == "+"}
        assert "c1" in contigs
        assert "c2" in contigs

    def test_empty_index_returns_empty(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        pkl = tmp_path / "short_index.pkl"
        with open(pkl, "wb") as fh:
            pickle.dump({}, fh)
        assert _run_kmer_scanner(fasta, pkl) == []


# ---------------------------------------------------------------------------
# scan() -- integration (minimap2 mocked)
# ---------------------------------------------------------------------------

class TestScan:
    def test_raises_if_mmi_missing(self, tmp_path):
        catalog_dir = tmp_path / "catalogs"
        catalog_dir.mkdir()
        (catalog_dir / "short_index.pkl").touch()
        with pytest.raises(FileNotFoundError, match="minimap2 index not found"):
            scan(tmp_path / "q.fasta", catalog_dir)

    def test_raises_if_pkl_missing(self, tmp_path):
        catalog_dir = tmp_path / "catalogs"
        catalog_dir.mkdir()
        (catalog_dir / "combined_long.mmi").touch()
        with pytest.raises(FileNotFoundError, match="k-mer index not found"):
            scan(tmp_path / "q.fasta", catalog_dir)

    def test_returns_combined_hits(self, tmp_path):
        kmer  = "ATCGATCG"
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\n" + kmer + "NNNNNNNNNN\n")
        catalog_dir = _make_catalog_dir(tmp_path, {kmer: ["short_cat"]})

        fake_paf = _make_paf_line(
            q_name="c1", q_start=0, q_end=50,
            t_name="long_cat", t_len=50, t_end=50,
            tags="de:f:0.01",
        )
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        mock_result.stdout = fake_paf + "\n"

        with patch("subprocess.run", return_value=mock_result):
            hits = scan(fasta, catalog_dir)

        sources = {h.source for h in hits}
        assert "minimap2" in sources
        assert "kmer" in sources

    def test_minimap2_error_raises_runtime_error(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\nACGT\n")
        catalog_dir = _make_catalog_dir(tmp_path, {})

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "fatal error"
        mock_result.stdout = ""

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(RuntimeError, match="minimap2 exited"):
                scan(fasta, catalog_dir)

    def test_minimap2_not_found_raises_file_not_found(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\nACGT\n")
        catalog_dir = _make_catalog_dir(tmp_path, {})

        with patch("subprocess.run", side_effect=FileNotFoundError):
            with pytest.raises(FileNotFoundError, match="minimap2 not found"):
                scan(fasta, catalog_dir)

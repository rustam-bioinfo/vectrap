"""Tests for vectrap.modules.homology_scanner.

The mappy-dependent paths are tested with a mocked mappy.Aligner so the
tests pass without a real .mmi index.  The k-mer scanner and identity/coverage
calculation logic are tested directly against synthetic data.
"""

import pickle
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from vectrap.modules.homology_scanner import (
    HomologyHit,
    _identity_from_hit,
    _kmer_scan_sequence,
    _run_kmer_scanner,
    _run_mappy,
    scan,
)
from vectrap.modules.utils import rev_comp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mappy_hit(
    ctg="catalog_entry",
    ctg_len=100,
    r_st=0,
    r_en=100,
    q_st=100,
    q_en=200,
    strand=1,
    mlen=98,
    blen=100,
) -> SimpleNamespace:
    """Build a fake mappy.Alignment-like object."""
    return SimpleNamespace(
        ctg=ctg, ctg_len=ctg_len,
        r_st=r_st, r_en=r_en,
        q_st=q_st, q_en=q_en,
        strand=strand,
        mlen=mlen, blen=blen,
    )


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

    def test_default_source_is_mappy(self):
        h = HomologyHit(
            contig="c1", start=0, end=10,
            strand="+", identity=1.0, coverage=1.0,
            catalog_id="cat1",
        )
        assert h.source == "mappy"

    def test_kmer_source(self):
        h = HomologyHit(
            contig="c1", start=0, end=10,
            strand="+", identity=1.0, coverage=1.0,
            catalog_id="cat1", source="kmer",
        )
        assert h.source == "kmer"


# ---------------------------------------------------------------------------
# _identity_from_hit
# ---------------------------------------------------------------------------

class TestIdentityFromHit:
    def test_normal(self):
        hit = _make_mappy_hit(mlen=95, blen=100)
        assert abs(_identity_from_hit(hit) - 0.95) < 1e-9

    def test_perfect(self):
        hit = _make_mappy_hit(mlen=100, blen=100)
        assert _identity_from_hit(hit) == 1.0

    def test_zero_blen(self):
        hit = _make_mappy_hit(mlen=0, blen=0)
        assert _identity_from_hit(hit) == 0.0


# ---------------------------------------------------------------------------
# _run_mappy
# ---------------------------------------------------------------------------

class TestRunMappy:
    def _make_aligner_mock(self, hits):
        aligner = MagicMock()
        aligner.__bool__ = lambda s: True
        aligner.map.return_value = iter(hits)
        return aligner

    def test_returns_hits_above_thresholds(self, tmp_path):
        fasta = tmp_path / "q.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        mmi = tmp_path / "idx.mmi"
        mmi.touch()

        fake_hit = _make_mappy_hit(mlen=98, blen=100, r_st=0, r_en=100, ctg_len=100)
        aligner_mock = self._make_aligner_mock([fake_hit])

        with patch("mappy.Aligner", return_value=aligner_mock):
            hits = _run_mappy(fasta, mmi, min_identity=0.90, min_coverage=0.80)

        assert len(hits) == 1
        assert hits[0].identity == pytest.approx(0.98)
        assert hits[0].coverage == pytest.approx(1.0)
        assert hits[0].source == "mappy"

    def test_filters_low_identity(self, tmp_path):
        fasta = tmp_path / "q.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        mmi = tmp_path / "idx.mmi"
        mmi.touch()

        fake_hit = _make_mappy_hit(mlen=80, blen=100)  # identity=0.80
        aligner_mock = self._make_aligner_mock([fake_hit])

        with patch("mappy.Aligner", return_value=aligner_mock):
            hits = _run_mappy(fasta, mmi, min_identity=0.90, min_coverage=0.80)

        assert hits == []

    def test_filters_low_coverage(self, tmp_path):
        fasta = tmp_path / "q.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        mmi = tmp_path / "idx.mmi"
        mmi.touch()

        # coverage = (r_en - r_st) / ctg_len = 50/100 = 0.50
        fake_hit = _make_mappy_hit(mlen=98, blen=100, r_st=0, r_en=50, ctg_len=100)
        aligner_mock = self._make_aligner_mock([fake_hit])

        with patch("mappy.Aligner", return_value=aligner_mock):
            hits = _run_mappy(fasta, mmi, min_identity=0.90, min_coverage=0.80)

        assert hits == []

    def test_reverse_strand(self, tmp_path):
        fasta = tmp_path / "q.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        mmi = tmp_path / "idx.mmi"
        mmi.touch()

        fake_hit = _make_mappy_hit(mlen=98, blen=100, r_st=0, r_en=100, ctg_len=100, strand=-1)
        aligner_mock = self._make_aligner_mock([fake_hit])

        with patch("mappy.Aligner", return_value=aligner_mock):
            hits = _run_mappy(fasta, mmi, min_identity=0.90, min_coverage=0.80)

        assert hits[0].strand == "-"

    def test_raises_if_aligner_fails(self, tmp_path):
        fasta = tmp_path / "q.fasta"
        fasta.write_text(">c1\nACGT\n")
        mmi = tmp_path / "idx.mmi"
        mmi.touch()

        bad_aligner = MagicMock()
        bad_aligner.__bool__ = lambda s: False

        with patch("mappy.Aligner", return_value=bad_aligner):
            with pytest.raises(RuntimeError, match="mappy could not load index"):
                _run_mappy(fasta, mmi, min_identity=0.90, min_coverage=0.80)


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
        kmer = "ATCGATCG"
        rc   = rev_comp(kmer)
        contig  = "NNNN" + rc + "NNNN"  # length=16
        seq_len = len(contig)
        klen    = len(kmer)

        hits = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        rev_hits = [h for h in hits if h.strand == "-"]
        assert len(rev_hits) == 1
        assert rev_hits[0].start == 4
        assert rev_hits[0].end   == 12

    def test_multiple_occurrences_found(self):
        kmer = "ATCG"
        contig = kmer + "NNNN" + kmer
        hits = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        fwd = [h for h in hits if h.strand == "+"]
        assert len(fwd) == 2

    def test_palindrome_not_double_counted(self):
        palindrome = "ACGT"
        assert rev_comp(palindrome) == palindrome
        hits = _kmer_scan_sequence("c1", palindrome * 3, {palindrome: ["cat1"]})
        assert "-" not in {h.strand for h in hits}

    def test_kmer_longer_than_contig_skipped(self):
        assert _kmer_scan_sequence("c1", "AAAA", {"A" * 100: ["cat1"]}) == []

    def test_multiple_catalog_ids(self):
        kmer = "ATCGATCG"
        hits = _kmer_scan_sequence("c1", kmer, {kmer: ["cat1", "cat2", "cat3"]})
        assert {h.catalog_id for h in hits if h.strand == "+"} == {"cat1", "cat2", "cat3"}

    def test_empty_index(self):
        assert _kmer_scan_sequence("c1", "ACGTACGT", {}) == []

    def test_rc_coordinates_roundtrip(self):
        kmer   = "GGTTAACC"
        rc     = rev_comp(kmer)
        contig = "TTTT" + rc + "TTTT"
        hits   = _kmer_scan_sequence("c1", contig, {kmer: ["cat1"]})
        rev_hits = [h for h in hits if h.strand == "-"]
        assert len(rev_hits) == 1
        assert contig[rev_hits[0].start:rev_hits[0].end] == rc


# ---------------------------------------------------------------------------
# _run_kmer_scanner
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
        assert {h.contig for h in hits if h.strand == "+"} == {"c1", "c2"}

    def test_empty_index(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\nACGTACGT\n")
        pkl = tmp_path / "short_index.pkl"
        with open(pkl, "wb") as fh:
            pickle.dump({}, fh)
        assert _run_kmer_scanner(fasta, pkl) == []


# ---------------------------------------------------------------------------
# scan() -- integration
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

        fake_hit = _make_mappy_hit(
            ctg="long_cat", ctg_len=50,
            r_st=0, r_en=50, q_st=0, q_en=50,
            mlen=49, blen=50,
        )
        aligner_mock = MagicMock()
        aligner_mock.__bool__ = lambda s: True
        aligner_mock.map.return_value = iter([fake_hit])

        with patch("mappy.Aligner", return_value=aligner_mock):
            hits = scan(fasta, catalog_dir)

        sources = {h.source for h in hits}
        assert "mappy" in sources
        assert "kmer" in sources

    def test_mappy_index_load_failure_raises(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">c1\nACGT\n")
        catalog_dir = _make_catalog_dir(tmp_path, {})

        bad_aligner = MagicMock()
        bad_aligner.__bool__ = lambda s: False

        with patch("mappy.Aligner", return_value=bad_aligner):
            with pytest.raises(RuntimeError, match="mappy could not load index"):
                scan(fasta, catalog_dir)

"""Tests for vectrap.cli.build_db utility functions.

These tests exercise the logic that can be tested without network access or
minimap2: download_catalogs (mocked), verify_manifest, and build_indexes
(short sequences only, skipping the minimap2 call).
"""

import gzip
import hashlib
import pickle
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vectrap.cli.build_db import (
    MIN_LEN,
    CATALOG_FILES,
    build_indexes,
    download_catalogs,
    verify_manifest,
)
from vectrap.modules.utils import write_fasta


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gz_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a gzipped FASTA to *path*."""
    with gzip.open(path, "wt") as fh:
        for hdr, seq in records:
            fh.write(f">{hdr}\n{seq}\n")


def _write_manifest(catalog_dir: Path, entries: list[tuple[str, int, float, str]]) -> None:
    """Write a catalog_manifest.tsv with the given entries.

    Each entry is (filename, n_seqs, size_mb, md5).
    """
    manifest = catalog_dir / "catalog_manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("filename\tn_seqs\tsize_mb\tmd5\n")
        for fname, n, size, md5 in entries:
            fh.write(f"{fname}\t{n}\t{size}\t{md5}\n")


# ---------------------------------------------------------------------------
# download_catalogs
# ---------------------------------------------------------------------------

class TestDownloadCatalogs:
    def test_creates_directory(self, tmp_path):
        catalog_dir = tmp_path / "catalogs" / "nested"
        with patch("urllib.request.urlretrieve") as mock_dl:
            mock_dl.side_effect = lambda url, dest: Path(dest).touch()
            download_catalogs(catalog_dir)
        assert catalog_dir.exists()

    def test_downloads_all_catalog_files(self, tmp_path):
        catalog_dir = tmp_path / "catalogs"
        downloaded = []
        def fake_dl(url, dest):
            Path(dest).touch()
            downloaded.append(Path(dest).name)
        with patch("urllib.request.urlretrieve", side_effect=fake_dl):
            download_catalogs(catalog_dir)
        for fname in CATALOG_FILES:
            assert fname in downloaded

    def test_skips_existing_files(self, tmp_path):
        catalog_dir = tmp_path / "catalogs"
        catalog_dir.mkdir()
        # Pre-create the first catalog file
        first = catalog_dir / CATALOG_FILES[0]
        first.touch()
        downloaded = []
        def fake_dl(url, dest):
            Path(dest).touch()
            downloaded.append(Path(dest).name)
        with patch("urllib.request.urlretrieve", side_effect=fake_dl):
            download_catalogs(catalog_dir)
        assert CATALOG_FILES[0] not in downloaded
        # All other files should have been downloaded
        for fname in CATALOG_FILES[1:]:
            assert fname in downloaded


# ---------------------------------------------------------------------------
# verify_manifest
# ---------------------------------------------------------------------------

class TestVerifyManifest:
    def test_passes_with_correct_checksums(self, tmp_path):
        fname = "catalog_CDS.fasta.gz"
        content = b"test content"
        (tmp_path / fname).write_bytes(content)
        md5 = hashlib.md5(content).hexdigest()
        _write_manifest(tmp_path, [(fname, 1, 0.001, md5)])
        # Should not raise or call sys.exit
        verify_manifest(tmp_path)

    def test_fails_on_checksum_mismatch(self, tmp_path):
        fname = "catalog_CDS.fasta.gz"
        (tmp_path / fname).write_bytes(b"real content")
        _write_manifest(tmp_path, [(fname, 1, 0.001, "deadbeef" * 4)])
        with pytest.raises(SystemExit):
            verify_manifest(tmp_path)

    def test_fails_on_missing_file(self, tmp_path):
        _write_manifest(tmp_path, [("catalog_CDS.fasta.gz", 1, 0.001, "abc123")])
        # File does not exist
        with pytest.raises(SystemExit):
            verify_manifest(tmp_path)

    def test_no_manifest_warns_but_does_not_exit(self, tmp_path, capsys):
        """If manifest is absent, a warning is printed but no SystemExit."""
        verify_manifest(tmp_path)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out

    def test_ignores_short_lines(self, tmp_path):
        manifest = tmp_path / "catalog_manifest.tsv"
        manifest.write_text("filename\tn_seqs\tsize_mb\tmd5\nbad_line\n")
        # Should not raise -- short lines are skipped
        verify_manifest(tmp_path)


# ---------------------------------------------------------------------------
# build_indexes
# ---------------------------------------------------------------------------

class TestBuildIndexes:
    def _make_catalog_dir(self, tmp_path: Path) -> Path:
        """Create a minimal catalog directory with one short and one long sequence."""
        catalog_dir = tmp_path / "catalogs"
        catalog_dir.mkdir()
        short_seq = "A" * (MIN_LEN - 1)   # goes to k-mer hash
        long_seq  = "G" * (MIN_LEN + 10)  # goes to minimap2 FASTA
        _gz_fasta(
            catalog_dir / "catalog_CDS.fasta.gz",
            [("short_seq", short_seq), ("long_seq", long_seq)],
        )
        return catalog_dir

    def test_short_index_pkl_created(self, tmp_path):
        catalog_dir = self._make_catalog_dir(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            build_indexes(catalog_dir)
        pkl = catalog_dir / "short_index.pkl"
        assert pkl.exists()

    def test_short_index_contains_correct_sequence(self, tmp_path):
        catalog_dir = self._make_catalog_dir(tmp_path)
        short_seq = "A" * (MIN_LEN - 1)
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            build_indexes(catalog_dir)
        pkl = catalog_dir / "short_index.pkl"
        with open(pkl, "rb") as fh:
            index = pickle.load(fh)
        assert short_seq in index
        assert "short_seq" in index[short_seq]

    def test_long_fasta_written(self, tmp_path):
        catalog_dir = self._make_catalog_dir(tmp_path)
        long_seq = "G" * (MIN_LEN + 10)
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            build_indexes(catalog_dir)
        from vectrap.modules.utils import read_fasta
        records = list(read_fasta(catalog_dir / "combined_long.fasta"))
        seqs = [s for _, s in records]
        assert long_seq in seqs

    def test_long_seq_not_in_short_index(self, tmp_path):
        catalog_dir = self._make_catalog_dir(tmp_path)
        long_seq = "G" * (MIN_LEN + 10)
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            build_indexes(catalog_dir)
        with open(catalog_dir / "short_index.pkl", "rb") as fh:
            index = pickle.load(fh)
        assert long_seq not in index

    def test_minimap2_failure_raises_sys_exit(self, tmp_path):
        catalog_dir = self._make_catalog_dir(tmp_path)
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "minimap2 error"
        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(SystemExit):
                build_indexes(catalog_dir)

    def test_missing_catalog_file_skipped(self, tmp_path):
        """Missing catalog files produce a warning but do not crash."""
        catalog_dir = tmp_path / "catalogs"
        catalog_dir.mkdir()
        # Don't create any catalog files
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        with patch("subprocess.run", return_value=mock_result):
            # Should complete without raising even with no catalog files
            build_indexes(catalog_dir)

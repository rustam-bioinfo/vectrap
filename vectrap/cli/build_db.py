"""Prepare VecTrap catalog indexes.

Installed as the `vectrap-build-db` command via pyproject.toml entry point.

Usage
-----
    # Download catalog files from Zenodo and build indexes
    vectrap-build-db --download

    # Use locally provided catalog files
    vectrap-build-db --catalog-dir /path/to/catalogs/

Outputs (written to the catalog directory)
------------------------------------------
    combined_long.fasta   sequences >= MIN_LEN bp
    combined_long.mmi     mappy/minimap2 index
    short_index.pkl       exact k-mer hash {seq: [id, ...]}
"""

import argparse
import hashlib
import pickle
import sys
import urllib.request
from pathlib import Path

import mappy

from vectrap.modules.utils import read_fasta, write_fasta

ZENODO_DOI  = "10.5281/zenodo.20844271"
ZENODO_BASE = "https://zenodo.org/records/20844271/files"

CATALOG_FILES = [
    "catalog_CDS.fasta.gz",
    "catalog_LTR.fasta.gz",
    "catalog_RBS.fasta.gz",
    "catalog_enhancer.fasta.gz",
    "catalog_misc_feature.fasta.gz",
    "catalog_misc_recomb.fasta.gz",
    "catalog_mobile_element.fasta.gz",
    "catalog_oriT.fasta.gz",
    "catalog_primer_bind.fasta.gz",
    "catalog_promoter.fasta.gz",
    "catalog_protein_bind.fasta.gz",
    "catalog_regulatory.fasta.gz",
    "catalog_rep_origin.fasta.gz",
    "catalog_terminator.fasta.gz",
    "catalog_manifest.tsv",
]

MIN_LEN = 50  # bp threshold between mappy aligner and k-mer hash


def _default_catalog_dir() -> Path:
    """Return the default catalog directory inside the installed package."""
    return Path(__file__).resolve().parents[2] / "vectrap" / "catalogs"


def download_catalogs(catalog_dir: Path) -> None:
    catalog_dir.mkdir(parents=True, exist_ok=True)
    for fname in CATALOG_FILES:
        dest = catalog_dir / fname
        if dest.exists():
            print(f"  already present, skipping: {fname}")
            continue
        url = f"{ZENODO_BASE}/{fname}?download=1"
        print(f"  downloading {fname} ...")
        urllib.request.urlretrieve(url, dest)
    print("Download complete.")


def verify_manifest(catalog_dir: Path) -> None:
    manifest_path = catalog_dir / "catalog_manifest.tsv"
    if not manifest_path.exists():
        print("WARNING: catalog_manifest.tsv not found, skipping checksum verification.")
        return
    errors = []
    with open(manifest_path) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            fname, _n, _s, expected_md5 = parts[0], parts[1], parts[2], parts[3]
            fpath = catalog_dir / fname
            if not fpath.exists():
                errors.append(f"  MISSING: {fname}")
                continue
            actual_md5 = hashlib.md5(fpath.read_bytes()).hexdigest()
            if actual_md5 != expected_md5:
                errors.append(f"  CHECKSUM MISMATCH: {fname}")
    if errors:
        print("Catalog verification FAILED:")
        for e in errors:
            print(e)
        sys.exit(1)
    print("Catalog verification passed.")


def build_indexes(catalog_dir: Path) -> None:
    long_records = []
    short_index: dict = {}

    for fname in CATALOG_FILES:
        if not fname.endswith(".fasta.gz"):
            continue
        fpath = catalog_dir / fname
        if not fpath.exists():
            print(f"  WARNING: {fname} not found, skipping.")
            continue
        for seq_id, seq in read_fasta(fpath):
            if len(seq) >= MIN_LEN:
                long_records.append((seq_id, seq))
            else:
                short_index.setdefault(seq, []).append(seq_id)

    print(f"  long sequences (>= {MIN_LEN} bp): {len(long_records):,}")
    print(f"  short sequences (<  {MIN_LEN} bp): {len(short_index):,} unique")

    # Write combined long FASTA
    long_fasta = catalog_dir / "combined_long.fasta"
    write_fasta(long_records, long_fasta)
    print(f"  wrote {long_fasta}")

    # Build mappy index (.mmi) -- no external binary required
    mmi_path = catalog_dir / "combined_long.mmi"
    print("  building mappy index ...")
    aligner = mappy.Aligner(str(long_fasta), preset="map-ont", best_n=10)
    if not aligner:
        print("ERROR: mappy could not index the long catalog FASTA.", file=sys.stderr)
        sys.exit(1)
    aligner.dump_index(str(mmi_path))
    print(f"  wrote {mmi_path}")

    # Write k-mer hash for short sequences
    pkl_path = catalog_dir / "short_index.pkl"
    with open(pkl_path, "wb") as fh:
        pickle.dump(short_index, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  wrote {pkl_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vectrap-build-db",
        description="Download and index the VecTrap sequence catalog.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--download",
        action="store_true",
        help="Download catalog files from Zenodo (DOI: 10.5281/zenodo.20844271) and build indexes.",
    )
    group.add_argument(
        "--catalog-dir",
        metavar="DIR",
        help="Path to a directory containing locally provided catalog FASTA files.",
    )
    args = parser.parse_args()

    if args.download:
        catalog_dir = _default_catalog_dir()
        print(f"Downloading catalog files to {catalog_dir} ...")
        download_catalogs(catalog_dir)
    else:
        catalog_dir = Path(args.catalog_dir).resolve()

    print("Verifying catalog files ...")
    verify_manifest(catalog_dir)

    print("Building indexes ...")
    build_indexes(catalog_dir)
    print("Done. Indexes are ready.")


if __name__ == "__main__":
    main()

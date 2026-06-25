"""Prepare VecTrap catalog indexes.

Usage
-----
    # Download catalog files from Zenodo and build indexes
    python db/build_db.py --download

    # Use locally provided catalog files
    python db/build_db.py --catalog-dir /path/to/catalogs/

Outputs
-------
    <catalog_dir>/combined_long.fasta   sequences >= BLAST_MIN_LEN bp
    <catalog_dir>/combined_long.mmi     minimap2 index
    <catalog_dir>/short_index.pkl       exact k-mer hash {seq: [id, ...]}
"""

import argparse
import hashlib
import pickle
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vectrap.modules.utils import read_fasta, write_fasta

ZENODO_DOI = "10.5281/zenodo.20844271"
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

MIN_LEN = 50  # bp threshold between minimap2 and k-mer hash


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
            fname, n_seqs, size_mb, expected_md5 = parts[0], parts[1], parts[2], parts[3]
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
    short_index = {}

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
    print(f"  short sequences (< {MIN_LEN} bp):  {len(short_index):,} unique")

    long_fasta = catalog_dir / "combined_long.fasta"
    write_fasta(long_records, long_fasta)
    print(f"  wrote {long_fasta}")

    mmi_path = catalog_dir / "combined_long.mmi"
    cmd = ["minimap2", "-d", str(mmi_path), str(long_fasta)]
    print(f"  building minimap2 index: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR: minimap2 failed:")
        print(result.stderr)
        sys.exit(1)
    print(f"  wrote {mmi_path}")

    pkl_path = catalog_dir / "short_index.pkl"
    with open(pkl_path, "wb") as fh:
        pickle.dump(short_index, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  wrote {pkl_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare VecTrap catalog indexes.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--download",
        action="store_true",
        help="Download catalog files from Zenodo and build indexes.",
    )
    group.add_argument(
        "--catalog-dir",
        metavar="DIR",
        help="Path to directory containing locally provided catalog files.",
    )
    args = parser.parse_args()

    if args.download:
        catalog_dir = Path(__file__).resolve().parent.parent / "vectrap" / "catalogs"
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

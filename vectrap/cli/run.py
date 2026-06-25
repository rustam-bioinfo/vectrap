"""Main VecTrap pipeline entry point.

Installed as the `vectrap` command via pyproject.toml entry point.

Usage
-----
    vectrap -i assembly.fasta -o results/ -c /path/to/catalogs/
"""

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vectrap",
        description="Detect and classify synthetic vector contamination in genome assemblies.",
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="FASTA",
        help="Input assembly FASTA file (plain or gzipped).",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )
    parser.add_argument(
        "-c", "--catalog-dir",
        required=True,
        metavar="DIR",
        help="Path to the catalog directory containing indexes built by vectrap-build-db.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.90,
        metavar="FLOAT",
        help="Minimum sequence identity for minimap2 hits (default: 0.90).",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.80,
        metavar="FLOAT",
        help="Minimum fraction of catalog sequence covered by a hit (default: 0.80).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )
    args = parser.parse_args()

    # pipeline logic will be wired here once scanner and scorer modules are complete
    raise NotImplementedError(
        "Pipeline modules are not yet implemented. "
        "See https://github.com/rustam-bioinfo/vectrap for development status."
    )


if __name__ == "__main__":
    main()

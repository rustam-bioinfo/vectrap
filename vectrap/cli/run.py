"""Main VecTrap pipeline entry point.

Installed as the `vectrap` command via pyproject.toml entry point.

Usage
-----
    vectrap -i assembly.fasta -o results/ -c /path/to/catalogs/
"""

import argparse
import csv
import sys
from pathlib import Path

from vectrap.modules.homology_scanner import scan
from vectrap.modules.scorer import summarize


# Default catalog directory bundled inside the package (populated by
# vectrap-build-db --download).
_DEFAULT_CATALOG_DIR = Path(__file__).resolve().parents[2] / "vectrap" / "catalogs"


def _write_hits_tsv(hits, output_dir: Path) -> Path:
    """Write all hits to vectrap_hits.tsv and return the path."""
    out_path = output_dir / "vectrap_hits.tsv"
    fieldnames = [
        "contig", "start", "end", "strand", "length",
        "identity", "coverage", "catalog_id", "source",
        "feature_type", "label", "tier", "confidence", "reasoning",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for h in hits:
            writer.writerow({
                "contig":       h.contig,
                "start":        h.start,
                "end":          h.end,
                "strand":       h.strand,
                "length":       h.length,
                "identity":     f"{h.identity:.4f}",
                "coverage":     f"{h.coverage:.4f}",
                "catalog_id":   h.catalog_id,
                "source":       h.source,
                "feature_type": h.feature_type,
                "label":        h.label,
                "tier":         h.tier,
                "confidence":   h.confidence,
                "reasoning":    h.reasoning,
            })
    return out_path


def _write_verdicts_tsv(summaries, output_dir: Path) -> Path:
    """Write per-contig hit summaries to vectrap_verdicts.tsv and return the path."""
    out_path = output_dir / "vectrap_verdicts.tsv"
    fieldnames = [
        "contig",
        "total_hits",
        "engineered_hits",
        "context_dependent_hits",
        "weak_hits",
        "unannotated_hits",
        "unique_feature_types",
        "unique_labels",
        "covered_bp",
        "top_label",
        "top_tier",
        "top_confidence",
        "evidence_summary",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow({
                "contig":                   s.contig,
                "total_hits":               s.total_hits,
                "engineered_hits":          s.engineered_hits,
                "context_dependent_hits":   s.context_dependent_hits,
                "weak_hits":                s.weak_hits,
                "unannotated_hits":         s.unannotated_hits,
                "unique_feature_types":     s.unique_feature_types,
                "unique_labels":            s.unique_labels,
                "covered_bp":               s.covered_bp,
                "top_label":                s.top_label,
                "top_tier":                 s.top_tier,
                "top_confidence":           s.top_confidence,
                "evidence_summary":         s.evidence_summary,
            })
    return out_path


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
        metavar="DIR",
        default=str(_DEFAULT_CATALOG_DIR),
        help=(
            "Path to the catalog directory containing indexes built by "
            "vectrap-build-db. Defaults to the bundled catalog directory "
            f"({_DEFAULT_CATALOG_DIR})."
        ),
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.90,
        metavar="FLOAT",
        help="Minimum sequence identity for mappy hits (default: 0.90).",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.80,
        metavar="FLOAT",
        help="Minimum fraction of catalog sequence covered by a hit (default: 0.80).",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        metavar="INT",
        help="Number of threads for mappy aligner (default: 4).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress messages to stderr.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )
    args = parser.parse_args()

    input_path  = Path(args.input).resolve()
    output_dir  = Path(args.output).resolve()
    catalog_dir = Path(args.catalog_dir).resolve()

    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    hits = scan(
        query_fasta=input_path,
        catalog_dir=catalog_dir,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        threads=args.threads,
        verbose=args.verbose,
    )

    hits_path = _write_hits_tsv(hits, output_dir)
    print(f"  hits    : {len(hits):,} written to {hits_path}")

    summaries = summarize(hits)
    verdicts_path = _write_verdicts_tsv(summaries, output_dir)
    print(f"  contigs : {len(summaries):,} with hits written to {verdicts_path}")
    print("Done.")


if __name__ == "__main__":
    main()

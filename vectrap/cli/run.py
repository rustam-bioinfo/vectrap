"""Main VecTrap pipeline entry point.

Installed as the `vectrap` command via pyproject.toml entry point.

Usage
-----
    # Single file
    vectrap -i assembly.fasta -o results/

    # Directory of FASTA files
    vectrap -i genomes/ -o results/

Output files (prefixed with sample name)
-----------------------------------------
    <sample>_hits.tsv
    <sample>_contig_verdicts.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

from vectrap.modules.homology_scanner import scan
from vectrap.modules.scorer import summarize
from vectrap.modules.utils import read_fasta


_DEFAULT_CATALOG_DIR = Path(__file__).resolve().parents[2] / "vectrap" / "catalogs"

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")


def _collect_inputs(path: Path) -> list[Path]:
    """Return a list of FASTA files to process.

    If *path* is a file, return ``[path]``.
    If *path* is a directory, return all FASTA files found in it (non-recursive).
    """
    if not path.exists():
        print(f"ERROR: input path not found: {path}", file=sys.stderr)
        sys.exit(1)

    if path.is_file():
        return [path]

    if path.is_dir():
        files = [
            f for f in sorted(path.iterdir())
            if f.is_file() and any(str(f).endswith(s) for s in _FASTA_SUFFIXES)
        ]
        if not files:
            print(
                f"ERROR: no FASTA files found in {path}\n"
                f"  Expected extensions: {', '.join(_FASTA_SUFFIXES)}",
                file=sys.stderr,
            )
            sys.exit(1)
        return files

    print(f"ERROR: input is neither a file nor a directory: {path}", file=sys.stderr)
    sys.exit(1)


def _sample_stem(fasta_path: Path) -> str:
    """Return the sample name by stripping all known FASTA suffixes."""
    name = fasta_path.name
    for s in _FASTA_SUFFIXES:
        if name.endswith(s):
            return name[: -len(s)]
    return fasta_path.stem


def _get_contig_lengths(fasta_path: Path) -> dict[str, int]:
    """Return a dict mapping contig name -> sequence length."""
    return {header.split()[0]: len(seq) for header, seq in read_fasta(fasta_path)}


def _write_hits_tsv(hits, output_dir: Path, sample: str) -> Path:
    out_path = output_dir / f"{sample}_hits.tsv"
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


def _write_verdicts_tsv(summaries, output_dir: Path, sample: str) -> Path:
    out_path = output_dir / f"{sample}_contig_verdicts.tsv"
    fieldnames = [
        "contig", "contig_length", "total_hits",
        "engineered_hits", "context_dependent_hits", "weak_hits", "unannotated_hits",
        "unique_feature_types", "unique_labels",
        "covered_bp", "covered_fraction",
        "top_label", "top_tier", "top_confidence",
        "evidence_summary",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow({
                "contig":                   s.contig,
                "contig_length":            s.contig_length,
                "total_hits":               s.total_hits,
                "engineered_hits":          s.engineered_hits,
                "context_dependent_hits":   s.context_dependent_hits,
                "weak_hits":                s.weak_hits,
                "unannotated_hits":         s.unannotated_hits,
                "unique_feature_types":     s.unique_feature_types,
                "unique_labels":            s.unique_labels,
                "covered_bp":               s.covered_bp,
                "covered_fraction":         f"{s.covered_fraction:.6f}",
                "top_label":                s.top_label,
                "top_tier":                 s.top_tier,
                "top_confidence":           s.top_confidence,
                "evidence_summary":         s.evidence_summary,
            })
    return out_path


def _process_one(
    input_path: Path,
    output_dir: Path,
    catalog_dir: Path,
    min_identity: float,
    min_coverage: float,
    threads: int,
    verbose: bool,
) -> None:
    """Run the full pipeline on a single FASTA file."""
    sample = _sample_stem(input_path)
    contig_lengths = _get_contig_lengths(input_path)

    hits = scan(
        query_fasta=input_path,
        catalog_dir=catalog_dir,
        min_identity=min_identity,
        min_coverage=min_coverage,
        threads=threads,
        verbose=verbose,
    )

    hits_path     = _write_hits_tsv(hits, output_dir, sample)
    summaries     = summarize(hits, contig_lengths)
    verdicts_path = _write_verdicts_tsv(summaries, output_dir, sample)

    print(f"  hits    : {len(hits):,} -> {hits_path.name}")
    print(f"  contigs : {len(summaries):,} with hits -> {verdicts_path.name}")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vectrap",
        description="Detect and classify synthetic vector contamination in genome assemblies.",
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="FASTA|DIR",
        help="Input FASTA file or directory of FASTA files.",
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

    fasta_files = _collect_inputs(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    if len(fasta_files) == 1:
        print(f"Processing: {fasta_files[0].name}")
        _process_one(
            input_path=fasta_files[0],
            output_dir=output_dir,
            catalog_dir=catalog_dir,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            threads=args.threads,
            verbose=args.verbose,
        )
    else:
        print(f"Found {len(fasta_files)} FASTA files in {input_path}")
        for i, fasta in enumerate(fasta_files, 1):
            print(f"[{i}/{len(fasta_files)}] {fasta.name}")
            _process_one(
                input_path=fasta,
                output_dir=output_dir,
                catalog_dir=catalog_dir,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                threads=args.threads,
                verbose=args.verbose,
            )

    print("Done.")


if __name__ == "__main__":
    main()

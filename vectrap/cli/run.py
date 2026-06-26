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
import time
from datetime import datetime
from pathlib import Path

from vectrap.modules.homology_scanner import scan
from vectrap.modules.scorer import summarize
from vectrap.modules.utils import read_fasta


_DEFAULT_CATALOG_DIR = Path(__file__).resolve().parents[2] / "vectrap" / "catalogs"

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")


def _ts() -> str:
    return datetime.now().strftime("[%H:%M:%S]")


def _elapsed(t0: float) -> str:
    secs = time.time() - t0
    if secs < 60:
        return f"{secs:.1f}s"
    m, s = divmod(int(secs), 60)
    return f"{m}m{s:02d}s"


def _collect_inputs(path: Path) -> list[Path]:
    """Return a list of FASTA files to process."""
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
    name = fasta_path.name
    for s in _FASTA_SUFFIXES:
        if name.endswith(s):
            return name[: -len(s)]
    return fasta_path.stem


def _get_contig_lengths(fasta_path: Path) -> dict[str, int]:
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
    min_len: int,
    best_n: int,
    threads: int,
    verbose: bool,
    label: str = "",
) -> None:
    """Run the full pipeline on a single FASTA file."""
    sample = _sample_stem(input_path)
    t0 = time.time()

    contig_lengths = _get_contig_lengths(input_path)

    hits = scan(
        query_fasta=input_path,
        catalog_dir=catalog_dir,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_len=min_len,
        best_n=best_n,
        threads=threads,
        verbose=verbose,
    )

    _write_hits_tsv(hits, output_dir, sample)
    summaries = summarize(hits, contig_lengths)
    _write_verdicts_tsv(summaries, output_dir, sample)

    prefix = f"{label} " if label else ""
    print(
        f"{_ts()} {prefix}{sample}  "
        f"hits={len(hits):,}  contigs_with_hits={len(summaries):,}  "
        f"elapsed={_elapsed(t0)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vectrap",
        description="Detect and classify synthetic vector contamination in genome assemblies.",
    )
    parser.add_argument("-i", "--input", required=True, metavar="FASTA|DIR",
                        help="Input FASTA file or directory of FASTA files.")
    parser.add_argument("-o", "--output", required=True, metavar="DIR",
                        help="Output directory. Created if it does not exist.")
    parser.add_argument("-c", "--catalog-dir", metavar="DIR",
                        default=str(_DEFAULT_CATALOG_DIR),
                        help=("Path to the catalog directory containing indexes built by "
                              "vectrap-build-db. Defaults to the bundled catalog directory "
                              f"({_DEFAULT_CATALOG_DIR})."))
    parser.add_argument("--min-identity", type=float, default=0.90, metavar="FLOAT",
                        help="Minimum sequence identity for mappy hits (default: 0.90).")
    parser.add_argument("--min-coverage", type=float, default=0.80, metavar="FLOAT",
                        help="Minimum fraction of catalog sequence covered by a hit (default: 0.80).")
    parser.add_argument("--min-len", type=int, default=50, metavar="INT",
                        help=("Minimum sequence length (bp) to use the mappy aligner; "
                              "shorter sequences use exact k-mer matching. "
                              "Must match the value used when building the catalog indexes (default: 50)."))
    parser.add_argument("--best-n", type=int, default=10, metavar="INT",
                        help="Maximum number of mappy alignments reported per query sequence (default: 10).")
    parser.add_argument("-t", "--threads", type=int, default=4, metavar="INT",
                        help="Number of threads for mappy aligner (default: 4).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print detailed per-stage progress messages to stderr.")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")
    args = parser.parse_args()

    input_path  = Path(args.input).resolve()
    output_dir  = Path(args.output).resolve()
    catalog_dir = Path(args.catalog_dir).resolve()

    fasta_files = _collect_inputs(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    n = len(fasta_files)
    if n > 1:
        print(f"{_ts()} Found {n} FASTA files in {input_path}")

    for i, fasta in enumerate(fasta_files, 1):
        label = f"[{i}/{n}]" if n > 1 else ""
        _process_one(
            input_path=fasta,
            output_dir=output_dir,
            catalog_dir=catalog_dir,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            min_len=args.min_len,
            best_n=args.best_n,
            threads=args.threads,
            verbose=args.verbose,
            label=label,
        )

    print(f"{_ts()} Done.")


if __name__ == "__main__":
    main()

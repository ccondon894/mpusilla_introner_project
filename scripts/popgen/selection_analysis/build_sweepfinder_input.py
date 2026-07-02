#!/usr/bin/env python3
"""Build SweepFinder2 allele-frequency and grid inputs from a Group 1 4D VCF."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Group 1 4D all-sites VCF (plain gzip).")
    parser.add_argument("--fai", required=True, help="Reference FAI for contig lengths.")
    parser.add_argument("--samples", required=True, help="Comma-separated Group 1 sample names.")
    parser.add_argument("--output-dir", required=True, help="Directory for per-contig freq/grid files.")
    parser.add_argument("--combined-freq", required=True, help="Genome-wide combined freq file for SF2 -f.")
    parser.add_argument("--contigs-tsv", required=True, help="Contig table: contig, length, freq_file, grid_file.")
    parser.add_argument("--grid-spacing", type=int, default=1000, help="Uniform grid spacing in bp.")
    parser.add_argument(
        "--min-called-haplotypes",
        type=int,
        default=8,
        help="Minimum called haplotypes required to retain a polymorphic site.",
    )
    parser.add_argument("--mating-contig", default="", help="Mating-type contig to exclude from grids.")
    parser.add_argument("--mating-start", type=int, default=0)
    parser.add_argument("--mating-end", type=int, default=0)
    parser.add_argument(
        "--exclude-contigs",
        default="",
        help="Comma-separated contig names to skip entirely (e.g. mating-type scaffold).",
    )
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def open_text(path: str | Path):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_field: str) -> int | None:
    gt = sample_field if len(sample_field) == 1 else sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|.", ""}:
        return None
    alt_count = 0
    for token in gt.replace("|", "/").split("/"):
        if token in {".", ""}:
            return None
        if token != "0":
            alt_count += 1
    return alt_count


def sanitize_contig(contig: str) -> str:
    return contig.replace("#", "_").replace("/", "_").replace(" ", "_")


def read_fai(fai_path: str | Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            contig, length, *_rest = line.rstrip("\n").split("\t")
            lengths[contig] = int(length)
    return lengths


def overlaps_interval(contig: str, pos: int, query_contig: str, start: int, end: int) -> bool:
    return contig == query_contig and start <= pos <= end


def write_freq_rows(rows: list[tuple[int, int, int, int]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("position\tx\tn\tfolded\n")
        for position, x, n, folded in rows:
            handle.write(f"{position}\t{x}\t{n}\t{folded}\n")


def write_grid_file(
    contig: str,
    length: int,
    grid_spacing: int,
    out_path: Path,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    positions = []
    pos = grid_spacing
    while pos <= length:
        if not overlaps_interval(contig, pos, mating_contig, mating_start, mating_end):
            positions.append(pos)
        pos += grid_spacing
    with open(out_path, "w", encoding="utf-8") as handle:
        for position in positions:
            handle.write(f"{position}\n")
    return len(positions)


def stream_vcf_polymorphisms(
    vcf_path: str | Path,
    samples: list[str],
    min_called: int,
) -> tuple[dict[str, list[tuple[int, int, int, int]]], dict[str, int]]:
    by_contig: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    stats = {
        "records_seen": 0,
        "records_used": 0,
        "skipped_missing": 0,
        "skipped_multiallelic": 0,
        "skipped_monomorphic": 0,
        "skipped_low_call_rate": 0,
        "skipped_no_sample_header": 0,
    }
    sample_cols: list[int] | None = None

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [sample for sample in samples if sample not in vcf_samples]
                if missing:
                    raise ValueError(f"VCF missing requested samples: {missing}")
                sample_cols = [9 + vcf_samples.index(sample) for sample in samples]
                continue
            if not line or line.startswith("#"):
                continue
            if sample_cols is None:
                stats["skipped_no_sample_header"] += 1
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sample_cols):
                continue
            stats["records_seen"] += 1

            alt = fields[4]
            if "," in alt or alt in {".", ""}:
                if alt not in {".", ""}:
                    stats["skipped_multiallelic"] += 1
                continue

            genotypes = []
            missing_gt = False
            for col_idx in sample_cols:
                gt = parse_gt(fields[col_idx])
                if gt is None:
                    missing_gt = True
                    break
                genotypes.append(gt)
            if missing_gt:
                stats["skipped_missing"] += 1
                continue

            n_called = len(genotypes)
            if n_called < min_called:
                stats["skipped_low_call_rate"] += 1
                continue

            alt_count = int(sum(genotypes))
            minor = min(alt_count, n_called - alt_count)
            if minor <= 0:
                stats["skipped_monomorphic"] += 1
                continue

            contig = fields[0]
            pos = int(fields[1])
            by_contig[contig].append((pos, minor, n_called, 1))
            stats["records_used"] += 1

    for contig in by_contig:
        by_contig[contig].sort(key=lambda row: row[0])
    return by_contig, stats


def main() -> None:
    args = parse_args()
    samples = parse_csv(args.samples)
    exclude_contigs = set(parse_csv(args.exclude_contigs))
    output_dir = Path(args.output_dir)
    freq_dir = output_dir / "freq"
    grid_dir = output_dir / "grid"
    freq_dir.mkdir(parents=True, exist_ok=True)
    grid_dir.mkdir(parents=True, exist_ok=True)

    contig_lengths = read_fai(args.fai)
    polymorphisms, vcf_stats = stream_vcf_polymorphisms(
        args.vcf, samples, args.min_called_haplotypes
    )

    combined_rows: list[tuple[int, int, int, int]] = []
    contig_table_rows = []

    for contig, length in sorted(contig_lengths.items(), key=lambda item: item[0]):
        if contig in exclude_contigs:
            continue
        rows = polymorphisms.get(contig, [])
        if not rows:
            continue

        safe = sanitize_contig(contig)
        freq_path = freq_dir / f"{safe}.freq"
        grid_path = grid_dir / f"{safe}.grid"
        write_freq_rows(rows, freq_path)
        n_grid = write_grid_file(
            contig,
            length,
            args.grid_spacing,
            grid_path,
            args.mating_contig,
            args.mating_start,
            args.mating_end,
        )
        if n_grid == 0:
            continue

        combined_rows.extend(rows)
        contig_table_rows.append(
            {
                "contig": contig,
                "length": length,
                "n_polymorphic_sites": len(rows),
                "n_grid_points": n_grid,
                "freq_file": str(freq_path),
                "grid_file": str(grid_path),
            }
        )

    combined_rows.sort(key=lambda row: (row[0], row[1]))
    write_freq_rows(combined_rows, Path(args.combined_freq))

    contigs_df_path = Path(args.contigs_tsv)
    contigs_df_path.parent.mkdir(parents=True, exist_ok=True)
    with open(contigs_df_path, "w", encoding="utf-8") as handle:
        handle.write("contig\tlength\tn_polymorphic_sites\tn_grid_points\tfreq_file\tgrid_file\n")
        for row in contig_table_rows:
            handle.write(
                f"{row['contig']}\t{row['length']}\t{row['n_polymorphic_sites']}\t"
                f"{row['n_grid_points']}\t{row['freq_file']}\t{row['grid_file']}\n"
            )

    print("VCF parse stats:")
    for key, value in vcf_stats.items():
        print(f"{key}\t{value}")
    print(f"Wrote {len(contig_table_rows)} contigs with polymorphic sites and grid points")
    print(f"Combined polymorphic sites: {len(combined_rows)}")


if __name__ == "__main__":
    main()

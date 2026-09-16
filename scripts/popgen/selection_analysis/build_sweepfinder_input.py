#!/usr/bin/env python3
"""Build filtered SweepFinder2 frequency, grid, and pyrho map inputs."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CANONICAL_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class PyrhoMap:
    starts: np.ndarray
    ends: np.ndarray
    rates: np.ndarray
    cumulative_morgans: np.ndarray
    map_start: int
    map_end: int
    usable_blocks: tuple[tuple[int, int], ...]
    masked_intervals: tuple[tuple[int, int, int, int, float], ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Group 1 all-sites VCF (plain gzip).")
    parser.add_argument("--fai", required=True, help="Reference FAI for contig lengths.")
    parser.add_argument("--pyrho-dir", required=True, help="Directory containing *.pyrho.out maps.")
    parser.add_argument(
        "--min-pyrho-rate",
        type=float,
        default=1e-20,
        help="Mask pyrho intervals at or below this numerical rate floor.",
    )
    parser.add_argument(
        "--pyrho-mask-buffer",
        type=int,
        default=50000,
        help="Physical bp masked on each side of a rate-floor interval.",
    )
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
    parser.add_argument(
        "--accepted-filters",
        default=".,PASS",
        help="Comma-separated VCF FILTER values to retain (default: . and PASS).",
    )
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def open_text(path: str | Path):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_field: str) -> tuple[int, int] | None:
    gt = sample_field if len(sample_field) == 1 else sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|.", ""}:
        return None
    alt_count = 0
    called = 0
    for token in gt.replace("|", "/").split("/"):
        if token in {".", ""}:
            return None
        if token not in {"0", "1"}:
            return None
        alt_count += int(token)
        called += 1
    return alt_count, called


def normalize_contig(contig: str) -> str:
    for token in contig.split("#"):
        if token.startswith("scaffold_"):
            return token
    return contig


def find_masked_intervals(
    starts: np.ndarray,
    ends: np.ndarray,
    rates: np.ndarray,
    threshold: float,
    buffer_bp: int,
) -> tuple[
    tuple[tuple[int, int, int, int, float], ...],
    tuple[tuple[int, int], ...],
]:
    """Find buffered rate-floor runs and all remaining continuous map blocks."""
    floor = rates <= threshold
    masked: list[tuple[int, int, int, int, float]] = []
    index = 0
    while index < len(rates):
        if not floor[index]:
            index += 1
            continue
        run_start_index = index
        while index + 1 < len(rates) and floor[index + 1]:
            index += 1
        run_end_index = index
        raw_start = int(starts[run_start_index])
        raw_end = int(ends[run_end_index])
        buffered_start = max(int(starts[0]), raw_start - buffer_bp)
        buffered_end = min(int(ends[-1]), raw_end + buffer_bp)
        masked.append(
            (
                raw_start,
                raw_end,
                buffered_start,
                buffered_end,
                float(rates[run_start_index : run_end_index + 1].min()),
            )
        )
        index += 1

    if not masked:
        return tuple(), ((int(starts[0]), int(ends[-1])),)

    merged: list[list[int]] = []
    for _raw_start, _raw_end, buffered_start, buffered_end, _min_rate in masked:
        if not merged or buffered_start > merged[-1][1]:
            merged.append([buffered_start, buffered_end])
        else:
            merged[-1][1] = max(merged[-1][1], buffered_end)

    usable_blocks: list[tuple[int, int]] = []
    cursor = int(starts[0])
    for masked_start, masked_end in merged:
        if cursor < masked_start:
            usable_blocks.append((cursor, masked_start))
        cursor = max(cursor, masked_end)
    if cursor < int(ends[-1]):
        usable_blocks.append((cursor, int(ends[-1])))

    if not usable_blocks:
        return tuple(masked), tuple()
    return tuple(masked), tuple(usable_blocks)


def load_pyrho_maps(
    pyrho_dir: str | Path,
    min_rate: float,
    mask_buffer: int,
) -> dict[str, PyrhoMap]:
    maps: dict[str, PyrhoMap] = {}
    for path in sorted(Path(pyrho_dir).glob("*.pyrho.out")):
        values = np.loadtxt(path, dtype=float)
        if values.ndim == 1:
            values = values.reshape(1, -1)
        if values.shape[1] < 3:
            raise ValueError(f"Malformed pyrho map: {path}")
        starts = values[:, 0].astype(np.int64)
        ends = values[:, 1].astype(np.int64)
        rates = values[:, 2].astype(float)
        order = np.lexsort((ends, starts))
        starts, ends, rates = starts[order], ends[order], rates[order]
        if np.any(ends <= starts) or np.any(rates < 0):
            raise ValueError(f"Invalid intervals or negative rates in {path}")
        if len(starts) > 1 and np.any(starts[1:] != ends[:-1]):
            raise ValueError(f"Non-contiguous pyrho intervals in {path}")
        interval_morgans = (
            (ends - starts).astype(np.longdouble)
            * rates.astype(np.longdouble)
        )
        cumulative = np.concatenate(
            [
                np.array([0.0], dtype=np.longdouble),
                np.cumsum(interval_morgans, dtype=np.longdouble),
            ]
        )
        masked_intervals, usable_blocks = find_masked_intervals(
            starts,
            ends,
            rates,
            min_rate,
            mask_buffer,
        )
        maps[path.name.removesuffix(".pyrho.out")] = PyrhoMap(
            starts=starts,
            ends=ends,
            rates=rates,
            cumulative_morgans=cumulative,
            map_start=int(starts[0]),
            map_end=int(ends[-1]),
            usable_blocks=usable_blocks,
            masked_intervals=masked_intervals,
        )
    if not maps:
        raise FileNotFoundError(f"No *.pyrho.out files in {pyrho_dir}")
    return maps


def cumulative_morgans_at(rate_map: PyrhoMap, coordinates: np.ndarray) -> np.ndarray:
    """Return integrated recombination probability from map start to coordinates."""
    x = np.asarray(coordinates, dtype=np.int64)
    result = np.zeros(len(x), dtype=np.longdouble)
    above = x >= rate_map.map_end
    result[above] = rate_map.cumulative_morgans[-1]
    inside = (x > rate_map.map_start) & (x < rate_map.map_end)
    if inside.any():
        xi = x[inside]
        index = np.searchsorted(rate_map.starts, xi, side="right") - 1
        index = np.clip(index, 0, len(rate_map.starts) - 1)
        covered = np.clip(
            xi - rate_map.starts[index],
            0,
            rate_map.ends[index] - rate_map.starts[index],
        )
        result[inside] = (
            rate_map.cumulative_morgans[index]
            + covered.astype(np.longdouble) * rate_map.rates[index]
        )
    return result


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


def write_recombination_rows(
    rows: list[tuple[int, int, int, int]],
    rate_map: PyrhoMap,
    out_path: Path,
) -> None:
    """Write cM since the preceding frequency site, as required by SF2."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    positions = np.asarray([row[0] for row in rows], dtype=np.int64)
    zero_based = positions - 1
    cumulative = cumulative_morgans_at(rate_map, zero_based)
    interval_cm = np.zeros(len(rows), dtype=np.longdouble)
    if len(rows) > 1:
        interval_cm[1:] = np.diff(cumulative) * 100.0
    if np.any(interval_cm < 0):
        raise ValueError("Integrated pyrho map produced a negative genetic distance")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("position\trate\n")
        for position, genetic_distance_cm in zip(positions, interval_cm):
            handle.write(f"{position}\t{float(genetic_distance_cm):.17g}\n")


def write_grid_file(
    contig: str,
    min_position: int,
    max_position: int,
    grid_spacing: int,
    out_path: Path,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    positions = []
    pos = max(grid_spacing, ((min_position + grid_spacing - 1) // grid_spacing) * grid_spacing)
    while pos <= max_position:
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
    accepted_filters: set[str],
    pyrho_maps: dict[str, PyrhoMap],
    exclude_contigs: set[str],
) -> tuple[dict[tuple[str, int], list[tuple[int, int, int, int]]], dict[str, int]]:
    by_contig_block: dict[tuple[str, int], list[tuple[int, int, int, int]]] = defaultdict(list)
    stats = {
        "records_seen": 0,
        "records_used": 0,
        "skipped_excluded_contig": 0,
        "skipped_filter": 0,
        "missing_genotypes": 0,
        "records_with_missing_genotypes": 0,
        "skipped_multiallelic": 0,
        "skipped_non_snp": 0,
        "skipped_monomorphic": 0,
        "skipped_low_call_rate": 0,
        "skipped_no_pyrho_map": 0,
        "skipped_outside_pyrho_map": 0,
        "skipped_masked_pyrho_interval": 0,
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

            contig = fields[0]
            if contig in exclude_contigs:
                stats["skipped_excluded_contig"] += 1
                continue

            if fields[6] not in accepted_filters:
                stats["skipped_filter"] += 1
                continue

            ref = fields[3].upper()
            alt = fields[4]
            if "," in alt or alt in {".", ""}:
                if alt not in {".", ""}:
                    stats["skipped_multiallelic"] += 1
                continue
            alt = alt.upper()
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in CANONICAL_BASES
                or alt not in CANONICAL_BASES
            ):
                stats["skipped_non_snp"] += 1
                continue

            map_key = normalize_contig(contig)
            rate_map = pyrho_maps.get(map_key)
            if rate_map is None:
                stats["skipped_no_pyrho_map"] += 1
                continue
            pos = int(fields[1])
            zero_based = pos - 1
            if zero_based < rate_map.map_start or zero_based >= rate_map.map_end:
                stats["skipped_outside_pyrho_map"] += 1
                continue
            block_index = next(
                (
                    index
                    for index, (block_start, block_end) in enumerate(rate_map.usable_blocks)
                    if block_start <= zero_based < block_end
                ),
                None,
            )
            if block_index is None:
                stats["skipped_masked_pyrho_interval"] += 1
                continue

            alt_count = 0
            n_called = 0
            record_has_missing = False
            for col_idx in sample_cols:
                parsed = parse_gt(fields[col_idx])
                if parsed is None:
                    stats["missing_genotypes"] += 1
                    record_has_missing = True
                    continue
                sample_alt, sample_called = parsed
                alt_count += sample_alt
                n_called += sample_called
            if record_has_missing:
                stats["records_with_missing_genotypes"] += 1
            if n_called < min_called:
                stats["skipped_low_call_rate"] += 1
                continue

            minor = min(alt_count, n_called - alt_count)
            if minor <= 0:
                stats["skipped_monomorphic"] += 1
                continue

            by_contig_block[(contig, block_index)].append((pos, minor, n_called, 1))
            stats["records_used"] += 1

    for key in by_contig_block:
        by_contig_block[key].sort(key=lambda row: row[0])
    return by_contig_block, stats


def main() -> None:
    args = parse_args()
    if args.min_pyrho_rate < 0:
        raise ValueError("--min-pyrho-rate must be nonnegative")
    if args.pyrho_mask_buffer < 0:
        raise ValueError("--pyrho-mask-buffer must be nonnegative")
    samples = parse_csv(args.samples)
    exclude_contigs = set(parse_csv(args.exclude_contigs))
    accepted_filters = set(parse_csv(args.accepted_filters))
    pyrho_maps = load_pyrho_maps(
        args.pyrho_dir,
        args.min_pyrho_rate,
        args.pyrho_mask_buffer,
    )
    output_dir = Path(args.output_dir)
    freq_dir = output_dir / "freq"
    grid_dir = output_dir / "grid"
    freq_dir.mkdir(parents=True, exist_ok=True)
    grid_dir.mkdir(parents=True, exist_ok=True)
    mask_summary_path = output_dir / "pyrho_masked_intervals.tsv"
    with open(mask_summary_path, "w", encoding="utf-8") as handle:
        handle.write(
            "contig\traw_start\traw_end\tbuffered_start\tbuffered_end\t"
            "minimum_rate\tusable_blocks\n"
        )
        for contig, rate_map in sorted(pyrho_maps.items()):
            for (
                raw_start,
                raw_end,
                buffered_start,
                buffered_end,
                min_rate,
            ) in rate_map.masked_intervals:
                usable_blocks = ";".join(
                    f"{start}-{end}" for start, end in rate_map.usable_blocks
                )
                handle.write(
                    f"{contig}\t{raw_start}\t{raw_end}\t{buffered_start}\t"
                    f"{buffered_end}\t{min_rate}\t{usable_blocks}\n"
                )

    contig_lengths = read_fai(args.fai)
    polymorphisms, vcf_stats = stream_vcf_polymorphisms(
        args.vcf,
        samples,
        args.min_called_haplotypes,
        accepted_filters,
        pyrho_maps,
        exclude_contigs,
    )

    combined_rows: list[tuple[int, int, int, int]] = []
    contig_table_rows = []

    for contig, length in sorted(contig_lengths.items(), key=lambda item: item[0]):
        if contig in exclude_contigs:
            continue
        map_key = normalize_contig(contig)
        if map_key not in pyrho_maps:
            continue
        rate_map = pyrho_maps[map_key]
        if not rate_map.usable_blocks:
            continue
        safe_contig = sanitize_contig(contig)
        multiple_blocks = len(rate_map.usable_blocks) > 1
        for block_index, _block in enumerate(rate_map.usable_blocks):
            rows = polymorphisms.get((contig, block_index), [])
            if not rows:
                continue
            safe = (
                f"{safe_contig}_block_{block_index + 1}"
                if multiple_blocks
                else safe_contig
            )
            freq_path = freq_dir / f"{safe}.freq"
            grid_path = grid_dir / f"{safe}.grid"
            recombination_path = output_dir / "recombination" / f"{safe}.rec"
            write_freq_rows(rows, freq_path)
            write_recombination_rows(rows, rate_map, recombination_path)
            n_grid = write_grid_file(
                contig,
                rows[0][0],
                rows[-1][0],
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
                    "recombination_file": str(recombination_path),
                    "analysis_start": rows[0][0],
                    "analysis_end": rows[-1][0],
                    "block_index": block_index + 1,
                }
            )

    combined_rows.sort(key=lambda row: (row[0], row[1]))
    write_freq_rows(combined_rows, Path(args.combined_freq))

    contigs_df_path = Path(args.contigs_tsv)
    contigs_df_path.parent.mkdir(parents=True, exist_ok=True)
    with open(contigs_df_path, "w", encoding="utf-8") as handle:
        handle.write(
            "contig\tlength\tn_polymorphic_sites\tn_grid_points\t"
            "freq_file\tgrid_file\trecombination_file\tanalysis_start\t"
            "analysis_end\tblock_index\n"
        )
        for row in contig_table_rows:
            handle.write(
                f"{row['contig']}\t{row['length']}\t{row['n_polymorphic_sites']}\t"
                f"{row['n_grid_points']}\t{row['freq_file']}\t{row['grid_file']}\t"
                f"{row['recombination_file']}\t{row['analysis_start']}\t"
                f"{row['analysis_end']}\t{row['block_index']}\n"
            )

    print("VCF parse stats:")
    for key, value in vcf_stats.items():
        print(f"{key}\t{value}")
    print(
        f"Wrote {len(contig_table_rows)} analysis blocks across "
        f"{len({row['contig'] for row in contig_table_rows})} contigs"
    )
    print(f"Combined polymorphic sites: {len(combined_rows)}")
    print(f"Accepted FILTER values: {','.join(sorted(accepted_filters))}")
    print("Recombination distances: integrated pyrho r converted to centimorgans")
    print(f"Masked pyrho rate threshold: <= {args.min_pyrho_rate}")
    print(f"Pyrho mask boundary buffer (bp): {args.pyrho_mask_buffer}")
    print("Masked pyrho intervals (0-based, half-open):")
    for contig, rate_map in sorted(pyrho_maps.items()):
        for raw_start, raw_end, buffered_start, buffered_end, min_rate in rate_map.masked_intervals:
            print(
                "pyrho_mask\t"
                f"{contig}\t{raw_start}\t{raw_end}\t{buffered_start}\t"
                f"{buffered_end}\t{min_rate}"
            )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Score RNA-seq splice-junction support and retention PSI for introners."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median


THRESHOLDS = (0, 1, 2, 5, 10)


@dataclass(frozen=True)
class ExpectedInterval:
    source: str
    start: int
    end: int


@dataclass(frozen=True)
class Junction:
    replicate: str
    contig: str
    start: int
    end: int
    score: int
    strand: str
    name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--replicates", required=True, nargs="+")
    parser.add_argument("--junction-beds", required=True, nargs="+", type=Path)
    parser.add_argument("--bams", required=True, nargs="+", type=Path)
    parser.add_argument("--per-locus", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--plot", required=True, type=Path)
    parser.add_argument("--flank-length", type=int, default=100)
    parser.add_argument("--max-loci", type=int, default=0)
    parser.add_argument("--min-retained-aligned-bp", type=int, default=8)
    parser.add_argument("--low-psi-threshold", type=float, default=0.10)
    parser.add_argument("--mating-contig", default="")
    parser.add_argument("--mating-start", type=int, default=0)
    parser.add_argument("--mating-end", type=int, default=0)
    return parser.parse_args()


def parse_splice_site(splice_site: str) -> tuple[str, int | None, int | None]:
    if not splice_site:
        return "", None, None

    donor_type = ""
    donor_offset = None
    donor_match = re.search(r"\b(GT|GC)@(\d+)", splice_site)
    if donor_match:
        donor_type = donor_match.group(1)
        donor_offset = int(donor_match.group(2))

    acceptor_offset = None
    acceptor_matches = list(re.finditer(r"\bAG@(\d+)", splice_site))
    if acceptor_matches:
        acceptor_offset = int(acceptor_matches[-1].group(1))

    return donor_type, donor_offset, acceptor_offset


def add_expected_interval(
    intervals: list[ExpectedInterval],
    seen: set[tuple[int, int, str]],
    source: str,
    start: int,
    end: int,
    body_start: int,
    body_end: int,
) -> None:
    if body_start <= start < end <= body_end:
        key = (start, end, source)
        if key not in seen:
            intervals.append(ExpectedInterval(source, start, end))
            seen.add(key)


def expected_intervals(
    body_start: int,
    body_end: int,
    splice_site: str,
) -> tuple[list[ExpectedInterval], str, int | None, int | None]:
    donor_type, donor_offset, acceptor_offset = parse_splice_site(splice_site)
    intervals: list[ExpectedInterval] = []
    seen: set[tuple[int, int, str]] = set()

    add_expected_interval(
        intervals, seen, "body", body_start, body_end, body_start, body_end
    )

    if donor_offset is not None and acceptor_offset is not None:
        plus_start = body_start + donor_offset
        plus_end = body_start + acceptor_offset + 2
        add_expected_interval(
            intervals,
            seen,
            "splice_site_plus",
            plus_start,
            plus_end,
            body_start,
            body_end,
        )

        rc_start = body_end - acceptor_offset - 2
        rc_end = body_end - donor_offset
        add_expected_interval(
            intervals,
            seen,
            "splice_site_reverse_complement",
            rc_start,
            rc_end,
            body_start,
            body_end,
        )

    return intervals, donor_type, donor_offset, acceptor_offset


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def should_exclude_mating_type(row: dict[str, str], args: argparse.Namespace) -> bool:
    if args.sample != "CCMP1545" or not args.mating_contig:
        return False
    return (
        row["contig"] == args.mating_contig
        and overlaps(
            int(row["start"]), int(row["end"]), args.mating_start, args.mating_end
        )
    )


def load_loci(args: argparse.Namespace) -> list[dict[str, str]]:
    loci: list[dict[str, str]] = []
    with args.matrix.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["sample"] != args.sample or row["presence"] != "1":
                continue
            if should_exclude_mating_type(row, args):
                continue
            start = int(row["start"])
            end = int(row["end"])
            body_start = start + args.flank_length
            body_end = end - args.flank_length
            if body_end <= body_start:
                continue
            row = dict(row)
            row["body_start"] = str(body_start)
            row["body_end"] = str(body_end)
            row["body_len"] = str(body_end - body_start)
            loci.append(row)
            if args.max_loci and len(loci) >= args.max_loci:
                break
    return loci


def parse_regtools_bed(path: Path, replicate: str) -> list[Junction]:
    junctions: list[Junction] = []
    if not path.exists():
        return junctions

    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            chrom_start = int(fields[1])
            chrom_end = int(fields[2])
            score = int(float(fields[4]))
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            if len(block_sizes) < 2:
                continue
            junction_start = chrom_start + block_sizes[0]
            junction_end = chrom_end - block_sizes[1]
            if junction_start >= junction_end:
                continue

            junctions.append(
                Junction(
                    replicate=replicate,
                    contig=fields[0],
                    start=junction_start,
                    end=junction_end,
                    score=score,
                    strand=fields[5],
                    name=fields[3],
                )
            )
    return junctions


def load_junctions(
    replicates: list[str],
    junction_beds: list[Path],
) -> dict[str, list[Junction]]:
    junctions_by_contig: dict[str, list[Junction]] = defaultdict(list)
    for replicate, path in zip(replicates, junction_beds):
        for junction in parse_regtools_bed(path, replicate):
            junctions_by_contig[junction.contig].append(junction)
    return junctions_by_contig


def boundary_distance(
    expected: ExpectedInterval, junction: Junction
) -> tuple[int, int, int, int]:
    start_delta = junction.start - expected.start
    end_delta = junction.end - expected.end
    max_abs_delta = max(abs(start_delta), abs(end_delta))
    sum_abs_delta = abs(start_delta) + abs(end_delta)
    return start_delta, end_delta, max_abs_delta, sum_abs_delta


def choose_best_expected(
    intervals: list[ExpectedInterval],
    junctions: list[Junction],
) -> tuple[
    ExpectedInterval,
    Junction | None,
    tuple[int | str, int | str, int | str, int | str],
]:
    if not junctions:
        return intervals[0], None, ("", "", "", "")

    best = None
    for expected in intervals:
        for junction in junctions:
            start_delta, end_delta, max_abs_delta, sum_abs_delta = boundary_distance(
                expected, junction
            )
            candidate = (
                max_abs_delta,
                sum_abs_delta,
                -junction.score,
                expected.source != "body",
                expected.source,
                expected.start,
                expected.end,
                junction.replicate,
                junction.name,
                expected,
                junction,
                (start_delta, end_delta, max_abs_delta, sum_abs_delta),
            )
            if best is None or candidate < best:
                best = candidate

    assert best is not None
    return best[9], best[10], best[11]


def junction_support_by_threshold(
    expected: ExpectedInterval,
    junctions: list[Junction],
) -> tuple[dict[int, int], dict[int, int]]:
    support = {threshold: 0 for threshold in THRESHOLDS}
    replicate_sets = {threshold: set() for threshold in THRESHOLDS}
    for junction in junctions:
        _, _, max_abs_delta, _ = boundary_distance(expected, junction)
        for threshold in THRESHOLDS:
            if max_abs_delta <= threshold:
                support[threshold] += junction.score
                if junction.score > 0:
                    replicate_sets[threshold].add(junction.replicate)
    return support, {threshold: len(values) for threshold, values in replicate_sets.items()}


def count_retained_signal(
    bam_path: Path,
    contig: str,
    start: int,
    end: int,
    min_aligned_bp: int,
) -> tuple[int, int]:
    import pysam

    retained_reads = 0
    retained_aligned_bases = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if contig not in bam.references:
            return 0, 0
        for read in bam.fetch(contig, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            aligned_bases = 0
            for block_start, block_end in read.get_blocks():
                if block_end <= start:
                    continue
                if block_start >= end:
                    break
                aligned_bases += max(0, min(block_end, end) - max(block_start, start))
            if aligned_bases >= min_aligned_bp:
                retained_reads += 1
                retained_aligned_bases += aligned_bases
    return retained_reads, retained_aligned_bases


def fmt_float(value: float | str, digits: int = 6) -> str:
    if value == "" or value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return f"{float(value):.{digits}f}"


def score_loci(
    loci: list[dict[str, str]],
    junctions_by_contig: dict[str, list[Junction]],
    args: argparse.Namespace,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    bam_by_replicate = dict(zip(args.replicates, args.bams))

    for locus in loci:
        contig = locus["contig"]
        body_start = int(locus["body_start"])
        body_end = int(locus["body_end"])
        body_len = int(locus["body_len"])
        splice_site = locus.get("splice_site", "")

        intervals, donor_type, donor_offset, acceptor_offset = expected_intervals(
            body_start, body_end, splice_site
        )
        contig_junctions = junctions_by_contig.get(contig, [])
        best_expected, best_junction, deltas = choose_best_expected(
            intervals, contig_junctions
        )
        threshold_support, threshold_replicates = junction_support_by_threshold(
            best_expected, contig_junctions
        )

        retained_by_rep = {}
        retained_bases_by_rep = {}
        for replicate, bam_path in bam_by_replicate.items():
            retained_reads, retained_bases = count_retained_signal(
                bam_path, contig, body_start, body_end, args.min_retained_aligned_bp
            )
            retained_by_rep[replicate] = retained_reads
            retained_bases_by_rep[replicate] = retained_bases

        retained_signal = sum(retained_by_rep.values())
        retained_aligned_bases = sum(retained_bases_by_rep.values())
        retained_mean_depth = retained_aligned_bases / body_len if body_len else 0.0
        psi_junction_signal = threshold_support[2]
        total_psi_signal = retained_signal + psi_junction_signal
        retention_psi = (
            retained_signal / total_psi_signal if total_psi_signal > 0 else ""
        )
        spliced_fraction = (
            psi_junction_signal / total_psi_signal if total_psi_signal > 0 else ""
        )

        high_conf = (
            psi_junction_signal > 0
            and threshold_replicates[2] >= 2
            and retention_psi != ""
            and float(retention_psi) <= args.low_psi_threshold
        )

        start_delta, end_delta, max_abs_delta, sum_abs_delta = deltas
        row = {
            "ortholog_id": locus["ortholog_id"],
            "sample": locus["sample"],
            "contig": contig,
            "raw_start": locus["start"],
            "raw_end": locus["end"],
            "body_start": str(body_start),
            "body_end": str(body_end),
            "body_len": str(body_len),
            "family": locus.get("family", ""),
            "gene": locus.get("gene", ""),
            "orientation": locus.get("orientation", ""),
            "matrix_splice_site": splice_site,
            "splice_donor_type": donor_type,
            "splice_donor_offset": "" if donor_offset is None else str(donor_offset),
            "splice_acceptor_offset": (
                "" if acceptor_offset is None else str(acceptor_offset)
            ),
            "expected_source": best_expected.source,
            "expected_splice_start": str(best_expected.start),
            "expected_splice_end": str(best_expected.end),
            "best_junction_replicate": (
                "" if best_junction is None else best_junction.replicate
            ),
            "best_junction_name": "" if best_junction is None else best_junction.name,
            "best_junction_start": (
                "" if best_junction is None else str(best_junction.start)
            ),
            "best_junction_end": (
                "" if best_junction is None else str(best_junction.end)
            ),
            "best_junction_score": (
                "" if best_junction is None else str(best_junction.score)
            ),
            "best_junction_strand": (
                "" if best_junction is None else best_junction.strand
            ),
            "start_delta": "" if start_delta == "" else str(start_delta),
            "end_delta": "" if end_delta == "" else str(end_delta),
            "max_abs_boundary_delta": "" if max_abs_delta == "" else str(max_abs_delta),
            "sum_abs_boundary_delta": "" if sum_abs_delta == "" else str(sum_abs_delta),
            "retained_signal": str(retained_signal),
            "retained_aligned_bases": str(retained_aligned_bases),
            "retained_mean_depth": fmt_float(retained_mean_depth),
            "psi_junction_signal_within_2bp": str(psi_junction_signal),
            "retention_psi": fmt_float(retention_psi),
            "spliced_fraction": fmt_float(spliced_fraction),
            "covered_for_psi": "yes" if total_psi_signal > 0 else "no",
            "low_retention_psi": (
                "yes"
                if retention_psi != ""
                and float(retention_psi) <= args.low_psi_threshold
                else "no"
            ),
            "high_confidence_efficient_splicing": "yes" if high_conf else "no",
        }

        for threshold in THRESHOLDS:
            label = "exact" if threshold == 0 else f"within_{threshold}bp"
            row[f"{label}_junction_support"] = str(threshold_support[threshold])
            row[f"{label}_replicate_count"] = str(threshold_replicates[threshold])

        rows.append(row)

    return rows


def write_per_locus(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def safe_median(values: list[float]) -> str:
    return fmt_float(median(values)) if values else ""


def summarize_group(group: str, rows: list[dict[str, str]]) -> dict[str, str]:
    n_loci = len(rows)
    covered = [row for row in rows if row["covered_for_psi"] == "yes"]
    psi_values = [float(row["retention_psi"]) for row in covered if row["retention_psi"]]
    retained_values = [float(row["retained_signal"]) for row in rows]
    junction_values = [float(row["psi_junction_signal_within_2bp"]) for row in rows]
    low_psi_count = sum(row["low_retention_psi"] == "yes" for row in rows)
    high_conf_count = sum(
        row["high_confidence_efficient_splicing"] == "yes" for row in rows
    )

    out = {
        "group": group,
        "n_loci": str(n_loci),
        "covered_for_psi_loci": str(len(covered)),
        "median_retention_psi": safe_median(psi_values),
        "median_retained_signal": safe_median(retained_values),
        "median_junction_support_within_2bp": safe_median(junction_values),
        "low_retention_psi_loci": str(low_psi_count),
        "low_retention_psi_fraction_of_covered": fmt_float(
            low_psi_count / len(covered) if covered else ""
        ),
        "high_confidence_efficient_splicing_loci": str(high_conf_count),
        "high_confidence_efficient_splicing_fraction": fmt_float(
            high_conf_count / n_loci if n_loci else ""
        ),
    }

    for threshold in THRESHOLDS:
        label = "exact" if threshold == 0 else f"within_{threshold}bp"
        key = f"{label}_junction_support"
        count = sum(int(row[key]) > 0 for row in rows)
        out[f"{label}_supported_loci"] = str(count)
        out[f"{label}_supported_fraction"] = fmt_float(count / n_loci if n_loci else "")

    any_count = sum(row["best_junction_start"] != "" for row in rows)
    out["any_junction_on_contig_loci"] = str(any_count)
    out["any_junction_on_contig_fraction"] = fmt_float(any_count / n_loci if n_loci else "")
    return out


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    groups: list[tuple[str, list[dict[str, str]]]] = [("all", rows)]
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_family[row["family"] or "NA"].append(row)
    groups.extend(
        (f"family_{family}", fam_rows)
        for family, fam_rows in sorted(by_family.items())
    )

    summary_rows = [summarize_group(group, group_rows) for group, group_rows in groups]
    fieldnames = list(summary_rows[0].keys()) if summary_rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


def write_plot(path: Path, sample: str, rows: list[dict[str, str]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    labels = ["exact", "within_1bp", "within_2bp", "within_5bp"]
    fractions = []
    for label in labels:
        key = f"{label}_junction_support"
        fractions.append(
            sum(int(row[key]) > 0 for row in rows) / len(rows) if rows else 0.0
        )

    axes[0].bar(["exact", "<=1 bp", "<=2 bp", "<=5 bp"], fractions, color="#4C78A8")
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Fraction of present introners")
    axes[0].set_title("RNA junction boundary support")

    psi_values = [float(row["retention_psi"]) for row in rows if row["retention_psi"]]
    if psi_values:
        axes[1].hist(
            psi_values, bins=30, range=(0, 1), color="#59A14F", edgecolor="white"
        )
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Retention PSI")
    axes[1].set_ylabel("Introner loci")
    axes[1].set_title("Retention PSI among covered loci")

    fig.suptitle(f"{sample} introner splice-boundary verification")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if len(args.replicates) != len(args.junction_beds):
        raise SystemExit("--replicates and --junction-beds must have the same length")
    if len(args.replicates) != len(args.bams):
        raise SystemExit("--replicates and --bams must have the same length")

    loci = load_loci(args)
    junctions_by_contig = load_junctions(args.replicates, args.junction_beds)
    rows = score_loci(loci, junctions_by_contig, args)
    write_per_locus(args.per_locus, rows)
    write_summary(args.summary, rows)
    write_plot(args.plot, args.sample, rows)


if __name__ == "__main__":
    main()

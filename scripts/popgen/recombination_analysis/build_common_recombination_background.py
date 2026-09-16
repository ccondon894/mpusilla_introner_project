#!/usr/bin/env python3
"""Build one exonic recombination background shared by three focal classes.

The focal classes are:
  1. all introner-centered windows,
  2. polymorphic introner-centered windows, and
  3. Group 1 fixed introner-centered windows, and
  4. Group 1 polymorphic canonical-intron-centered windows.

Overlap among focal classes is retained. A candidate exonic background window
is retained only when it is disjoint from the union of all focal windows and
from the CCMP1545 mating-type region. Recombination-rate values are carried
through without applying plotting-only rate cutoffs.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_GROUP1 = [
    "CCMP1545",
    "RCC114",
    "RCC1614",
    "RCC1698",
    "RCC2482",
    "RCC373",
    "RCC465",
    "RCC629",
    "RCC692",
    "RCC693",
    "RCC833",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-introner-windows", required=True)
    parser.add_argument("--polymorphic-introner-windows", required=True)
    parser.add_argument("--introner-genotype-matrix", required=True)
    parser.add_argument("--canonical-intron-matrix", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument(
        "--pyrho-dir",
        required=True,
        help="Directory containing CCMP1545 *.pyrho.out recombination maps.",
    )
    parser.add_argument(
        "--group1-samples",
        default=",".join(DEFAULT_GROUP1),
        help="Comma-separated Group 1 samples used to define polymorphism.",
    )
    parser.add_argument("--window-size", type=int, default=5000)
    parser.add_argument("--merge-distance", type=int, default=10000)
    parser.add_argument(
        "--mating-contig", default="CCMP1545#0#scaffold_2"
    )
    parser.add_argument("--mating-start", type=int, default=49808)
    parser.add_argument("--mating-end", type=int, default=1730591)
    parser.add_argument("--output-background", required=True)
    parser.add_argument("--output-excluded-background", required=True)
    parser.add_argument("--output-fixed-introner-focals", required=True)
    parser.add_argument("--output-canonical-focals", required=True)
    parser.add_argument("--output-comparison-windows", required=True)
    parser.add_argument("--output-audit", required=True)
    return parser.parse_args()


def normalize_contig(contig: str) -> str:
    """Normalize graph-style names such as CCMP1545#0#scaffold_2."""
    return str(contig).split("#")[-1]


def intervals_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Return overlap for zero-based, half-open intervals."""
    return start1 < end2 and start2 < end1


def parse_gtf_attribute(attributes: str, key: str) -> str | None:
    pattern = rf"(?:^|;\s*){re.escape(key)} \"([^\"]+)\""
    match = re.search(pattern, attributes)
    return match.group(1) if match else None


def load_transcript_to_gene(gtf_path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            transcript_id = parse_gtf_attribute(fields[8], "transcript_id")
            gene_id = parse_gtf_attribute(fields[8], "gene_id")
            if transcript_id and gene_id:
                mapping[transcript_id] = gene_id
    return mapping


def load_gtf_gene_exons(
    gtf_path: str,
) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = defaultdict(
        lambda: {"contig": None, "exons": []}
    )
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            gene_id = parse_gtf_attribute(fields[8], "gene_id")
            if gene_id is None:
                continue
            genes[gene_id]["contig"] = normalize_contig(fields[0])
            genes[gene_id]["exons"].append(
                (int(fields[3]) - 1, int(fields[4]))
            )

    for gene_data in genes.values():
        merged: list[list[int]] = []
        for start, end in sorted(gene_data["exons"]):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        gene_data["exons"] = [(start, end) for start, end in merged]
    return dict(genes)


def overlaps_mating_region(
    contig: str,
    start: int,
    end: int,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> bool:
    return (
        normalize_contig(contig) == normalize_contig(mating_contig)
        and intervals_overlap(start, end, mating_start, mating_end)
    )


def load_existing_focals(
    path: str,
    focal_class: str,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, int]:
    frame = pd.read_csv(path, sep="\t")
    required = {"contig", "start", "end", "type"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    frame = frame.loc[frame["type"].eq("introner_containing")].copy()
    frame["contig"] = frame["contig"].map(normalize_contig)
    frame["start"] = pd.to_numeric(frame["start"], errors="raise").astype(int)
    frame["end"] = pd.to_numeric(frame["end"], errors="raise").astype(int)
    frame["focal_class"] = focal_class

    mating_mask = frame.apply(
        lambda row: overlaps_mating_region(
            row["contig"],
            int(row["start"]),
            int(row["end"]),
            mating_contig,
            mating_start,
            mating_end,
        ),
        axis=1,
    )
    n_mating = int(mating_mask.sum())
    frame = frame.loc[~mating_mask].reset_index(drop=True)
    return frame, n_mating


def load_polymorphic_canonical_loci(
    matrix_path: str,
    group1_samples: list[str],
    transcript_to_gene: dict[str, str],
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = pd.read_csv(matrix_path, sep="\t")
    required = {
        "gene_id",
        "contig",
        "ref_start",
        "ref_end",
        "intron_index",
        *group1_samples,
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{matrix_path} is missing columns: {sorted(missing)}")

    genotypes = frame[group1_samples].apply(pd.to_numeric, errors="raise")
    complete = genotypes.ne(3).all(axis=1)
    polymorphic = genotypes.eq(1).any(axis=1) & genotypes.eq(2).any(axis=1)
    selected = frame.loc[complete & polymorphic].copy()

    # The matrix's gene_id field contains transcript IDs (the source catalog is
    # transcript-based), so retain that identity explicitly and map it to the
    # actual GTF gene_id.
    selected = selected.rename(columns={"gene_id": "transcript_id"})
    selected["gene_id"] = selected["transcript_id"].map(transcript_to_gene)
    selected["gene_id"] = selected["gene_id"].fillna(selected["transcript_id"])
    selected["contig"] = selected["contig"].map(normalize_contig)

    # Matrix reference coordinates are 1-based inclusive. Convert to the
    # zero-based, half-open convention used by recombination windows.
    selected["start"] = selected["ref_start"].astype(int) - 1
    selected["end"] = selected["ref_end"].astype(int)
    selected["locus_id"] = (
        selected["transcript_id"].astype(str)
        + ":"
        + selected["intron_index"].astype(str)
    )

    mating_mask = selected.apply(
        lambda row: overlaps_mating_region(
            row["contig"],
            int(row["start"]),
            int(row["end"]),
            mating_contig,
            mating_start,
            mating_end,
        ),
        axis=1,
    )
    n_mating = int(mating_mask.sum())
    selected = selected.loc[~mating_mask].copy()

    before_dedup = len(selected)
    selected = selected.drop_duplicates(
        ["contig", "start", "end", "gene_id"]
    ).reset_index(drop=True)
    audit = {
        "canonical_matrix_loci": len(frame),
        "canonical_complete_group1_loci": int(complete.sum()),
        "canonical_polymorphic_group1_loci": int((complete & polymorphic).sum()),
        "canonical_polymorphic_mating_loci_excluded": n_mating,
        "canonical_polymorphic_duplicate_loci_removed": before_dedup - len(selected),
        "canonical_polymorphic_unique_loci": len(selected),
    }
    return selected, audit


def load_fixed_group1_introner_loci(
    matrix_path: str,
    group1_samples: list[str],
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    required = {
        "ortholog_id",
        "sample",
        "contig",
        "start",
        "end",
        "presence",
        "group1_present_count",
        "group1_absent_count",
        "group1_missing_count",
        "group1_callable_count",
        "group1_pattern",
        "within_group_status",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{matrix_path} is missing columns: {sorted(missing)}")

    reference = frame.loc[frame["sample"].eq("CCMP1545")].copy()
    fixed_mask = (
        reference["presence"].eq(1)
        & reference["group1_present_count"].eq(len(group1_samples))
        & reference["group1_absent_count"].eq(0)
        & reference["group1_missing_count"].eq(0)
        & reference["group1_callable_count"].eq(len(group1_samples))
        & reference["group1_pattern"].eq("fixed_present")
    )
    fixed = reference.loc[fixed_mask].copy()
    fixed["contig"] = fixed["contig"].map(normalize_contig)
    fixed["start"] = fixed["start"].astype(int)
    fixed["end"] = fixed["end"].astype(int)

    mating_mask = fixed.apply(
        lambda row: overlaps_mating_region(
            row["contig"],
            int(row["start"]),
            int(row["end"]),
            mating_contig,
            mating_start,
            mating_end,
        ),
        axis=1,
    )
    n_mating = int(mating_mask.sum())
    fixed = fixed.loc[~mating_mask].copy()

    before_dedup = len(fixed)
    fixed = fixed.drop_duplicates(
        ["ortholog_id", "contig", "start", "end"]
    ).reset_index(drop=True)
    audit = {
        "fixed_introner_reference_loci_before_mating_filter": int(fixed_mask.sum()),
        "fixed_introner_mating_loci_excluded": n_mating,
        "fixed_introner_duplicate_loci_removed": before_dedup - len(fixed),
        "fixed_introner_unique_loci": len(fixed),
        "fixed_introner_low_identity_loci": int(
            fixed["within_group_status"].eq("low_identity").sum()
        ),
    }
    return fixed, audit


def cluster_fixed_introner_loci(
    loci: pd.DataFrame,
    gene_exons: dict[str, dict[str, object]],
    window_size: int,
    merge_distance: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    loci_by_contig: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in loci.to_dict("records"):
        loci_by_contig[str(record["contig"])].append(record)

    rows: list[dict[str, object]] = []
    assigned_orthologs: set[str] = set()
    n_gene_locus_assignments = 0
    for gene_id, gene_data in gene_exons.items():
        contig = str(gene_data["contig"])
        exons = gene_data["exons"]
        matched: list[dict[str, object]] = []
        for record in loci_by_contig.get(contig, []):
            if any(
                intervals_overlap(
                    int(record["start"]), int(record["end"]), exon_start, exon_end
                )
                for exon_start, exon_end in exons
            ):
                matched.append(record)
                assigned_orthologs.add(str(record["ortholog_id"]))
                n_gene_locus_assignments += 1
        if not matched:
            continue

        matched.sort(key=lambda record: (int(record["start"]), int(record["end"])))
        current = [matched[0]]
        clusters: list[list[dict[str, object]]] = []
        for record in matched[1:]:
            if int(record["start"]) - int(current[-1]["end"]) <= merge_distance:
                current.append(record)
            else:
                clusters.append(current)
                current = [record]
        clusters.append(current)

        for cluster in clusters:
            cluster_start = int(cluster[0]["start"])
            cluster_end = int(cluster[-1]["end"])
            midpoint = (cluster_start + cluster_end) // 2
            window_start = max(0, midpoint - window_size // 2)
            rows.append(
                {
                    "contig": contig,
                    "start": window_start,
                    "end": window_start + window_size,
                    "focal_class": "fixed_introner",
                    "cluster_midpoint": midpoint,
                    "cluster_locus_start": cluster_start,
                    "cluster_locus_end": cluster_end,
                    "n_loci": len(cluster),
                    "gene_id": gene_id,
                    "ortholog_ids": ",".join(
                        sorted({str(record["ortholog_id"]) for record in cluster})
                    ),
                }
            )

    focal = pd.DataFrame(rows).sort_values(
        ["contig", "start", "end", "gene_id"]
    ).reset_index(drop=True)
    focal.insert(
        0,
        "focal_id",
        [f"fixed_introner_{index:04d}" for index in range(1, len(focal) + 1)],
    )
    audit = {
        "fixed_introner_loci_assigned_to_exonic_gene": len(assigned_orthologs),
        "fixed_introner_loci_not_assigned_to_exonic_gene": len(loci) - len(assigned_orthologs),
        "fixed_introner_gene_locus_assignments": n_gene_locus_assignments,
        "fixed_introner_focal_windows": len(focal),
        "fixed_introner_multi_locus_focal_windows": int(focal["n_loci"].gt(1).sum()),
        "fixed_introner_max_loci_per_focal_window": int(focal["n_loci"].max()),
    }
    return focal, audit


def cluster_canonical_loci(
    loci: pd.DataFrame, window_size: int, merge_distance: int
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for (contig, gene_id), group in loci.groupby(["contig", "gene_id"]):
        records = group.sort_values(["start", "end"]).to_dict("records")
        current = [records[0]]
        clusters: list[list[dict[str, object]]] = []

        for record in records[1:]:
            if int(record["start"]) - int(current[-1]["end"]) <= merge_distance:
                current.append(record)
            else:
                clusters.append(current)
                current = [record]
        clusters.append(current)

        for cluster in clusters:
            cluster_start = int(cluster[0]["start"])
            cluster_end = int(cluster[-1]["end"])
            midpoint = (cluster_start + cluster_end) // 2
            window_start = max(0, midpoint - window_size // 2)
            window_end = window_start + window_size
            rows.append(
                {
                    "contig": contig,
                    "start": window_start,
                    "end": window_end,
                    "focal_class": "polymorphic_canonical_intron",
                    "cluster_midpoint": midpoint,
                    "cluster_locus_start": cluster_start,
                    "cluster_locus_end": cluster_end,
                    "n_loci": len(cluster),
                    "gene_id": gene_id,
                    "transcript_ids": ",".join(
                        sorted({str(record["transcript_id"]) for record in cluster})
                    ),
                    "locus_ids": ",".join(
                        sorted({str(record["locus_id"]) for record in cluster})
                    ),
                }
            )

    focal = pd.DataFrame(rows).sort_values(
        ["contig", "start", "end", "gene_id"]
    ).reset_index(drop=True)
    focal.insert(
        0,
        "focal_id",
        [f"polymorphic_canonical_{index:04d}" for index in range(1, len(focal) + 1)],
    )
    return focal


class IntervalLookup:
    def __init__(self, frame: pd.DataFrame):
        self.by_contig: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for row in frame.itertuples(index=False):
            self.by_contig[str(row.contig)].append((int(row.start), int(row.end)))
        for contig in self.by_contig:
            self.by_contig[contig].sort()

    def overlaps(self, contig: str, start: int, end: int) -> bool:
        for other_start, other_end in self.by_contig.get(str(contig), []):
            if other_start >= end:
                break
            if intervals_overlap(start, end, other_start, other_end):
                return True
        return False


def count_rows_overlapping(left: pd.DataFrame, right: pd.DataFrame) -> int:
    lookup = IntervalLookup(right)
    return sum(
        lookup.overlaps(str(row.contig), int(row.start), int(row.end))
        for row in left.itertuples(index=False)
    )


def count_internal_overlap_rows(frame: pd.DataFrame) -> int:
    hit_rows: set[int] = set()
    for _contig, group in frame.groupby("contig"):
        records = sorted(
            (int(row.start), int(row.end), int(index))
            for index, row in group.iterrows()
        )
        for position, (start, end, index) in enumerate(records):
            for other_start, other_end, other_index in records[position + 1 :]:
                if other_start >= end:
                    break
                if intervals_overlap(start, end, other_start, other_end):
                    hit_rows.update((index, other_index))
    return len(hit_rows)


def load_pyrho_maps(pyrho_dir: str) -> dict[str, pd.DataFrame]:
    maps: dict[str, pd.DataFrame] = {}
    for path in sorted(Path(pyrho_dir).glob("*.pyrho.out")):
        contig = normalize_contig(path.name.removesuffix(".pyrho.out"))
        maps[contig] = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["start", "end", "rate"],
        )
    if not maps:
        raise ValueError(f"No *.pyrho.out files found in {pyrho_dir}")
    return maps


def weighted_mean_rate(
    start: int, end: int, pyrho_windows: pd.DataFrame
) -> float:
    overlapping = pyrho_windows.loc[
        pyrho_windows["end"].gt(start) & pyrho_windows["start"].lt(end)
    ]
    if overlapping.empty:
        return np.nan

    overlap_starts = np.maximum(overlapping["start"].to_numpy(), start)
    overlap_ends = np.minimum(overlapping["end"].to_numpy(), end)
    overlap_lengths = overlap_ends - overlap_starts
    total_length = overlap_lengths.sum()
    if total_length <= 0:
        return np.nan
    return float(np.sum(overlapping["rate"].to_numpy() * overlap_lengths) / total_length)


def add_recombination_rates(
    focal_windows: pd.DataFrame, pyrho_maps: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    focal_windows = focal_windows.copy()
    rates: list[float] = []
    for row in focal_windows.itertuples(index=False):
        recombination_map = pyrho_maps.get(str(row.contig))
        if recombination_map is None:
            rates.append(np.nan)
            continue
        rates.append(
            weighted_mean_rate(int(row.start), int(row.end), recombination_map)
        )
    focal_windows["rate"] = rates
    return focal_windows


def build_common_background(
    all_window_path: str,
    all_focals: pd.DataFrame,
    polymorphic_focals: pd.DataFrame,
    fixed_focals: pd.DataFrame,
    canonical_focals: pd.DataFrame,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(all_window_path, sep="\t")
    required = {"contig", "start", "end", "type", "rate"}
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"{all_window_path} is missing columns: {sorted(missing)}")

    candidates = source.loc[source["type"].eq("non_introner_containing")].copy()
    candidates["contig"] = candidates["contig"].map(normalize_contig)
    candidates["start"] = candidates["start"].astype(int)
    candidates["end"] = candidates["end"].astype(int)

    lookups = {
        "overlap_all_introner_focal": IntervalLookup(all_focals),
        "overlap_polymorphic_introner_focal": IntervalLookup(polymorphic_focals),
        "overlap_fixed_introner_focal": IntervalLookup(fixed_focals),
        "overlap_polymorphic_canonical_focal": IntervalLookup(canonical_focals),
    }
    for column, lookup in lookups.items():
        candidates[column] = [
            lookup.overlaps(str(row.contig), int(row.start), int(row.end))
            for row in candidates.itertuples(index=False)
        ]

    candidates["overlap_mating_type"] = [
        overlaps_mating_region(
            str(row.contig),
            int(row.start),
            int(row.end),
            mating_contig,
            mating_start,
            mating_end,
        )
        for row in candidates.itertuples(index=False)
    ]

    exclusion_columns = [*lookups, "overlap_mating_type"]
    excluded_mask = candidates[exclusion_columns].any(axis=1)
    common = candidates.loc[~excluded_mask].copy()
    excluded = candidates.loc[excluded_mask].copy()

    common = common.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    excluded = excluded.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    common.insert(
        0,
        "background_id",
        [f"common_exonic_background_{index:04d}" for index in range(1, len(common) + 1)],
    )
    excluded.insert(
        0,
        "background_id",
        [f"excluded_exonic_background_{index:04d}" for index in range(1, len(excluded) + 1)],
    )
    return common, excluded


def write_audit(path: str, rows: list[tuple[str, object]]) -> None:
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(path, sep="\t", index=False)


def standardize_comparison_class(
    frame: pd.DataFrame, comparison_class: str, id_prefix: str
) -> pd.DataFrame:
    standardized = frame.copy().reset_index(drop=True)
    standardized["comparison_class"] = comparison_class
    standardized["window_id"] = [
        f"{id_prefix}_{index:04d}" for index in range(1, len(standardized) + 1)
    ]
    if "gene_id" not in standardized:
        standardized["gene_id"] = pd.NA
    standardized["window_length"] = standardized["end"] - standardized["start"]
    return standardized[
        [
            "comparison_class",
            "window_id",
            "contig",
            "start",
            "end",
            "window_length",
            "rate",
            "gene_id",
        ]
    ]


def main() -> None:
    args = parse_args()
    if args.window_size <= 0 or args.merge_distance < 0:
        raise ValueError("window-size must be positive and merge-distance nonnegative")

    group1_samples = [
        sample.strip() for sample in args.group1_samples.split(",") if sample.strip()
    ]
    transcript_to_gene = load_transcript_to_gene(args.gtf)
    gene_exons = load_gtf_gene_exons(args.gtf)

    all_focals, all_mating_excluded = load_existing_focals(
        args.all_introner_windows,
        "all_introner",
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )
    polymorphic_focals, polymorphic_mating_excluded = load_existing_focals(
        args.polymorphic_introner_windows,
        "polymorphic_introner",
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )
    fixed_loci, fixed_locus_audit = load_fixed_group1_introner_loci(
        args.introner_genotype_matrix,
        group1_samples,
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )
    fixed_focals, fixed_focal_audit = cluster_fixed_introner_loci(
        fixed_loci,
        gene_exons,
        args.window_size,
        args.merge_distance,
    )
    canonical_loci, canonical_audit = load_polymorphic_canonical_loci(
        args.canonical_intron_matrix,
        group1_samples,
        transcript_to_gene,
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )
    canonical_focals = cluster_canonical_loci(
        canonical_loci, args.window_size, args.merge_distance
    )
    pyrho_maps = load_pyrho_maps(args.pyrho_dir)
    fixed_focals = add_recombination_rates(fixed_focals, pyrho_maps)
    canonical_focals = add_recombination_rates(canonical_focals, pyrho_maps)

    common, excluded = build_common_background(
        args.all_introner_windows,
        all_focals,
        polymorphic_focals,
        fixed_focals,
        canonical_focals,
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )

    focal_union = pd.concat(
        [
            all_focals[["contig", "start", "end"]],
            polymorphic_focals[["contig", "start", "end"]],
            fixed_focals[["contig", "start", "end"]],
            canonical_focals[["contig", "start", "end"]],
        ],
        ignore_index=True,
    )
    if count_rows_overlapping(common, focal_union):
        raise RuntimeError("Common background still overlaps at least one focal window")
    if count_internal_overlap_rows(common):
        raise RuntimeError("Common background contains internally overlapping windows")

    comparison_windows = pd.concat(
        [
            standardize_comparison_class(
                common, "common_exonic_background", "common_background"
            ),
            standardize_comparison_class(
                all_focals, "all_introner", "all_introner"
            ),
            standardize_comparison_class(
                polymorphic_focals,
                "polymorphic_introner",
                "polymorphic_introner",
            ),
            standardize_comparison_class(
                fixed_focals,
                "fixed_introner",
                "fixed_introner",
            ),
            standardize_comparison_class(
                canonical_focals,
                "polymorphic_canonical_intron",
                "polymorphic_canonical",
            ),
        ],
        ignore_index=True,
    )

    audit_rows: list[tuple[str, object]] = [
        ("window_size", args.window_size),
        ("merge_distance", args.merge_distance),
        ("group1_samples", ",".join(group1_samples)),
        ("all_introner_focal_windows", len(all_focals)),
        ("all_introner_mating_windows_excluded", all_mating_excluded),
        ("polymorphic_introner_focal_windows", len(polymorphic_focals)),
        (
            "polymorphic_introner_mating_windows_excluded",
            polymorphic_mating_excluded,
        ),
        *fixed_locus_audit.items(),
        *fixed_focal_audit.items(),
        (
            "fixed_introner_focals_missing_rates",
            int(fixed_focals["rate"].isna().sum()),
        ),
        (
            "fixed_introner_focal_rates_below_1e-14",
            int((fixed_focals["rate"] < 1e-14).sum()),
        ),
        (
            "fixed_introner_focal_rates_above_5e-11",
            int((fixed_focals["rate"] > 5e-11).sum()),
        ),
        *canonical_audit.items(),
        ("polymorphic_canonical_focal_windows", len(canonical_focals)),
        (
            "polymorphic_canonical_multi_locus_focal_windows",
            int(canonical_focals["n_loci"].gt(1).sum()),
        ),
        (
            "polymorphic_canonical_max_loci_per_focal_window",
            int(canonical_focals["n_loci"].max()),
        ),
        (
            "polymorphic_canonical_focals_missing_rates",
            int(canonical_focals["rate"].isna().sum()),
        ),
        (
            "polymorphic_canonical_focal_rates_below_1e-14",
            int((canonical_focals["rate"] < 1e-14).sum()),
        ),
        (
            "polymorphic_canonical_focal_rates_above_5e-11",
            int((canonical_focals["rate"] > 5e-11).sum()),
        ),
        (
            "canonical_focals_overlapping_all_introner_focals",
            count_rows_overlapping(canonical_focals, all_focals),
        ),
        (
            "fixed_focals_overlapping_all_introner_focals",
            count_rows_overlapping(fixed_focals, all_focals),
        ),
        (
            "fixed_focals_overlapping_polymorphic_introner_focals",
            count_rows_overlapping(fixed_focals, polymorphic_focals),
        ),
        (
            "fixed_focals_overlapping_polymorphic_canonical_focals",
            count_rows_overlapping(fixed_focals, canonical_focals),
        ),
        (
            "canonical_focals_overlapping_polymorphic_introner_focals",
            count_rows_overlapping(canonical_focals, polymorphic_focals),
        ),
        (
            "polymorphic_introner_focals_overlapping_all_introner_focals",
            count_rows_overlapping(polymorphic_focals, all_focals),
        ),
        ("candidate_exonic_background_windows", len(common) + len(excluded)),
        (
            "candidate_background_overlapping_all_introner_focal",
            int(excluded["overlap_all_introner_focal"].sum()),
        ),
        (
            "candidate_background_overlapping_polymorphic_introner_focal",
            int(excluded["overlap_polymorphic_introner_focal"].sum()),
        ),
        (
            "candidate_background_overlapping_fixed_introner_focal",
            int(excluded["overlap_fixed_introner_focal"].sum()),
        ),
        (
            "candidate_background_overlapping_polymorphic_canonical_focal",
            int(excluded["overlap_polymorphic_canonical_focal"].sum()),
        ),
        (
            "candidate_background_overlapping_mating_type",
            int(excluded["overlap_mating_type"].sum()),
        ),
        ("excluded_background_windows_union", len(excluded)),
        ("common_background_windows", len(common)),
        (
            "common_background_full_5kb_windows",
            int((common["end"] - common["start"]).eq(args.window_size).sum()),
        ),
        (
            "common_background_shorter_windows",
            int((common["end"] - common["start"]).lt(args.window_size).sum()),
        ),
        (
            "common_background_median_window_length",
            float((common["end"] - common["start"]).median()),
        ),
        ("common_background_internal_overlap_rows", count_internal_overlap_rows(common)),
        ("common_background_focal_overlap_rows", count_rows_overlapping(common, focal_union)),
        ("common_background_missing_rates", int(pd.to_numeric(common["rate"], errors="coerce").isna().sum())),
        ("common_background_rates_below_1e-14", int((pd.to_numeric(common["rate"], errors="coerce") < 1e-14).sum())),
        ("common_background_rates_above_5e-11", int((pd.to_numeric(common["rate"], errors="coerce") > 5e-11).sum())),
    ]

    for output_path in [
        args.output_background,
        args.output_excluded_background,
        args.output_fixed_introner_focals,
        args.output_canonical_focals,
        args.output_comparison_windows,
        args.output_audit,
    ]:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    common.to_csv(args.output_background, sep="\t", index=False)
    excluded.to_csv(args.output_excluded_background, sep="\t", index=False)
    fixed_focals.to_csv(args.output_fixed_introner_focals, sep="\t", index=False)
    canonical_focals.to_csv(args.output_canonical_focals, sep="\t", index=False)
    comparison_windows.to_csv(
        args.output_comparison_windows, sep="\t", index=False
    )
    write_audit(args.output_audit, audit_rows)

    print(f"All-introner focal windows: {len(all_focals)}")
    print(f"Polymorphic-introner focal windows: {len(polymorphic_focals)}")
    print(f"Fixed-introner focal windows: {len(fixed_focals)}")
    print(f"Polymorphic canonical-intron focal windows: {len(canonical_focals)}")
    print(f"Candidate exonic background windows: {len(common) + len(excluded)}")
    print(f"Excluded background windows: {len(excluded)}")
    print(f"Common background windows: {len(common)}")
    print(f"Wrote {args.output_background}")
    print(f"Wrote {args.output_excluded_background}")
    print(f"Wrote {args.output_fixed_introner_focals}")
    print(f"Wrote {args.output_canonical_focals}")
    print(f"Wrote {args.output_comparison_windows}")
    print(f"Wrote {args.output_audit}")


if __name__ == "__main__":
    main()

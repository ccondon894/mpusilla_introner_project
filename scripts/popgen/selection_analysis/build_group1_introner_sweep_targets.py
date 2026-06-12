#!/usr/bin/env python3
"""Build CCMP1545-anchored Group 1 introner sweep targets."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--group1-samples", required=True)
    parser.add_argument("--group2-samples", required=True)
    parser.add_argument("--reference-sample", default="CCMP1545")
    parser.add_argument(
        "--target-mode",
        choices=[
            "group1_only_polymorphic",
            "group1_fixed_present_group2_absent",
            "group1_fixed_present_any_group2",
        ],
        default="group1_only_polymorphic",
        help="Introner target class to retain.",
    )
    parser.add_argument(
        "--accepted-within-statuses",
        default="consistent,singleton",
        help="Comma-separated within_group_status values to keep.",
    )
    parser.add_argument("--mating-contig", required=True)
    parser.add_argument("--mating-start", type=int, required=True)
    parser.add_argument("--mating-end", type=int, required=True)
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def frequency_category(present_count: int, n_samples: int) -> str:
    absent_count = n_samples - present_count
    if present_count == n_samples:
        return "fixed_present"
    if absent_count == n_samples:
        return "fixed_absent"
    if present_count == 1:
        return "singleton_present"
    if absent_count == 1:
        return "singleton_absent"
    if present_count == 2:
        return "doubleton_present"
    if absent_count == 2:
        return "doubleton_absent"
    if present_count <= n_samples // 3:
        return "low_frequency"
    if present_count >= (2 * n_samples) // 3:
        return "high_frequency"
    return "intermediate_frequency"


def overlaps_interval(
    contig: str, start: int, end: int, query_contig: str, query_start: int, query_end: int
) -> bool:
    return contig == query_contig and start <= query_end and end >= query_start


def build_targets(
    genotype_matrix: str | Path,
    group1_samples: list[str],
    group2_samples: list[str],
    reference_sample: str,
    target_mode: str,
    accepted_within_statuses: set[str],
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    df = pd.read_csv(genotype_matrix, sep="\t")
    required = {"ortholog_id", "sample", "presence", "contig", "start", "end"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Genotype matrix missing required columns: {missing}")

    if "within_group_status" not in df.columns:
        df["within_group_status"] = ""

    targets: list[dict[str, object]] = []
    counts = {
        "orthologs_seen": int(df["ortholog_id"].nunique()),
        "group1_complete_target_class_accepted": 0,
        "excluded_group2_present": 0,
        "excluded_no_reference_anchor": 0,
        "excluded_mating_type": 0,
        "targets_retained": 0,
    }

    for ortholog_id, group in df.groupby("ortholog_id", sort=True):
        g1 = group[group["sample"].isin(group1_samples)]
        g1_by_sample = g1.drop_duplicates("sample").set_index("sample")
        if not set(group1_samples).issubset(g1_by_sample.index):
            continue

        g1_presence = g1_by_sample.loc[group1_samples, "presence"]
        if not g1_presence.isin([1, 2]).all():
            continue
        present_count = int((g1_presence == 1).sum())
        absent_count = int((g1_presence == 2).sum())

        if target_mode == "group1_only_polymorphic":
            if present_count == 0 or absent_count == 0:
                continue
        elif target_mode in {"group1_fixed_present_group2_absent", "group1_fixed_present_any_group2"}:
            if present_count != len(group1_samples) or absent_count != 0:
                continue
        else:
            raise ValueError(f"Unknown target mode: {target_mode}")

        g2 = group[group["sample"].isin(group2_samples)]
        group2_present_count = int((g2["presence"] == 1).sum())
        if target_mode in {"group1_only_polymorphic", "group1_fixed_present_group2_absent"}:
            if group2_present_count > 0:
                counts["excluded_group2_present"] += 1
                continue

        if target_mode == "group1_fixed_present_group2_absent":
            group2_presence = g2["presence"]
            if not group2_presence.empty and not group2_presence.isin([2, 3]).all():
                counts["excluded_group2_present"] += 1
                continue

        if present_count == 0:
            continue

        within_values = {
            str(value)
            for value in group.get("within_group_status", pd.Series(dtype=object)).dropna().unique()
            if str(value)
        }
        if within_values and not (within_values & accepted_within_statuses):
            continue

        counts["group1_complete_target_class_accepted"] += 1

        ref_rows = group[group["sample"] == reference_sample]
        if ref_rows.empty:
            counts["excluded_no_reference_anchor"] += 1
            continue
        ref = ref_rows.iloc[0]
        try:
            start = int(ref["start"])
            end = int(ref["end"])
        except (TypeError, ValueError):
            counts["excluded_no_reference_anchor"] += 1
            continue
        contig = str(ref["contig"])
        if start < 0 or end < 0 or not contig:
            counts["excluded_no_reference_anchor"] += 1
            continue
        if end < start:
            start, end = end, start

        if overlaps_interval(contig, start, end, mating_contig, mating_start, mating_end):
            counts["excluded_mating_type"] += 1
            continue

        present_samples = [sample for sample in group1_samples if int(g1_by_sample.loc[sample, "presence"]) == 1]
        absent_samples = [sample for sample in group1_samples if int(g1_by_sample.loc[sample, "presence"]) == 2]
        midpoint = (start + end) // 2
        targets.append(
            {
                "ortholog_id": ortholog_id,
                "contig": contig,
                "start": start,
                "end": end,
                "midpoint": midpoint,
                "length_bp": end - start + 1,
                "reference_sample": reference_sample,
                "reference_presence": int(ref["presence"]) if pd.notna(ref["presence"]) else None,
                "group1_present_count": present_count,
                "group1_absent_count": absent_count,
                "group1_n_samples": len(group1_samples),
                "group1_present_frequency": present_count / len(group1_samples),
                "frequency_category": frequency_category(present_count, len(group1_samples)),
                "present_samples": ",".join(present_samples),
                "absent_samples": ",".join(absent_samples),
                "group2_present_count": group2_present_count,
                "group2_samples_checked": ",".join(group2_samples),
                "within_group_status": ";".join(sorted(within_values)) if within_values else "",
                "primary_set": target_mode,
            }
        )

    targets_df = pd.DataFrame(targets)
    if not targets_df.empty:
        targets_df = targets_df.sort_values(["contig", "midpoint", "ortholog_id"]).reset_index(drop=True)
    counts["targets_retained"] = int(len(targets_df))
    return targets_df, counts


def write_bed(targets: pd.DataFrame, output_bed: str | Path) -> None:
    columns = ["contig", "bed_start", "end", "ortholog_id", "group1_present_count", "midpoint"]
    if targets.empty:
        Path(output_bed).write_text("")
        return
    bed = targets.copy()
    bed["bed_start"] = (bed["start"].astype(int) - 1).clip(lower=0)
    bed[columns].to_csv(output_bed, sep="\t", header=False, index=False)


def main() -> None:
    args = parse_args()
    output_tsv = Path(args.output_tsv)
    output_bed = Path(args.output_bed)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    targets, counts = build_targets(
        genotype_matrix=args.genotype_matrix,
        group1_samples=parse_csv(args.group1_samples),
        group2_samples=parse_csv(args.group2_samples),
        reference_sample=args.reference_sample,
        target_mode=args.target_mode,
        accepted_within_statuses=set(parse_csv(args.accepted_within_statuses)),
        mating_contig=args.mating_contig,
        mating_start=args.mating_start,
        mating_end=args.mating_end,
    )
    targets.to_csv(output_tsv, sep="\t", index=False)
    write_bed(targets, output_bed)

    for key, value in counts.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()

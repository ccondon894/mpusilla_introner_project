#!/usr/bin/env python3
"""Build a clean RCC1749-specific exonic background for Group 2 introners."""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group2-windows", required=True)
    parser.add_argument(
        "--mating-contig", default="RCC1749#0#intronerless_contig_28"
    )
    parser.add_argument("--mating-start", type=int, default=25000)
    parser.add_argument("--mating-end", type=int, default=2148000)
    parser.add_argument("--output-background", required=True)
    parser.add_argument("--output-excluded-background", required=True)
    parser.add_argument("--output-focals", required=True)
    parser.add_argument("--output-comparison-windows", required=True)
    parser.add_argument("--output-audit", required=True)
    return parser.parse_args()


def normalize_contig(contig: str) -> str:
    return str(contig).split("#")[-1]


def intervals_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 < end2 and start2 < end1


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


def write_audit(path: str, rows: list[tuple[str, object]]) -> None:
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(
        path, sep="\t", index=False
    )


def main() -> None:
    args = parse_args()
    # Preserve the source rate strings. Converting them to floats and writing
    # them back can collapse nearly identical values into artificial ties,
    # which can shift a rank statistic even when the numerical difference is
    # far below biological precision.
    frame = pd.read_csv(args.group2_windows, sep="\t", dtype={"rate": "string"})
    required = {"contig", "start", "end", "type", "rate"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Group 2 window table missing columns: {sorted(missing)}")

    frame = frame.copy()
    frame["contig"] = frame["contig"].map(normalize_contig)
    frame["start"] = frame["start"].astype(int)
    frame["end"] = frame["end"].astype(int)
    mating_contig = normalize_contig(args.mating_contig)
    mating_mask = frame["contig"].eq(mating_contig) & frame["start"].lt(args.mating_end) & frame["end"].gt(args.mating_start - 1)

    focal_type = "introner_containing"
    background_type = "non_introner_containing"
    focal_mating_excluded = int((mating_mask & frame["type"].eq(focal_type)).sum())
    background_mating_excluded = int(
        (mating_mask & frame["type"].eq(background_type)).sum()
    )

    focals = frame.loc[
        frame["type"].eq(focal_type) & ~mating_mask
    ].copy().reset_index(drop=True)
    candidates = frame.loc[frame["type"].eq(background_type)].copy()
    focal_lookup = IntervalLookup(focals)
    candidates["overlap_mating_type"] = candidates["contig"].eq(mating_contig) & candidates["start"].lt(args.mating_end) & candidates["end"].gt(args.mating_start - 1)
    candidates["overlap_population2_focal"] = [
        focal_lookup.overlaps(str(row.contig), int(row.start), int(row.end))
        for row in candidates.itertuples(index=False)
    ]

    excluded_mask = candidates[
        ["overlap_mating_type", "overlap_population2_focal"]
    ].any(axis=1)
    background = candidates.loc[~excluded_mask].copy()
    excluded = candidates.loc[excluded_mask].copy()

    background = background.sort_values(["contig", "start", "end"]).reset_index(
        drop=True
    )
    excluded = excluded.sort_values(["contig", "start", "end"]).reset_index(
        drop=True
    )
    focals = focals.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    background.insert(
        0,
        "background_id",
        [f"population2_exonic_background_{i:04d}" for i in range(1, len(background) + 1)],
    )
    excluded.insert(
        0,
        "background_id",
        [f"excluded_population2_background_{i:04d}" for i in range(1, len(excluded) + 1)],
    )
    focals.insert(
        0,
        "focal_id",
        [f"population2_introner_{i:04d}" for i in range(1, len(focals) + 1)],
    )

    remaining_focal_collisions = sum(
        focal_lookup.overlaps(str(row.contig), int(row.start), int(row.end))
        for row in background.itertuples(index=False)
    )
    internal_overlap_rows = count_internal_overlap_rows(background)
    if remaining_focal_collisions:
        raise RuntimeError("Population 2 background still overlaps focal windows")
    if internal_overlap_rows:
        raise RuntimeError("Population 2 background contains internal overlaps")

    comparison = pd.concat([background, focals], ignore_index=True, sort=False)
    background_lengths = background["end"] - background["start"]
    audit_rows = [
        ("input_group2_windows", len(frame)),
        ("population2_focal_windows", len(focals)),
        ("population2_focal_mating_windows_excluded", focal_mating_excluded),
        ("population2_background_candidates_before_mating_filter", len(candidates)),
        ("population2_background_mating_windows_excluded", background_mating_excluded),
        (
            "population2_background_candidates_after_mating_filter",
            len(candidates) - background_mating_excluded,
        ),
        (
            "population2_background_overlapping_focal",
            int(
                (
                    candidates["overlap_population2_focal"]
                    & ~candidates["overlap_mating_type"]
                ).sum()
            ),
        ),
        ("population2_excluded_background_union", len(excluded)),
        ("population2_clean_background_windows", len(background)),
        (
            "population2_clean_background_full_5kb_windows",
            int(background_lengths.eq(5000).sum()),
        ),
        (
            "population2_clean_background_shorter_windows",
            int(background_lengths.lt(5000).sum()),
        ),
        (
            "population2_clean_background_median_window_length",
            float(background_lengths.median()),
        ),
        ("population2_clean_background_internal_overlap_rows", internal_overlap_rows),
        (
            "population2_clean_background_focal_overlap_rows",
            remaining_focal_collisions,
        ),
        (
            "population2_clean_background_missing_rates",
            int(pd.to_numeric(background["rate"], errors="coerce").isna().sum()),
        ),
        (
            "population2_focal_missing_rates",
            int(pd.to_numeric(focals["rate"], errors="coerce").isna().sum()),
        ),
    ]

    for path in [
        args.output_background,
        args.output_excluded_background,
        args.output_focals,
        args.output_comparison_windows,
        args.output_audit,
    ]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    background.to_csv(args.output_background, sep="\t", index=False)
    excluded.to_csv(args.output_excluded_background, sep="\t", index=False)
    focals.to_csv(args.output_focals, sep="\t", index=False)
    comparison.to_csv(args.output_comparison_windows, sep="\t", index=False)
    write_audit(args.output_audit, audit_rows)

    print(f"Population 2 focal windows: {len(focals)}")
    print(f"Population 2 clean background windows: {len(background)}")
    print(f"Excluded background windows: {len(excluded)}")
    print(f"Wrote {args.output_background}")
    print(f"Wrote {args.output_excluded_background}")
    print(f"Wrote {args.output_focals}")
    print(f"Wrote {args.output_comparison_windows}")
    print(f"Wrote {args.output_audit}")


if __name__ == "__main__":
    main()

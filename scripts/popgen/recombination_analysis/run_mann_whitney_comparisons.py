#!/usr/bin/env python3
"""Run unfiltered recombination-rate Mann-Whitney tests for five focal classes."""

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group1-comparison-windows", required=True)
    parser.add_argument("--group2-windows", required=True)
    parser.add_argument(
        "--group2-exclude-contig",
        default="intronerless_contig_28",
        help="Mating-type contig removed from Group 2 before testing.",
    )
    parser.add_argument("--group2-mating-start", type=int, default=25000)
    parser.add_argument("--group2-mating-end", type=int, default=2148000)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def test_row(
    comparison: str,
    population: str,
    focal_class: str,
    background_class: str,
    focal: pd.Series,
    background: pd.Series,
) -> dict[str, object]:
    focal = pd.to_numeric(focal, errors="coerce").dropna()
    background = pd.to_numeric(background, errors="coerce").dropna()
    if focal.empty or background.empty:
        raise ValueError(f"Insufficient data for {comparison}")

    statistic, p_value = mannwhitneyu(
        focal, background, alternative="two-sided"
    )
    median_focal = float(focal.median())
    median_background = float(background.median())
    fold_change = median_focal / median_background
    rank_biserial = 2 * statistic / (len(focal) * len(background)) - 1
    return {
        "comparison": comparison,
        "population": population,
        "focal_class": focal_class,
        "background_class": background_class,
        "n_focal": len(focal),
        "n_background": len(background),
        "median_focal_rate": median_focal,
        "median_background_rate": median_background,
        "median_fold_change": fold_change,
        "median_elevation_percent": (fold_change - 1) * 100,
        "u_statistic_focal": statistic,
        "rank_biserial_correlation": rank_biserial,
        "p_value_raw": p_value,
    }


def main() -> None:
    args = parse_args()
    group1 = pd.read_csv(args.group1_comparison_windows, sep="\t")
    group2 = pd.read_csv(args.group2_windows, sep="\t")

    required_group1 = {"comparison_class", "rate"}
    required_group2 = {"type", "rate"}
    if missing := required_group1 - set(group1.columns):
        raise ValueError(f"Group 1 table missing columns: {sorted(missing)}")
    if missing := required_group2 - set(group2.columns):
        raise ValueError(f"Group 2 table missing columns: {sorted(missing)}")

    if "contig" not in group2.columns:
        raise ValueError("Group 2 table missing contig column")
    normalized_group2_contigs = group2["contig"].astype(str).str.split("#").str[-1]
    excluded_contig = str(args.group2_exclude_contig).split("#")[-1]
    mating_mask = normalized_group2_contigs.eq(excluded_contig) & group2.start.lt(args.group2_mating_end) & group2.end.gt(args.group2_mating_start - 1)
    n_group2_mating_rows = int(mating_mask.sum())
    group2 = group2.loc[~mating_mask].copy()

    common_background = group1.loc[
        group1["comparison_class"].eq("common_exonic_background"), "rate"
    ]
    rows = []
    for label, focal_class in [
        ("Population 1 all introners", "all_introner"),
        ("Population 1 fixed introners", "fixed_introner"),
        ("Population 1 polymorphic introners", "polymorphic_introner"),
    ]:
        rows.append(
            test_row(
                label,
                "Population 1",
                focal_class,
                "common_exonic_background",
                group1.loc[group1["comparison_class"].eq(focal_class), "rate"],
                common_background,
            )
        )

    rows.append(
        test_row(
            "Population 2 all introners",
            "Population 2",
            "all_introner",
            "population2_exonic_background",
            group2.loc[group2["type"].eq("introner_containing"), "rate"],
            group2.loc[group2["type"].eq("non_introner_containing"), "rate"],
        )
    )
    rows.append(
        test_row(
            "Population 1 polymorphic canonical introns",
            "Population 1",
            "polymorphic_canonical_intron",
            "common_exonic_background",
            group1.loc[
                group1["comparison_class"].eq("polymorphic_canonical_intron"),
                "rate",
            ],
            common_background,
        )
    )

    results = pd.DataFrame(rows)
    n_tests = len(results)
    results["p_value_bonferroni"] = (
        results["p_value_raw"] * n_tests
    ).clip(upper=1.0)
    results["significant_bonferroni_0.05"] = results[
        "p_value_bonferroni"
    ].lt(0.05)
    results["alternative"] = "two-sided"
    results["rate_filter"] = "finite values only; no min/max cutoff"
    results["n_tests_bonferroni"] = n_tests
    results["group2_mating_rows_excluded"] = n_group2_mating_rows

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False)
    print(results.to_string(index=False))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

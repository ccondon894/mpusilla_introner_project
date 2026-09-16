#!/usr/bin/env python3
"""Compare lengths of fixed and polymorphic introners at unique ortholog loci."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats


GROUPS = {
    "group1": (
        "CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482", "RCC373",
        "RCC465", "RCC629", "RCC692", "RCC693", "RCC833",
    ),
    "group2": ("RCC1749", "RCC3052"),
}
SAMPLE_GROUP = {sample: group for group, samples in GROUPS.items() for sample in samples}
RNA_SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
MT_CONTIG = "CCMP1545#0#scaffold_2"
MT_START = 49_808
MT_END = 1_730_591
FIXED_COLOR = "#326b77"
POLYMORPHIC_COLOR = "#80ae9a"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix", type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
    )
    parser.add_argument(
        "--boundary-dir", type=Path,
        default=Path("results/expression/splice_junctions/introner_boundary_support"),
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("results/expression/functional/fixed_vs_polymorphic_introner_lengths"),
    )
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--bootstraps", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=1545)
    return parser.parse_args()


def bh(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(ranked, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def mt_orthologs(matrix: pd.DataFrame) -> set[str]:
    ccmp = matrix[matrix["sample"].eq("CCMP1545")]
    overlap = (
        ccmp["contig"].eq(MT_CONTIG)
        & ccmp["start"].lt(MT_END)
        & ccmp["end"].gt(MT_START)
    )
    return set(ccmp.loc[overlap, "ortholog_id"])


def classify(row: pd.Series, group: str) -> str:
    if row[f"{group}_pattern"] == "fixed_present":
        return "fixed"
    present = int(row[f"{group}_present_count"])
    callable_count = int(row[f"{group}_callable_count"])
    n_samples = int(row[f"{group}_n_samples"])
    if callable_count == n_samples and 0 < present < n_samples:
        return "polymorphic"
    return "excluded"


def parse_splice_span(value: object) -> float:
    if pd.isna(value):
        return math.nan
    text = str(value)
    donor = re.search(r"\b(?:GT|GC)@(\d+)", text)
    acceptors = list(re.finditer(r"\bAG@(\d+)", text))
    if donor is None or not acceptors:
        return math.nan
    return float(int(acceptors[-1].group(1)) + 2 - int(donor.group(1)))


def present_rows(matrix: pd.DataFrame, group: str, samples: tuple[str, ...]) -> pd.DataFrame:
    excluded_mt = mt_orthologs(matrix)
    data = matrix.loc[
        matrix["sample"].isin(samples)
        & matrix["presence"].eq(1)
        & ~matrix["ortholog_id"].isin(excluded_mt)
    ].copy()
    data["locus_class"] = data.apply(classify, axis=1, group=group)
    return data[data["locus_class"].isin(["fixed", "polymorphic"])]


def build_body_lengths(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, samples in GROUPS.items():
        data = present_rows(matrix, group, samples)
        data["length"] = data["end"] - data["start"] - 200
        if data["length"].le(0).any():
            raise ValueError(f"Non-positive introner body length found in {group}.")
        collapsed = (
            data.groupby(["ortholog_id", "locus_class", "family"])
            .agg(
                length=("length", "median"),
                minimum_present_copy_length=("length", "min"),
                maximum_present_copy_length=("length", "max"),
                present_copy_length_sd=("length", "std"),
                n_present_copies=("length", "size"),
            )
            .reset_index()
        )
        collapsed.insert(0, "group", group)
        collapsed["source"] = "body_length"
        rows.append(collapsed)
    return pd.concat(rows, ignore_index=True)


def build_annotated_spans(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, samples in GROUPS.items():
        data = present_rows(matrix, group, samples)
        data["body_length"] = data["end"] - data["start"] - 200
        data["length"] = data["splice_site"].map(parse_splice_span)
        data = data[
            data["length"].notna()
            & data["length"].gt(0)
            & data["length"].le(data["body_length"])
        ]
        collapsed = (
            data.groupby(["ortholog_id", "locus_class", "family"])
            .agg(
                length=("length", "median"),
                minimum_present_copy_length=("length", "min"),
                maximum_present_copy_length=("length", "max"),
                present_copy_length_sd=("length", "std"),
                n_present_copies=("length", "size"),
            )
            .reset_index()
        )
        collapsed.insert(0, "group", group)
        collapsed["source"] = "annotated_splice_span"
        rows.append(collapsed)
    return pd.concat(rows, ignore_index=True)


def build_rna_junction_lengths(
    matrix: pd.DataFrame, boundary_dir: Path
) -> pd.DataFrame:
    excluded_mt = mt_orthologs(matrix)
    join_columns = [
        "ortholog_id", "family",
        "group1_pattern", "group1_present_count", "group1_callable_count", "group1_n_samples",
        "group2_pattern", "group2_present_count", "group2_callable_count", "group2_n_samples",
    ]
    rows = []
    for sample in RNA_SAMPLES:
        boundary = pd.read_csv(
            boundary_dir / f"{sample}.per_locus.tsv", sep="\t", low_memory=False
        )
        boundary = boundary.loc[
            boundary["best_junction_strand"].isin(["+", "-"])
            & boundary["max_abs_boundary_delta"].le(2)
            & ~boundary["ortholog_id"].isin(excluded_mt)
        ].copy()
        genotype = matrix.loc[
            matrix["sample"].eq(sample), join_columns
        ].drop_duplicates("ortholog_id")
        boundary = boundary.merge(
            genotype, on="ortholog_id", how="left", suffixes=("", "_gt")
        )
        group = SAMPLE_GROUP[sample]
        boundary["locus_class"] = boundary.apply(classify, axis=1, group=group)
        boundary = boundary[boundary["locus_class"].isin(["fixed", "polymorphic"])]
        boundary["length"] = (
            boundary["best_junction_end"] - boundary["best_junction_start"]
        )
        boundary["group"] = group
        boundary["source"] = "rna_supported_junction_length"
        rows.append(
            boundary[[
                "group", "sample", "ortholog_id", "locus_class", "family", "length",
                "best_junction_score", "max_abs_boundary_delta",
                "within_2bp_replicate_count", "source",
            ]]
        )
    occurrences = pd.concat(rows, ignore_index=True)
    collapsed = (
        occurrences.sort_values(
            ["group", "ortholog_id", "best_junction_score"],
            ascending=[True, True, False],
        )
        .drop_duplicates(["group", "ortholog_id"], keep="first")
        .reset_index(drop=True)
    )
    collapsed["minimum_present_copy_length"] = collapsed["length"]
    collapsed["maximum_present_copy_length"] = collapsed["length"]
    collapsed["present_copy_length_sd"] = math.nan
    collapsed["n_present_copies"] = 1
    return collapsed


def families_with_both_classes(data: pd.DataFrame) -> list[int]:
    table = pd.crosstab(data["family"], data["locus_class"])
    for column in ("fixed", "polymorphic"):
        if column not in table:
            table[column] = 0
    return table.index[(table["fixed"] > 0) & (table["polymorphic"] > 0)].tolist()


def stratified_permutation(
    data: pd.DataFrame, rng: np.random.Generator, n_permutations: int
) -> tuple[float, float, float, float]:
    values = data["length"].to_numpy(float)
    labels = data["locus_class"].eq("polymorphic").to_numpy()
    observed = float(values[labels].mean() - values[~labels].mean())
    indices = [index.to_numpy() for _, index in data.groupby("family").groups.items()]
    null = np.empty(n_permutations)
    for iteration in range(n_permutations):
        permuted = labels.copy()
        for index in indices:
            permuted[index] = rng.permutation(permuted[index])
        null[iteration] = values[permuted].mean() - values[~permuted].mean()
    p_two = (np.count_nonzero(np.abs(null) >= abs(observed)) + 1) / (
        n_permutations + 1
    )
    p_shorter = (np.count_nonzero(null <= observed) + 1) / (n_permutations + 1)
    return observed, float(p_two), float(p_shorter), float(null.mean())


def bootstrap_ci(
    fixed: np.ndarray,
    polymorphic: np.ndarray,
    rng: np.random.Generator,
    n_bootstraps: int,
) -> tuple[float, float]:
    differences = np.empty(n_bootstraps)
    for iteration in range(n_bootstraps):
        differences[iteration] = (
            rng.choice(polymorphic, len(polymorphic), replace=True).mean()
            - rng.choice(fixed, len(fixed), replace=True).mean()
        )
    return tuple(np.quantile(differences, [0.025, 0.975]))


def compare_lengths(
    datasets: list[pd.DataFrame],
    rng: np.random.Generator,
    n_permutations: int,
    n_bootstraps: int,
) -> pd.DataFrame:
    rows = []
    for source_data in datasets:
        source = str(source_data["source"].iloc[0])
        for group, data in source_data.groupby("group", sort=False):
            eligible_families = families_with_both_classes(data)
            tested = data[data["family"].isin(eligible_families)].copy().reset_index(drop=True)
            fixed = tested.loc[tested["locus_class"].eq("fixed"), "length"].to_numpy(float)
            poly = tested.loc[
                tested["locus_class"].eq("polymorphic"), "length"
            ].to_numpy(float)
            observed, permutation_p, p_shorter, null_mean = stratified_permutation(
                tested, rng, n_permutations
            )
            ci_low, ci_high = bootstrap_ci(fixed, poly, rng, n_bootstraps)
            mann_whitney = stats.mannwhitneyu(poly, fixed, alternative="two-sided")
            tested["is_polymorphic"] = tested["locus_class"].eq("polymorphic").astype(int)
            regression = smf.ols(
                "np.log(length) ~ is_polymorphic + C(family)", data=tested
            ).fit(cov_type="HC3")
            coefficient = float(regression.params["is_polymorphic"])
            coefficient_ci = regression.conf_int().loc["is_polymorphic"]
            rows.append(
                {
                    "source": source,
                    "group": group,
                    "families_tested": ";".join(map(str, eligible_families)),
                    "fixed_n": len(fixed),
                    "polymorphic_n": len(poly),
                    "fixed_mean": float(fixed.mean()),
                    "polymorphic_mean": float(poly.mean()),
                    "fixed_median": float(np.median(fixed)),
                    "polymorphic_median": float(np.median(poly)),
                    "polymorphic_minus_fixed_mean_bp": observed,
                    "bootstrap_ci_low_bp": ci_low,
                    "bootstrap_ci_high_bp": ci_high,
                    "polymorphic_minus_fixed_median_bp": float(
                        np.median(poly) - np.median(fixed)
                    ),
                    "family_stratified_permutation_p_two_sided": permutation_p,
                    "family_stratified_permutation_p_polymorphic_shorter": p_shorter,
                    "permutation_null_mean_bp": null_mean,
                    "mannwhitney_u": float(mann_whitney.statistic),
                    "mannwhitney_p_two_sided": float(mann_whitney.pvalue),
                    "family_adjusted_log_length_coefficient": coefficient,
                    "family_adjusted_percent_difference": 100 * math.expm1(coefficient),
                    "family_adjusted_percent_ci_low": 100 * math.expm1(float(coefficient_ci.iloc[0])),
                    "family_adjusted_percent_ci_high": 100 * math.expm1(float(coefficient_ci.iloc[1])),
                    "family_adjusted_regression_p": float(regression.pvalues["is_polymorphic"]),
                    "permutations": n_permutations,
                    "bootstraps": n_bootstraps,
                }
            )
    results = pd.DataFrame(rows)
    results["permutation_fdr_bh"] = bh(
        results["family_stratified_permutation_p_two_sided"]
    )
    results["regression_fdr_bh"] = bh(results["family_adjusted_regression_p"])
    return results


def summarize(data: pd.DataFrame) -> pd.DataFrame:
    return (
        data.groupby(["source", "group", "locus_class"])
        .agg(
            n_loci=("ortholog_id", "size"),
            mean_length=("length", "mean"),
            median_length=("length", "median"),
            sd_length=("length", "std"),
            minimum_length=("length", "min"),
            maximum_length=("length", "max"),
            loci_variable_among_present_copies=(
                "present_copy_length_sd", lambda x: int(x.fillna(0).gt(0).sum())
            ),
        )
        .reset_index()
    )


def family_summary(body: pd.DataFrame) -> pd.DataFrame:
    return (
        body.groupby(["group", "family", "locus_class"])
        .agg(
            n_loci=("ortholog_id", "size"),
            mean_length=("length", "mean"),
            median_length=("length", "median"),
            sd_length=("length", "std"),
        )
        .reset_index()
    )


def plot_lengths(body: pd.DataFrame, output: Path) -> None:
    rng = np.random.default_rng(1545)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, group in zip(axes, GROUPS):
        data = body[body["group"].eq(group)]
        arrays = [
            data.loc[data["locus_class"].eq(label), "length"].to_numpy()
            for label in ("fixed", "polymorphic")
        ]
        violin = ax.violinplot(arrays, positions=[0, 1], widths=0.8, showextrema=False)
        for patch, color in zip(violin["bodies"], [FIXED_COLOR, POLYMORPHIC_COLOR]):
            patch.set_facecolor(color)
            patch.set_edgecolor(color)
            patch.set_alpha(0.45)
        boxes = ax.boxplot(
            arrays, positions=[0, 1], widths=0.25, showfliers=False,
            patch_artist=True, medianprops={"color": "black"},
        )
        for patch, color in zip(boxes["boxes"], [FIXED_COLOR, POLYMORPHIC_COLOR]):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)
        for position, values, color in zip(
            [0, 1], arrays, [FIXED_COLOR, POLYMORPHIC_COLOR]
        ):
            sample = rng.choice(values, min(750, len(values)), replace=False)
            ax.scatter(
                rng.normal(position, 0.045, len(sample)), sample,
                color=color, s=7, alpha=0.12, linewidth=0,
            )
        ax.set_title(group.replace("group", "Group "))
        ax.set_xticks(
            [0, 1],
            [f"Fixed\n(n={len(arrays[0])})", f"Polymorphic\n(n={len(arrays[1])})"],
        )
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Median introner body length across present copies (bp)")
    fig.suptitle("Fixed and polymorphic introner lengths")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_results(path: Path, tests: pd.DataFrame) -> None:
    body = tests[tests["source"].eq("body_length")]
    lines = [
        "# Fixed versus polymorphic introner lengths",
        "",
        "The primary analysis uses the median body length across present copies for",
        "each independent ortholog locus. Body length excludes the 100-bp matrix flank",
        "on each side. Mating-type-region orthologs are excluded, populations are tested",
        "separately, and label permutations are restricted within introner family.",
        "",
        "## Primary body-length results",
        "",
        "| Group | Fixed n | Poly n | Fixed median | Poly median | Poly - fixed mean | 95% bootstrap CI | Family-stratified P | BH FDR | Adjusted % difference | Regression P |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in body.itertuples():
        lines.append(
            f"| {row.group} | {row.fixed_n} | {row.polymorphic_n} | "
            f"{row.fixed_median:.1f} bp | {row.polymorphic_median:.1f} bp | "
            f"{row.polymorphic_minus_fixed_mean_bp:.2f} bp | "
            f"{row.bootstrap_ci_low_bp:.2f} to {row.bootstrap_ci_high_bp:.2f} bp | "
            f"{row.family_stratified_permutation_p_two_sided:.4g} | "
            f"{row.permutation_fdr_bh:.4g} | "
            f"{row.family_adjusted_percent_difference:.2f}% | "
            f"{row.family_adjusted_regression_p:.4g} |"
        )
    lines.extend([
        "", "## Boundary-definition sensitivity analyses", "",
        "| Source | Group | Fixed n | Poly n | Fixed median | Poly median | Family-stratified P | Adjusted % difference | Regression P |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in tests[~tests["source"].eq("body_length")].itertuples():
        lines.append(
            f"| {row.source} | {row.group} | {row.fixed_n} | {row.polymorphic_n} | "
            f"{row.fixed_median:.1f} bp | {row.polymorphic_median:.1f} bp | "
            f"{row.family_stratified_permutation_p_two_sided:.4g} | "
            f"{row.family_adjusted_percent_difference:.2f}% | "
            f"{row.family_adjusted_regression_p:.4g} |"
        )
    lines.extend([
        "", "## Notes", "",
        "- `fixed` means fixed-present in the focal population; `polymorphic` requires complete population callability and both present and absent samples.",
        "- Families occurring in only one class are summarized but omitted from family-stratified tests.",
        "- Annotated splice span is calculated from the GT/GC donor offset through the AG acceptor offset.",
        "- RNA-supported junction length uses the strongest stranded junction within 2 bp of the expected boundary for each population-level ortholog.",
        "- Group 2 RNA-supported polymorphic loci are limited to variants present in RCC1749 because RCC3052 lacks RNA-seq data.",
    ])
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    if args.permutations < 1 or args.bootstraps < 1:
        raise ValueError("Permutation and bootstrap counts must be positive.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    body = build_body_lengths(matrix)
    annotated = build_annotated_spans(matrix)
    rna = build_rna_junction_lengths(matrix, args.boundary_dir)
    rng = np.random.default_rng(args.seed)
    tests = compare_lengths(
        [body, annotated, rna], rng, args.permutations, args.bootstraps
    )

    body.to_csv(args.outdir / "locus_body_lengths.tsv", sep="\t", index=False)
    annotated.to_csv(args.outdir / "annotated_splice_spans.tsv", sep="\t", index=False)
    rna.to_csv(args.outdir / "rna_supported_junction_lengths.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "length_tests.tsv", sep="\t", index=False)
    pd.concat([body, annotated, rna], ignore_index=True).pipe(summarize).to_csv(
        args.outdir / "length_summary.tsv", sep="\t", index=False
    )
    family_summary(body).to_csv(
        args.outdir / "body_length_family_summary.tsv", sep="\t", index=False
    )
    plot_lengths(body, args.outdir / "introner_length_comparison.png")
    write_results(args.outdir / "RESULTS.md", tests)
    print(f"Body-length loci: {len(body)}")
    print(f"Annotated-splice loci: {len(annotated)}")
    print(f"RNA-supported loci: {len(rna)}")
    print(f"Wrote results to {args.outdir}")


if __name__ == "__main__":
    main()

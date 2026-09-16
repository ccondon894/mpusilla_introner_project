#!/usr/bin/env python3
"""Compare splicing efficiency for polymorphic and fixed-present introners."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
SAMPLE_TO_GROUP = {"CCMP1545": "group1", "RCC1614": "group1", "RCC1749": "group2"}
GROUP_N = {"group1": 11, "group2": 2}
MATING_CONTIG = "CCMP1545#0#scaffold_2"
MATING_START = 49808
MATING_END = 1730591
PRIMARY_CLASSES = ("fixed_present", "polymorphic")
BINARY_ENDPOINTS = (
    "within_2bp_supported",
    "low_retention_psi_flag",
    "high_confidence_efficient_splicing_flag",
)
CONTINUOUS_ENDPOINTS = ("retention_psi", "spliced_fraction")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
        help="Final introner genotype matrix.",
    )
    parser.add_argument(
        "--psi-dir",
        type=Path,
        default=Path("results/expression/splice_junctions/introner_boundary_support"),
        help="Directory containing SAMPLE.per_locus.tsv PSI tables.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/expression/functional/introner_splicing_efficiency"),
        help="Output directory.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=10000,
        help="Number of stratified label permutations per endpoint.",
    )
    parser.add_argument(
        "--polymorphic-min-total-present",
        type=int,
        default=0,
        help="Minimum total strains with presence for polymorphic loci; 0 disables.",
    )
    parser.add_argument(
        "--polymorphic-max-total-present",
        type=int,
        default=0,
        help="Maximum total strains with presence for polymorphic loci; 0 disables.",
    )
    parser.add_argument("--seed", type=int, default=1545)
    return parser.parse_args()


def load_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t")
    numeric_cols = [
        "start",
        "end",
        "presence",
        "group1_present_count",
        "group1_callable_count",
        "group1_n_samples",
        "group2_present_count",
        "group2_callable_count",
        "group2_n_samples",
    ]
    for col in numeric_cols:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")
    total_present = (
        matrix["presence"].eq(1).groupby(matrix["ortholog_id"]).sum().astype(int)
    )
    matrix["total_present_count"] = matrix["ortholog_id"].map(total_present).astype(int)
    return matrix


def mating_type_orthologs(matrix: pd.DataFrame) -> set[str]:
    ccmp = matrix[matrix["sample"] == "CCMP1545"].copy()
    mask = (
        (ccmp["contig"] == MATING_CONTIG)
        & (ccmp["start"] < MATING_END)
        & (MATING_START < ccmp["end"])
    )
    return set(ccmp.loc[mask, "ortholog_id"])


def total_present_in_range(row: pd.Series, args: argparse.Namespace) -> bool:
    total_present = int(row["total_present_count"])
    if args.polymorphic_min_total_present and total_present < args.polymorphic_min_total_present:
        return False
    if args.polymorphic_max_total_present and total_present > args.polymorphic_max_total_present:
        return False
    return True


def classify_genotype(row: pd.Series, group: str, args: argparse.Namespace) -> str:
    pattern = row[f"{group}_pattern"]
    present = int(row[f"{group}_present_count"])
    callable_count = int(row[f"{group}_callable_count"])
    n_samples = int(row[f"{group}_n_samples"])

    if pattern == "fixed_present":
        return "fixed_present"
    if callable_count == n_samples and 0 < present < n_samples:
        if not total_present_in_range(row, args):
            return "excluded_polymorphic_outside_total_present_range"
        return "polymorphic"
    if 0 < present < callable_count:
        if not total_present_in_range(row, args):
            return "excluded_partial_callable_outside_total_present_range"
        return "polymorphic_partial_callable"
    return f"excluded_{pattern or 'unclassified'}"


def load_classified_rows(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = load_matrix(args.matrix)
    mt_orthologs = mating_type_orthologs(matrix)
    genotype_cols = [
        "ortholog_id",
        "sample",
        "presence",
        "within_group_status",
        "cross_group_status",
        "group1_present_count",
        "group1_absent_count",
        "group1_missing_count",
        "group1_callable_count",
        "group1_n_samples",
        "group1_callability",
        "group1_pattern",
        "group2_present_count",
        "group2_absent_count",
        "group2_missing_count",
        "group2_callable_count",
        "group2_n_samples",
        "group2_callability",
        "group2_pattern",
        "total_present_count",
    ]

    rows = []
    audit_rows = []
    for sample in SAMPLES:
        group = SAMPLE_TO_GROUP[sample]
        psi_path = args.psi_dir / f"{sample}.per_locus.tsv"
        psi = pd.read_csv(psi_path, sep="\t")
        geno = matrix.loc[
            (matrix["sample"] == sample) & (matrix["presence"] == 1), genotype_cols
        ].copy()
        merged = psi.merge(geno, on=["sample", "ortholog_id"], how="left", validate="1:1")
        merged["focal_group"] = group
        merged["mating_type_ortholog_excluded"] = merged["ortholog_id"].isin(mt_orthologs)
        merged["analysis_class"] = merged.apply(
            lambda row: classify_genotype(row, group, args)
            if not pd.isna(row["presence"])
            else "excluded_no_genotype_row",
            axis=1,
        )

        audit = (
            merged.groupby(["sample", "focal_group", "analysis_class"], dropna=False)
            .agg(
                n_loci=("ortholog_id", "size"),
                covered_for_psi_loci=("covered_for_psi", lambda x: (x == "yes").sum()),
                mating_type_excluded_loci=(
                    "mating_type_ortholog_excluded",
                    lambda x: int(x.sum()),
                ),
            )
            .reset_index()
        )
        audit_rows.append(audit)

        filtered = merged.loc[
            ~merged["mating_type_ortholog_excluded"]
            & merged["analysis_class"].isin(
                PRIMARY_CLASSES + ("polymorphic_partial_callable",)
            )
        ].copy()
        rows.append(filtered)

    out = pd.concat(rows, ignore_index=True)
    audit = pd.concat(audit_rows, ignore_index=True)
    add_analysis_columns(out)
    return out, audit


def add_analysis_columns(df: pd.DataFrame) -> None:
    numeric_cols = [
        "retention_psi",
        "spliced_fraction",
        "retained_signal",
        "psi_junction_signal_within_2bp",
        "within_2bp_junction_support",
        "body_len",
        "retained_mean_depth",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["covered_for_psi_flag"] = df["covered_for_psi"] == "yes"
    df["within_2bp_supported"] = df["within_2bp_junction_support"] > 0
    df["low_retention_psi_flag"] = df["low_retention_psi"] == "yes"
    df["high_confidence_efficient_splicing_flag"] = (
        df["high_confidence_efficient_splicing"] == "yes"
    )
    df["total_psi_signal"] = (
        df["retained_signal"].fillna(0) + df["psi_junction_signal_within_2bp"].fillna(0)
    )
    df["log1p_total_psi_signal"] = np.log1p(df["total_psi_signal"])
    df["analysis_set"] = np.where(
        df["analysis_class"].isin(PRIMARY_CLASSES), "primary", "sensitivity"
    )


def safe_rate(num: float, den: float) -> float:
    return num / den if den else math.nan


def summarize_by_class(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = df.groupby(["analysis_set", "sample", "focal_group", "analysis_class"])
    for keys, group in grouped:
        analysis_set, sample, focal_group, analysis_class = keys
        covered = group[group["covered_for_psi_flag"]]
        row = {
            "analysis_set": analysis_set,
            "sample": sample,
            "focal_group": focal_group,
            "analysis_class": analysis_class,
            "n_loci": len(group),
            "covered_for_psi_loci": len(covered),
            "median_retention_psi": covered["retention_psi"].median(),
            "mean_retention_psi": covered["retention_psi"].mean(),
            "median_spliced_fraction": covered["spliced_fraction"].median(),
            "mean_spliced_fraction": covered["spliced_fraction"].mean(),
            "within_2bp_supported_loci": int(group["within_2bp_supported"].sum()),
            "within_2bp_supported_fraction": safe_rate(
                int(group["within_2bp_supported"].sum()), len(group)
            ),
            "low_retention_psi_loci": int(group["low_retention_psi_flag"].sum()),
            "low_retention_psi_fraction_of_covered": safe_rate(
                int(group["low_retention_psi_flag"].sum()), len(covered)
            ),
            "high_confidence_efficient_splicing_loci": int(
                group["high_confidence_efficient_splicing_flag"].sum()
            ),
            "high_confidence_efficient_splicing_fraction": safe_rate(
                int(group["high_confidence_efficient_splicing_flag"].sum()), len(group)
            ),
            "median_log1p_total_psi_signal": group["log1p_total_psi_signal"].median(),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def mean_difference(df: pd.DataFrame, endpoint: str) -> float:
    fixed = df.loc[df["analysis_class"] == "fixed_present", endpoint].dropna()
    poly = df.loc[df["analysis_class"] == "polymorphic", endpoint].dropna()
    return poly.mean() - fixed.mean()


def stratified_permutation(
    df: pd.DataFrame, endpoint: str, n_perm: int, rng: np.random.Generator
) -> dict[str, float | int | str]:
    data = df[df["analysis_class"].isin(PRIMARY_CLASSES)].copy().reset_index(drop=True)
    if endpoint in CONTINUOUS_ENDPOINTS:
        data = data[data["covered_for_psi_flag"] & data[endpoint].notna()]
    else:
        data = data[data[endpoint].notna()]
    data = data.reset_index(drop=True)

    observed = mean_difference(data, endpoint)
    if data["analysis_class"].nunique() < 2 or math.isnan(observed):
        return {
            "observed_polymorphic_minus_fixed": math.nan,
            "permutation_p_two_sided": math.nan,
            "permutations": 0,
            "null_mean": math.nan,
            "null_sd": math.nan,
        }

    strata = [
        idx.to_numpy()
        for _, idx in data.groupby(["sample", "family"], dropna=False).groups.items()
    ]
    labels = data["analysis_class"].to_numpy().copy()
    values = data[endpoint].to_numpy()
    valid_perms = []
    for _ in range(n_perm):
        perm_labels = labels.copy()
        for idx in strata:
            perm_labels[idx] = rng.permutation(perm_labels[idx])
        fixed_values = values[perm_labels == "fixed_present"]
        poly_values = values[perm_labels == "polymorphic"]
        if len(fixed_values) and len(poly_values):
            valid_perms.append(float(np.nanmean(poly_values) - np.nanmean(fixed_values)))

    null = np.asarray(valid_perms, dtype=float)
    exceed = np.sum(np.abs(null) >= abs(observed))
    p_value = (exceed + 1) / (len(null) + 1)
    return {
        "observed_polymorphic_minus_fixed": observed,
        "permutation_p_two_sided": p_value,
        "permutations": len(null),
        "null_mean": float(np.mean(null)) if len(null) else math.nan,
        "null_sd": float(np.std(null, ddof=1)) if len(null) > 1 else math.nan,
    }


def fisher_test(df: pd.DataFrame, endpoint: str) -> dict[str, float | int]:
    fixed = df[df["analysis_class"] == "fixed_present"]
    poly = df[df["analysis_class"] == "polymorphic"]
    fixed_yes = int(fixed[endpoint].sum())
    poly_yes = int(poly[endpoint].sum())
    table = [[poly_yes, len(poly) - poly_yes], [fixed_yes, len(fixed) - fixed_yes]]
    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
    return {
        "polymorphic_yes": poly_yes,
        "polymorphic_total": len(poly),
        "fixed_yes": fixed_yes,
        "fixed_total": len(fixed),
        "fisher_odds_ratio": odds_ratio,
        "fisher_p_two_sided": p_value,
    }


def continuous_test(df: pd.DataFrame, endpoint: str) -> dict[str, float | int]:
    covered = df[df["covered_for_psi_flag"]]
    fixed = covered.loc[covered["analysis_class"] == "fixed_present", endpoint].dropna()
    poly = covered.loc[covered["analysis_class"] == "polymorphic", endpoint].dropna()
    if len(fixed) == 0 or len(poly) == 0:
        return {
            "polymorphic_n": len(poly),
            "fixed_n": len(fixed),
            "mannwhitney_u": math.nan,
            "mannwhitney_p_two_sided": math.nan,
        }
    result = stats.mannwhitneyu(poly, fixed, alternative="two-sided")
    return {
        "polymorphic_n": len(poly),
        "fixed_n": len(fixed),
        "mannwhitney_u": result.statistic,
        "mannwhitney_p_two_sided": result.pvalue,
    }


def bootstrap_ci_difference(
    df: pd.DataFrame,
    endpoint: str,
    rng: np.random.Generator,
    n_boot: int = 5000,
) -> tuple[float, float]:
    data = df[df["analysis_class"].isin(PRIMARY_CLASSES)].copy()
    if endpoint in CONTINUOUS_ENDPOINTS:
        data = data[data["covered_for_psi_flag"] & data[endpoint].notna()]
    fixed = data.loc[data["analysis_class"] == "fixed_present", endpoint].to_numpy()
    poly = data.loc[data["analysis_class"] == "polymorphic", endpoint].to_numpy()
    if len(fixed) == 0 or len(poly) == 0:
        return math.nan, math.nan
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        diffs[i] = (
            rng.choice(poly, size=len(poly), replace=True).mean()
            - rng.choice(fixed, size=len(fixed), replace=True).mean()
        )
    return tuple(np.quantile(diffs, [0.025, 0.975]))


def statistical_tests(
    df: pd.DataFrame, n_perm: int, rng: np.random.Generator
) -> pd.DataFrame:
    primary = df[df["analysis_set"] == "primary"].copy()
    rows = []
    analysis_groups = [("all_samples", primary)]
    analysis_groups.extend((sample, sample_df) for sample, sample_df in primary.groupby("sample"))

    for label, subset in analysis_groups:
        for endpoint in BINARY_ENDPOINTS:
            perm = stratified_permutation(subset, endpoint, n_perm, rng)
            fisher = fisher_test(subset, endpoint)
            ci_low, ci_high = bootstrap_ci_difference(subset, endpoint, rng)
            rows.append(
                {
                    "comparison_scope": label,
                    "endpoint": endpoint,
                    "endpoint_type": "binary",
                    **perm,
                    **fisher,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                }
            )
        for endpoint in CONTINUOUS_ENDPOINTS:
            perm = stratified_permutation(subset, endpoint, n_perm, rng)
            mw = continuous_test(subset, endpoint)
            ci_low, ci_high = bootstrap_ci_difference(subset, endpoint, rng)
            rows.append(
                {
                    "comparison_scope": label,
                    "endpoint": endpoint,
                    "endpoint_type": "continuous_secondary",
                    **perm,
                    **mw,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                }
            )
    return pd.DataFrame(rows)


def save_bar_plot(summary: pd.DataFrame, outdir: Path) -> None:
    primary = summary[summary["analysis_set"] == "primary"].copy()
    metrics = [
        ("within_2bp_supported_fraction", "Within-2bp junction support"),
        (
            "high_confidence_efficient_splicing_fraction",
            "High-confidence efficient splicing",
        ),
        ("low_retention_psi_fraction_of_covered", "Low retention PSI"),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(12, 4), sharey=False)
    colors = {"fixed_present": "#4C78A8", "polymorphic": "#F58518"}
    x = np.arange(len(SAMPLES))
    width = 0.36
    for ax, (metric, title) in zip(axes, metrics):
        for offset, analysis_class in [(-width / 2, "fixed_present"), (width / 2, "polymorphic")]:
            values = []
            for sample in SAMPLES:
                row = primary[
                    (primary["sample"] == sample)
                    & (primary["analysis_class"] == analysis_class)
                ]
                values.append(float(row[metric].iloc[0]) if len(row) else np.nan)
            ax.bar(
                x + offset,
                values,
                width=width,
                label=analysis_class.replace("_", " "),
                color=colors[analysis_class],
            )
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES, rotation=30, ha="right")
        ax.set_ylim(0, max(0.15, np.nanmax(primary[metric]) * 1.25))
        ax.set_ylabel("Fraction")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "splicing_efficiency_binary_fractions.png", dpi=200)
    plt.close(fig)


def save_continuous_plot(df: pd.DataFrame, outdir: Path) -> None:
    primary = df[
        (df["analysis_set"] == "primary") & df["covered_for_psi_flag"]
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=False)
    colors = {"fixed_present": "#4C78A8", "polymorphic": "#F58518"}
    rng = np.random.default_rng(1545)
    for ax, endpoint in zip(axes, CONTINUOUS_ENDPOINTS):
        labels = []
        positions = []
        values_by_pos = []
        colors_by_pos = []
        pos = 0
        for sample in SAMPLES:
            for analysis_class in PRIMARY_CLASSES:
                values = primary.loc[
                    (primary["sample"] == sample)
                    & (primary["analysis_class"] == analysis_class),
                    endpoint,
                ].dropna()
                labels.append(f"{sample}\n{analysis_class.replace('_', ' ')}")
                positions.append(pos)
                values_by_pos.append(values.to_numpy())
                colors_by_pos.append(colors[analysis_class])
                jitter = rng.normal(loc=pos, scale=0.035, size=min(len(values), 500))
                plot_values = values.sample(
                    n=min(len(values), 500), random_state=1545
                ).to_numpy()
                ax.scatter(
                    jitter[: len(plot_values)],
                    plot_values,
                    s=7,
                    alpha=0.18,
                    color=colors[analysis_class],
                    linewidths=0,
                )
                pos += 1
            pos += 0.6
        box = ax.boxplot(
            values_by_pos,
            positions=positions,
            widths=0.42,
            patch_artist=True,
            showfliers=False,
        )
        for patch, color in zip(box["boxes"], colors_by_pos):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
        ax.axhline(1.0, color="0.35", linestyle=":", linewidth=1)
        ax.set_title(endpoint.replace("_", " ").title())
        ax.set_ylabel(endpoint.replace("_", " "))
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "splicing_efficiency_psi_distributions.png", dpi=200)
    plt.close(fig)


def save_family_effect_plot(df: pd.DataFrame, outdir: Path) -> None:
    primary = df[df["analysis_set"] == "primary"].copy()
    rows = []
    for (sample, family), subset in primary.groupby(["sample", "family"], dropna=False):
        if set(subset["analysis_class"]) != set(PRIMARY_CLASSES):
            continue
        diff = mean_difference(subset, "within_2bp_supported")
        rows.append(
            {
                "sample": sample,
                "family": str(family),
                "diff": diff,
                "n": len(subset),
            }
        )
    effects = pd.DataFrame(rows)
    if effects.empty:
        return
    effects = effects.sort_values(["sample", "family"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    labels = [f"{r.sample} F{r.family}" for r in effects.itertuples()]
    x = np.arange(len(effects))
    ax.bar(x, effects["diff"], color="#54A24B")
    ax.axhline(0, color="0.25", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Polymorphic - fixed fraction")
    ax.set_title("Within-2bp support effect by sample and family")
    fig.tight_layout()
    fig.savefig(outdir / "family_stratified_within2bp_effects.png", dpi=200)
    plt.close(fig)


def validate_outputs(df: pd.DataFrame, matrix: pd.DataFrame) -> None:
    primary = df[df["analysis_set"] == "primary"]
    if set(primary["analysis_class"]) != set(PRIMARY_CLASSES):
        raise ValueError("Primary output must contain fixed_present and polymorphic rows.")
    if primary["mating_type_ortholog_excluded"].any():
        raise ValueError("Mating-type ortholog rows remain after filtering.")
    for sample, group in SAMPLE_TO_GROUP.items():
        sample_rows = primary[primary["sample"] == sample]
        poly = sample_rows[sample_rows["analysis_class"] == "polymorphic"]
        expected_n = GROUP_N[group]
        if not (
            (poly[f"{group}_callable_count"] == expected_n)
            & (poly[f"{group}_present_count"] > 0)
            & (poly[f"{group}_present_count"] < expected_n)
        ).all():
            raise ValueError(f"Polymorphic class validation failed for {sample}.")
        fixed = sample_rows[sample_rows["analysis_class"] == "fixed_present"]
        if not (fixed[f"{group}_pattern"] == "fixed_present").all():
            raise ValueError(f"Fixed-present class validation failed for {sample}.")


def main() -> None:
    args = parse_args()
    if args.permutations < 1:
        raise ValueError("--permutations must be at least 1.")
    args.outdir.mkdir(parents=True, exist_ok=True)

    classified, audit = load_classified_rows(args)
    matrix = load_matrix(args.matrix)
    validate_outputs(classified, matrix)

    summary = summarize_by_class(classified)
    rng = np.random.default_rng(args.seed)
    tests = statistical_tests(classified, args.permutations, rng)

    classified.to_csv(args.outdir / "per_locus_classified.tsv", sep="\t", index=False)
    audit.to_csv(args.outdir / "class_assignment_audit.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "sample_class_summary.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "statistical_tests.tsv", sep="\t", index=False)

    save_bar_plot(summary, args.outdir)
    save_continuous_plot(classified, args.outdir)
    save_family_effect_plot(classified, args.outdir)

    print(f"Wrote classified loci: {args.outdir / 'per_locus_classified.tsv'}")
    print(f"Wrote class summary: {args.outdir / 'sample_class_summary.tsv'}")
    print(f"Wrote statistical tests: {args.outdir / 'statistical_tests.tsv'}")


if __name__ == "__main__":
    main()

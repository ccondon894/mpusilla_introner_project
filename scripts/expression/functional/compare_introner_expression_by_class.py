#!/usr/bin/env python3
"""Compare expression of genes carrying polymorphic versus fixed introners."""

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
PRIMARY_GENE_CLASSES = ("fixed_only", "polymorphic")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("results/expression/isoform_analysis/isoform_introner_data_filtered.csv"),
        help="Gene-strain expression table with mean_expression and log_expression.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/expression/functional/introner_expression_by_class"),
    )
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--bootstraps", type=int, default=5000)
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
    matrix["gene"] = matrix["gene"].fillna("")
    total_present = (
        matrix["presence"].eq(1).groupby(matrix["ortholog_id"]).sum().astype(int)
    )
    matrix["total_present_count"] = matrix["ortholog_id"].map(total_present).astype(int)
    return matrix


def mating_type_orthologs(matrix: pd.DataFrame) -> set[str]:
    ccmp = matrix[matrix["sample"] == "CCMP1545"]
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


def classify_locus(row: pd.Series, group: str, args: argparse.Namespace) -> str:
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


def build_locus_table(
    matrix: pd.DataFrame, args: argparse.Namespace
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mt_orthologs = mating_type_orthologs(matrix)
    rows = []
    audit_rows = []
    for sample in SAMPLES:
        group = SAMPLE_TO_GROUP[sample]
        sample_rows = matrix[
            (matrix["sample"] == sample)
            & (matrix["presence"] == 1)
            & (matrix["gene"] != "")
        ].copy()
        sample_rows["focal_group"] = group
        sample_rows["mating_type_ortholog_excluded"] = sample_rows["ortholog_id"].isin(
            mt_orthologs
        )
        sample_rows["locus_class"] = sample_rows.apply(
            lambda row: classify_locus(row, group, args), axis=1
        )
        audit_rows.append(
            sample_rows.groupby(["sample", "focal_group", "locus_class"], dropna=False)
            .agg(
                n_loci=("ortholog_id", "size"),
                n_genes=("gene", "nunique"),
                mating_type_excluded_loci=(
                    "mating_type_ortholog_excluded",
                    lambda x: int(x.sum()),
                ),
            )
            .reset_index()
        )
        rows.append(sample_rows[~sample_rows["mating_type_ortholog_excluded"]])

    loci = pd.concat(rows, ignore_index=True)
    audit = pd.concat(audit_rows, ignore_index=True)
    return loci, audit


def classify_gene_locus_classes(classes: pd.Series) -> str:
    values = set(classes)
    if "polymorphic" in values:
        return "polymorphic"
    if "excluded_polymorphic_outside_total_present_range" in values:
        return "excluded_common_polymorphic"
    if "fixed_present" in values:
        return "fixed_only"
    if "excluded_partial_callable_outside_total_present_range" in values:
        return "excluded_common_partial_callable"
    if "polymorphic_partial_callable" in values:
        return "polymorphic_partial_callable"
    return "excluded_other"


def build_gene_table(loci: pd.DataFrame, expression_path: Path) -> pd.DataFrame:
    gene_rows = (
        loci.groupby(["sample", "focal_group", "gene"], dropna=False)
        .agg(
            gene_class=("locus_class", classify_gene_locus_classes),
            n_fixed_present_loci=("locus_class", lambda x: int((x == "fixed_present").sum())),
            n_polymorphic_loci=("locus_class", lambda x: int((x == "polymorphic").sum())),
            n_partial_callable_loci=(
                "locus_class",
                lambda x: int((x == "polymorphic_partial_callable").sum()),
            ),
            families=("family", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)}))),
            n_loci=("ortholog_id", "size"),
        )
        .reset_index()
        .rename(columns={"sample": "strain", "gene": "gene_id"})
    )

    expr = pd.read_csv(expression_path)
    expr_cols = [
        "gene_id",
        "strain",
        "n_isoforms",
        "cds_length",
        "mean_expression",
        "log_expression",
        "log_cds_length",
    ]
    gene_rows = gene_rows.merge(expr[expr_cols], on=["gene_id", "strain"], how="left")
    gene_rows["has_expression"] = gene_rows["mean_expression"].notna()
    gene_rows["log2_mean_cpm_plus1"] = np.log2(gene_rows["mean_expression"] + 1)

    expr_quantiles = []
    for strain, strain_expr in expr.groupby("strain"):
        expressed = strain_expr[["gene_id", "mean_expression"]].copy()
        expressed["expression_percentile"] = expressed["mean_expression"].rank(
            pct=True, method="average"
        )
        expressed["bottom_decile"] = expressed["expression_percentile"] <= 0.10
        expressed["bottom_quartile"] = expressed["expression_percentile"] <= 0.25
        expressed["strain"] = strain
        expr_quantiles.append(
            expressed[["gene_id", "strain", "expression_percentile", "bottom_decile", "bottom_quartile"]]
        )
    quantiles = pd.concat(expr_quantiles, ignore_index=True)
    return gene_rows.merge(quantiles, on=["gene_id", "strain"], how="left")


def build_locus_expression_table(loci: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    expr_cols = [
        "gene_id",
        "strain",
        "mean_expression",
        "log_expression",
        "log2_mean_cpm_plus1",
        "expression_percentile",
        "bottom_decile",
        "bottom_quartile",
        "has_expression",
    ]
    out = loci.rename(columns={"sample": "strain", "gene": "gene_id"}).merge(
        genes[expr_cols].drop_duplicates(["gene_id", "strain"]),
        on=["gene_id", "strain"],
        how="left",
    )
    out["analysis_set"] = np.where(
        out["locus_class"].isin(("fixed_present", "polymorphic")),
        "primary_locus_sensitivity",
        "partial_callable_sensitivity",
    )
    return out


def summarize_gene_classes(genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in genes.groupby(["strain", "focal_group", "gene_class"], dropna=False):
        strain, focal_group, gene_class = keys
        expressed = group[group["has_expression"]]
        rows.append(
            {
                "strain": strain,
                "focal_group": focal_group,
                "gene_class": gene_class,
                "n_genes": len(group),
                "genes_with_expression": len(expressed),
                "median_mean_cpm": expressed["mean_expression"].median(),
                "mean_mean_cpm": expressed["mean_expression"].mean(),
                "median_log2_cpm_plus1": expressed["log2_mean_cpm_plus1"].median(),
                "mean_log2_cpm_plus1": expressed["log2_mean_cpm_plus1"].mean(),
                "median_expression_percentile": expressed["expression_percentile"].median(),
                "bottom_decile_genes": int(expressed["bottom_decile"].sum()),
                "bottom_decile_fraction": safe_rate(
                    int(expressed["bottom_decile"].sum()), len(expressed)
                ),
                "bottom_quartile_genes": int(expressed["bottom_quartile"].sum()),
                "bottom_quartile_fraction": safe_rate(
                    int(expressed["bottom_quartile"].sum()), len(expressed)
                ),
                "median_cds_length": expressed["cds_length"].median(),
                "median_n_isoforms": expressed["n_isoforms"].median(),
            }
        )
    return pd.DataFrame(rows)


def safe_rate(num: int, den: int) -> float:
    return num / den if den else math.nan


def mean_diff(data: pd.DataFrame, column: str) -> float:
    fixed = data.loc[data["gene_class"] == "fixed_only", column].dropna()
    poly = data.loc[data["gene_class"] == "polymorphic", column].dropna()
    return float(poly.mean() - fixed.mean())


def bootstrap_ci(
    data: pd.DataFrame,
    column: str,
    rng: np.random.Generator,
    n_boot: int,
) -> tuple[float, float]:
    fixed = data.loc[data["gene_class"] == "fixed_only", column].dropna().to_numpy()
    poly = data.loc[data["gene_class"] == "polymorphic", column].dropna().to_numpy()
    if len(fixed) == 0 or len(poly) == 0:
        return math.nan, math.nan
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        diffs[i] = (
            rng.choice(poly, size=len(poly), replace=True).mean()
            - rng.choice(fixed, size=len(fixed), replace=True).mean()
        )
    return tuple(np.quantile(diffs, [0.025, 0.975]))


def stratified_permutation(
    data: pd.DataFrame,
    column: str,
    rng: np.random.Generator,
    n_perm: int,
) -> dict[str, float | int]:
    data = data[
        data["gene_class"].isin(PRIMARY_GENE_CLASSES) & data[column].notna()
    ].copy()
    data = data.reset_index(drop=True)
    observed = mean_diff(data, column)
    if data["gene_class"].nunique() < 2 or math.isnan(observed):
        return {
            "observed_polymorphic_minus_fixed": math.nan,
            "permutation_p_two_sided": math.nan,
            "permutation_p_polymorphic_lower": math.nan,
            "permutation_p_polymorphic_higher": math.nan,
            "permutations": 0,
            "null_mean": math.nan,
            "null_sd": math.nan,
        }

    labels = data["gene_class"].to_numpy().copy()
    values = data[column].to_numpy()
    strata = [idx.to_numpy() for _, idx in data.groupby("strain").groups.items()]
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm_labels = labels.copy()
        for idx in strata:
            perm_labels[idx] = rng.permutation(perm_labels[idx])
        fixed = values[perm_labels == "fixed_only"]
        poly = values[perm_labels == "polymorphic"]
        null[i] = poly.mean() - fixed.mean()

    p_two = (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    p_lower = (np.sum(null <= observed) + 1) / (len(null) + 1)
    p_higher = (np.sum(null >= observed) + 1) / (len(null) + 1)
    return {
        "observed_polymorphic_minus_fixed": observed,
        "permutation_p_two_sided": p_two,
        "permutation_p_polymorphic_lower": p_lower,
        "permutation_p_polymorphic_higher": p_higher,
        "permutations": len(null),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
    }


def mann_whitney(data: pd.DataFrame, column: str) -> dict[str, float | int]:
    fixed = data.loc[data["gene_class"] == "fixed_only", column].dropna()
    poly = data.loc[data["gene_class"] == "polymorphic", column].dropna()
    if len(fixed) == 0 or len(poly) == 0:
        return {
            "fixed_n": len(fixed),
            "polymorphic_n": len(poly),
            "mannwhitney_u": math.nan,
            "mannwhitney_p_two_sided": math.nan,
            "mannwhitney_p_polymorphic_lower": math.nan,
        }
    two = stats.mannwhitneyu(poly, fixed, alternative="two-sided")
    lower = stats.mannwhitneyu(poly, fixed, alternative="less")
    return {
        "fixed_n": len(fixed),
        "polymorphic_n": len(poly),
        "mannwhitney_u": float(two.statistic),
        "mannwhitney_p_two_sided": float(two.pvalue),
        "mannwhitney_p_polymorphic_lower": float(lower.pvalue),
    }


def fisher_binary(data: pd.DataFrame, column: str) -> dict[str, float | int]:
    fixed = data[data["gene_class"] == "fixed_only"]
    poly = data[data["gene_class"] == "polymorphic"]
    fixed_yes = int(fixed[column].sum())
    poly_yes = int(poly[column].sum())
    table = [[poly_yes, len(poly) - poly_yes], [fixed_yes, len(fixed) - fixed_yes]]
    odds, p_two = stats.fisher_exact(table, alternative="two-sided")
    _, p_greater = stats.fisher_exact(table, alternative="greater")
    return {
        "fixed_yes": fixed_yes,
        "fixed_total": len(fixed),
        "polymorphic_yes": poly_yes,
        "polymorphic_total": len(poly),
        "fisher_odds_ratio": float(odds),
        "fisher_p_two_sided": float(p_two),
        "fisher_p_polymorphic_enriched": float(p_greater),
    }


def statistical_tests(
    genes: pd.DataFrame,
    rng: np.random.Generator,
    n_perm: int,
    n_boot: int,
) -> pd.DataFrame:
    primary = genes[
        genes["gene_class"].isin(PRIMARY_GENE_CLASSES) & genes["has_expression"]
    ].copy()
    scopes = [("all_samples", primary)]
    scopes.extend((strain, sub) for strain, sub in primary.groupby("strain"))
    rows = []

    for scope, data in scopes:
        for column in ("log2_mean_cpm_plus1", "expression_percentile"):
            perm = stratified_permutation(data, column, rng, n_perm)
            mw = mann_whitney(data, column)
            ci_low, ci_high = bootstrap_ci(data, column, rng, n_boot)
            rows.append(
                {
                    "comparison_scope": scope,
                    "endpoint": column,
                    "endpoint_type": "continuous",
                    **perm,
                    **mw,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                }
            )
        for column in ("bottom_decile", "bottom_quartile"):
            as_float = data.copy()
            as_float[column] = as_float[column].astype(float)
            perm = stratified_permutation(as_float, column, rng, n_perm)
            fisher = fisher_binary(data, column)
            ci_low, ci_high = bootstrap_ci(as_float, column, rng, n_boot)
            rows.append(
                {
                    "comparison_scope": scope,
                    "endpoint": column,
                    "endpoint_type": "binary_low_expression",
                    **perm,
                    **fisher,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                }
            )
    return pd.DataFrame(rows)


def save_distribution_plot(genes: pd.DataFrame, outdir: Path) -> None:
    primary = genes[
        genes["gene_class"].isin(PRIMARY_GENE_CLASSES) & genes["has_expression"]
    ].copy()
    colors = {"fixed_only": "#4C78A8", "polymorphic": "#F58518"}
    rng = np.random.default_rng(1545)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=False)
    for ax, endpoint, ylabel in [
        (axes[0], "log2_mean_cpm_plus1", "log2(mean CPM + 1)"),
        (axes[1], "expression_percentile", "Expression percentile"),
    ]:
        positions = []
        labels = []
        values_by_pos = []
        colors_by_pos = []
        pos = 0
        for sample in SAMPLES:
            for gene_class in PRIMARY_GENE_CLASSES:
                values = primary.loc[
                    (primary["strain"] == sample) & (primary["gene_class"] == gene_class),
                    endpoint,
                ].dropna()
                positions.append(pos)
                labels.append(f"{sample}\n{gene_class.replace('_', ' ')}")
                values_by_pos.append(values.to_numpy())
                colors_by_pos.append(colors[gene_class])
                if len(values):
                    plot_values = values.sample(
                        n=min(len(values), 500), random_state=1545
                    ).to_numpy()
                    jitter = rng.normal(loc=pos, scale=0.04, size=len(plot_values))
                    ax.scatter(
                        jitter,
                        plot_values,
                        s=7,
                        alpha=0.18,
                        linewidths=0,
                        color=colors[gene_class],
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
        ax.set_ylabel(ylabel)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "gene_expression_distributions.png", dpi=200)
    plt.close(fig)


def save_low_expression_plot(summary: pd.DataFrame, outdir: Path) -> None:
    primary = summary[summary["gene_class"].isin(PRIMARY_GENE_CLASSES)].copy()
    colors = {"fixed_only": "#4C78A8", "polymorphic": "#F58518"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    x = np.arange(len(SAMPLES))
    width = 0.36
    for ax, metric, title in [
        (axes[0], "bottom_decile_fraction", "Bottom expression decile"),
        (axes[1], "bottom_quartile_fraction", "Bottom expression quartile"),
    ]:
        for offset, gene_class in [(-width / 2, "fixed_only"), (width / 2, "polymorphic")]:
            values = []
            for sample in SAMPLES:
                row = primary[
                    (primary["strain"] == sample) & (primary["gene_class"] == gene_class)
                ]
                values.append(float(row[metric].iloc[0]) if len(row) else np.nan)
            ax.bar(
                x + offset,
                values,
                width=width,
                color=colors[gene_class],
                label=gene_class.replace("_", " "),
            )
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES, rotation=30, ha="right")
        ax.set_ylabel("Fraction of genes")
        ax.set_ylim(0, max(0.35, np.nanmax(primary[metric]) * 1.25))
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "low_expression_fraction_by_class.png", dpi=200)
    plt.close(fig)


def validate_outputs(genes: pd.DataFrame) -> None:
    primary = genes[genes["gene_class"].isin(PRIMARY_GENE_CLASSES)]
    if primary["mating_type_ortholog_excluded"].any() if "mating_type_ortholog_excluded" in primary.columns else False:
        raise ValueError("Mating-type rows remain in primary genes.")
    for sample, group in SAMPLE_TO_GROUP.items():
        sample_genes = primary[primary["strain"] == sample]
        if not set(PRIMARY_GENE_CLASSES).issubset(set(sample_genes["gene_class"])):
            raise ValueError(f"Missing primary classes for {sample}.")
        expected_n = GROUP_N[group]
        poly = sample_genes[sample_genes["gene_class"] == "polymorphic"]
        if (poly["n_polymorphic_loci"] <= 0).any():
            raise ValueError(f"Polymorphic gene with no polymorphic loci in {sample}.")
        sample_loci = poly["n_fixed_present_loci"].notna()
        if not sample_loci.all():
            raise ValueError(f"Missing locus counts for {sample}.")
        if expected_n not in (2, 11):
            raise ValueError(f"Unexpected group size for {sample}.")


def main() -> None:
    args = parse_args()
    if args.permutations < 1 or args.bootstraps < 1:
        raise ValueError("--permutations and --bootstraps must be at least 1.")
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrix = load_matrix(args.matrix)
    loci, audit = build_locus_table(matrix, args)
    genes = build_gene_table(loci, args.expression)
    locus_expr = build_locus_expression_table(loci, genes)
    validate_outputs(genes)

    summary = summarize_gene_classes(genes)
    rng = np.random.default_rng(args.seed)
    tests = statistical_tests(genes, rng, args.permutations, args.bootstraps)

    genes.to_csv(args.outdir / "gene_class_expression.tsv", sep="\t", index=False)
    locus_expr.to_csv(args.outdir / "locus_class_expression.tsv", sep="\t", index=False)
    audit.to_csv(args.outdir / "class_assignment_audit.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "sample_gene_class_summary.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "statistical_tests.tsv", sep="\t", index=False)

    save_distribution_plot(genes, args.outdir)
    save_low_expression_plot(summary, args.outdir)

    print(f"Wrote gene-level expression table: {args.outdir / 'gene_class_expression.tsv'}")
    print(f"Wrote summary: {args.outdir / 'sample_gene_class_summary.tsv'}")
    print(f"Wrote tests: {args.outdir / 'statistical_tests.tsv'}")


if __name__ == "__main__":
    main()

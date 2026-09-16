#!/usr/bin/env python3
"""Compare fixed-only and polymorphic-introner genes using RNA-seq counts.

This exploratory analysis deliberately starts from the complete featureCounts
matrix rather than from a long-read/isoform-qualified gene table.  Expression
comparisons are performed separately for each strain.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.genmod.cov_struct import Exchangeable


SAMPLE_REPLICATES = {
    "CCMP1545": ("834_A1", "834_A2", "834_A3", "834_B1"),
    "RCC1614": ("1614_A1", "1614_A2", "1614_A3", "1614_A4"),
    "RCC1749": ("1749_A1", "1749_A2", "1749_A4", "1749_B1"),
}
SAMPLE_TO_GROUP = {
    "CCMP1545": "group1",
    "RCC1614": "group1",
    "RCC1749": "group2",
}
MT_CONTIG = "CCMP1545#0#scaffold_2"
MT_START = 49_808
MT_END = 1_730_591
FIXED_COLOR = "#326b77"
POLYMORPHIC_COLOR = "#80ae9a"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=Path("results/expression/counts/merged_counts_matrix.csv"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/expression/functional/introner_expression_by_class_rnaseq_counts"),
    )
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--bootstraps", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=1545)
    return parser.parse_args()


def load_genotypes(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t")
    numeric = [
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
    for column in numeric:
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce")
    matrix["gene"] = matrix["gene"].fillna("")
    return matrix


def mating_type_orthologs(matrix: pd.DataFrame) -> set[str]:
    ccmp = matrix[matrix["sample"].eq("CCMP1545")]
    overlaps = (
        ccmp["contig"].eq(MT_CONTIG)
        & ccmp["start"].lt(MT_END)
        & ccmp["end"].gt(MT_START)
    )
    return set(ccmp.loc[overlaps, "ortholog_id"])


def classify_locus(row: pd.Series, group: str) -> str:
    pattern = str(row[f"{group}_pattern"])
    present = int(row[f"{group}_present_count"])
    callable_count = int(row[f"{group}_callable_count"])
    n_samples = int(row[f"{group}_n_samples"])
    if pattern == "fixed_present":
        return "fixed_present"
    if callable_count == n_samples and 0 < present < n_samples:
        return "polymorphic"
    if 0 < present < callable_count:
        return "polymorphic_partial_callable"
    return f"uncertain_{pattern or 'unclassified'}"


def classify_gene(classes: pd.Series) -> str:
    values = set(classes)
    if "polymorphic" in values:
        return "polymorphic"
    # The control is deliberately strict: every present introner assigned to
    # the gene must be confidently fixed-present in the focal population.
    if values == {"fixed_present"}:
        return "fixed_only"
    return "excluded_uncertain"


def build_gene_classes(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mt_orthologs = mating_type_orthologs(matrix)
    gene_tables = []
    audit_tables = []
    for strain, group in SAMPLE_TO_GROUP.items():
        loci = matrix.loc[
            matrix["sample"].eq(strain)
            & matrix["presence"].eq(1)
            & matrix["gene"].ne("")
            & ~matrix["ortholog_id"].isin(mt_orthologs)
        ].copy()
        loci["locus_class"] = loci.apply(classify_locus, axis=1, group=group)
        audit = (
            loci.groupby("locus_class")
            .agg(n_loci=("ortholog_id", "size"), n_genes=("gene", "nunique"))
            .reset_index()
        )
        audit.insert(0, "strain", strain)
        audit_tables.append(audit)
        genes = (
            loci.groupby("gene")
            .agg(
                gene_class=("locus_class", classify_gene),
                n_fixed_present_loci=(
                    "locus_class",
                    lambda x: int(x.eq("fixed_present").sum()),
                ),
                n_polymorphic_loci=(
                    "locus_class",
                    lambda x: int(x.eq("polymorphic").sum()),
                ),
                n_uncertain_loci=(
                    "locus_class",
                    lambda x: int((~x.isin(["fixed_present", "polymorphic"])).sum()),
                ),
                n_present_introner_loci=("ortholog_id", "size"),
            )
            .reset_index()
            .rename(columns={"gene": "gene_id"})
        )
        genes.insert(0, "strain", strain)
        genes.insert(1, "focal_group", group)
        gene_tables.append(genes)
    return pd.concat(gene_tables, ignore_index=True), pd.concat(audit_tables, ignore_index=True)


def normalize_counts(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"Geneid", "Length"} | {
        replicate for replicates in SAMPLE_REPLICATES.values() for replicate in replicates
    }
    missing = required - set(counts.columns)
    if missing:
        raise ValueError(f"Count matrix is missing columns: {sorted(missing)}")
    if counts["Geneid"].duplicated().any():
        raise ValueError("Count matrix contains duplicate Geneid values.")
    if counts["Length"].isna().any() or counts["Length"].le(0).any():
        raise ValueError("All genes require a positive featureCounts Length.")

    counts = counts.copy().rename(columns={"Geneid": "gene_id", "Length": "length_bp"})
    library_rows = []
    for strain, replicates in SAMPLE_REPLICATES.items():
        for replicate in replicates:
            raw = pd.to_numeric(counts[replicate], errors="raise").astype(float)
            library_size = float(raw.sum())
            counts[f"{replicate}_cpm"] = raw / library_size * 1_000_000
            rate = raw / counts["length_bp"]
            rate_sum = float(rate.sum())
            counts[f"{replicate}_tpm"] = rate / rate_sum * 1_000_000
            library_rows.append(
                {
                    "strain": strain,
                    "replicate": replicate,
                    "library_size": int(library_size),
                    "tpm_denominator": rate_sum,
                }
            )
    return counts, pd.DataFrame(library_rows)


def build_expression_table(
    classes: pd.DataFrame, counts: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_expression = []
    candidate_expression = []
    for strain, replicates in SAMPLE_REPLICATES.items():
        strain_all = counts[["gene_id", "length_bp", *replicates]].copy()
        cpm_columns = [f"{replicate}_cpm" for replicate in replicates]
        tpm_columns = [f"{replicate}_tpm" for replicate in replicates]
        strain_all["strain"] = strain
        strain_all["mean_cpm"] = counts[cpm_columns].mean(axis=1)
        strain_all["mean_tpm"] = counts[tpm_columns].mean(axis=1)
        strain_all["log2_mean_cpm_plus1"] = np.log2(strain_all["mean_cpm"] + 1)
        strain_all["log2_mean_tpm_plus1"] = np.log2(strain_all["mean_tpm"] + 1)
        strain_all["expression_percentile"] = strain_all["mean_tpm"].rank(
            pct=True, method="average"
        )
        strain_all["bottom_decile"] = strain_all["expression_percentile"].le(0.10)
        strain_all["bottom_quartile"] = strain_all["expression_percentile"].le(0.25)
        all_expression.append(strain_all)

        strain_classes = classes[
            classes["strain"].eq(strain)
            & classes["gene_class"].isin(["fixed_only", "polymorphic"])
        ]
        merged = strain_classes.merge(strain_all, on=["strain", "gene_id"], how="left")
        if merged["mean_tpm"].isna().any():
            missing = merged.loc[merged["mean_tpm"].isna(), "gene_id"].head().tolist()
            raise ValueError(f"Eligible {strain} genes missing from count matrix: {missing}")
        candidate_expression.append(merged)
    return pd.concat(candidate_expression, ignore_index=True), pd.concat(all_expression, ignore_index=True)


def bootstrap_difference(
    fixed: np.ndarray,
    polymorphic: np.ndarray,
    rng: np.random.Generator,
    n_bootstraps: int,
) -> tuple[float, float]:
    differences = np.empty(n_bootstraps)
    for index in range(n_bootstraps):
        differences[index] = (
            rng.choice(polymorphic, len(polymorphic), replace=True).mean()
            - rng.choice(fixed, len(fixed), replace=True).mean()
        )
    return tuple(np.quantile(differences, [0.025, 0.975]))


def permutation_test(
    values: np.ndarray,
    labels: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int,
) -> tuple[float, float]:
    observed = values[labels].mean() - values[~labels].mean()
    null = np.empty(n_permutations)
    for index in range(n_permutations):
        permuted = rng.permutation(labels)
        null[index] = values[permuted].mean() - values[~permuted].mean()
    p_value = (np.count_nonzero(np.abs(null) >= abs(observed)) + 1) / (
        n_permutations + 1
    )
    return float(observed), float(p_value)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def continuous_tests(
    genes: pd.DataFrame,
    rng: np.random.Generator,
    n_permutations: int,
    n_bootstraps: int,
) -> pd.DataFrame:
    rows = []
    for strain, group in genes.groupby("strain", sort=False):
        for endpoint in ("log2_mean_tpm_plus1", "expression_percentile"):
            fixed = group.loc[group["gene_class"].eq("fixed_only"), endpoint].to_numpy()
            poly = group.loc[group["gene_class"].eq("polymorphic"), endpoint].to_numpy()
            values = group[endpoint].to_numpy()
            labels = group["gene_class"].eq("polymorphic").to_numpy()
            difference, permutation_p = permutation_test(
                values, labels, rng, n_permutations
            )
            mann_whitney = stats.mannwhitneyu(poly, fixed, alternative="two-sided")
            ci_low, ci_high = bootstrap_difference(
                fixed, poly, rng, n_bootstraps
            )
            cliffs_delta = 2 * float(mann_whitney.statistic) / (len(poly) * len(fixed)) - 1
            rows.append(
                {
                    "strain": strain,
                    "endpoint": endpoint,
                    "fixed_n": len(fixed),
                    "polymorphic_n": len(poly),
                    "fixed_mean": float(fixed.mean()),
                    "polymorphic_mean": float(poly.mean()),
                    "polymorphic_minus_fixed": difference,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "permutation_p_two_sided": permutation_p,
                    "mannwhitney_u": float(mann_whitney.statistic),
                    "mannwhitney_p_two_sided": float(mann_whitney.pvalue),
                    "cliffs_delta": cliffs_delta,
                    "permutations": n_permutations,
                    "bootstraps": n_bootstraps,
                }
            )
    results = pd.DataFrame(rows)
    results["permutation_fdr_bh"] = results.groupby("endpoint")[
        "permutation_p_two_sided"
    ].transform(benjamini_hochberg)
    results["mannwhitney_fdr_bh"] = results.groupby("endpoint")[
        "mannwhitney_p_two_sided"
    ].transform(benjamini_hochberg)
    return results


def binary_tests(genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for strain, group in genes.groupby("strain", sort=False):
        for endpoint in ("bottom_decile", "bottom_quartile"):
            fixed = group[group["gene_class"].eq("fixed_only")]
            poly = group[group["gene_class"].eq("polymorphic")]
            fixed_yes = int(fixed[endpoint].sum())
            poly_yes = int(poly[endpoint].sum())
            table = [
                [poly_yes, len(poly) - poly_yes],
                [fixed_yes, len(fixed) - fixed_yes],
            ]
            odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
            rows.append(
                {
                    "strain": strain,
                    "endpoint": endpoint,
                    "fixed_n": len(fixed),
                    "polymorphic_n": len(poly),
                    "fixed_low_n": fixed_yes,
                    "polymorphic_low_n": poly_yes,
                    "fixed_low_fraction": fixed_yes / len(fixed),
                    "polymorphic_low_fraction": poly_yes / len(poly),
                    "odds_ratio_polymorphic": float(odds_ratio),
                    "fisher_p_two_sided": float(p_value),
                }
            )
    results = pd.DataFrame(rows)
    # These are secondary, correlated threshold tests; correct across all six
    # strain-by-threshold comparisons rather than separately by threshold.
    results["fisher_fdr_bh"] = benjamini_hochberg(results["fisher_p_two_sided"])
    return results


def fit_count_gee(
    genes: pd.DataFrame,
    counts: pd.DataFrame,
    libraries: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    library_sizes = libraries.set_index("replicate")["library_size"].to_dict()
    raw_columns = [rep for reps in SAMPLE_REPLICATES.values() for rep in reps]
    raw = counts[["gene_id", "length_bp", *raw_columns]]
    for strain, replicates in SAMPLE_REPLICATES.items():
        candidate = genes.loc[
            genes["strain"].eq(strain), ["gene_id", "gene_class", "length_bp"]
        ].merge(
            raw, on=["gene_id", "length_bp"], how="left", validate="one_to_one"
        )
        long = candidate.melt(
            id_vars=["gene_id", "gene_class", "length_bp"],
            value_vars=list(replicates),
            var_name="replicate",
            value_name="raw_count",
        )
        long["is_polymorphic"] = long["gene_class"].eq("polymorphic").astype(int)
        long["offset"] = (
            long["replicate"].map(library_sizes).astype(float).map(math.log)
            + (long["length_bp"] / 1_000).map(math.log)
        )
        try:
            model = sm.GEE.from_formula(
                "raw_count ~ is_polymorphic + C(replicate)",
                groups="gene_id",
                data=long,
                offset=long["offset"],
                cov_struct=Exchangeable(),
                family=sm.families.NegativeBinomial(alpha=1.0),
            ).fit(maxiter=100)
            coefficient = float(model.params["is_polymorphic"])
            standard_error = float(model.bse["is_polymorphic"])
            p_value = float(model.pvalues["is_polymorphic"])
            ci = model.conf_int().loc["is_polymorphic"]
            rows.append(
                {
                    "strain": strain,
                    "gene_clusters": int(long["gene_id"].nunique()),
                    "count_observations": len(long),
                    "log_rate_ratio_polymorphic": coefficient,
                    "standard_error": standard_error,
                    "rate_ratio_polymorphic": math.exp(coefficient),
                    "rate_ratio_ci_low": math.exp(float(ci.iloc[0])),
                    "rate_ratio_ci_high": math.exp(float(ci.iloc[1])),
                    "p_value": p_value,
                    "converged": bool(model.converged),
                }
            )
        except Exception as error:  # preserve a readable failure in exploratory output
            rows.append(
                {
                    "strain": strain,
                    "gene_clusters": int(long["gene_id"].nunique()),
                    "count_observations": len(long),
                    "log_rate_ratio_polymorphic": math.nan,
                    "standard_error": math.nan,
                    "rate_ratio_polymorphic": math.nan,
                    "rate_ratio_ci_low": math.nan,
                    "rate_ratio_ci_high": math.nan,
                    "p_value": math.nan,
                    "converged": False,
                    "error": str(error),
                }
            )
    results = pd.DataFrame(rows)
    finite = results["p_value"].notna()
    results.loc[finite, "fdr_bh"] = benjamini_hochberg(results.loc[finite, "p_value"])
    return results


def summarize(genes: pd.DataFrame) -> pd.DataFrame:
    return (
        genes.groupby(["strain", "gene_class"])
        .agg(
            n_genes=("gene_id", "size"),
            median_tpm=("mean_tpm", "median"),
            mean_tpm=("mean_tpm", "mean"),
            median_log2_tpm_plus1=("log2_mean_tpm_plus1", "median"),
            mean_log2_tpm_plus1=("log2_mean_tpm_plus1", "mean"),
            median_expression_percentile=("expression_percentile", "median"),
            zero_expression_genes=("mean_tpm", lambda x: int(x.eq(0).sum())),
            median_length_bp=("length_bp", "median"),
            median_present_introner_loci=("n_present_introner_loci", "median"),
        )
        .reset_index()
    )


def plot_expression(genes: pd.DataFrame, output: Path) -> None:
    rng = np.random.default_rng(1545)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharey=True)
    for ax, strain in zip(axes, SAMPLE_REPLICATES):
        subset = genes[genes["strain"].eq(strain)]
        arrays = [
            subset.loc[subset["gene_class"].eq(gene_class), "log2_mean_tpm_plus1"].to_numpy()
            for gene_class in ("fixed_only", "polymorphic")
        ]
        violin = ax.violinplot(arrays, positions=[0, 1], showextrema=False, widths=0.8)
        for body, color in zip(violin["bodies"], [FIXED_COLOR, POLYMORPHIC_COLOR]):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.45)
        box = ax.boxplot(
            arrays,
            positions=[0, 1],
            widths=0.25,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": "black"},
        )
        for patch, color in zip(box["boxes"], [FIXED_COLOR, POLYMORPHIC_COLOR]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for position, values, color in zip(
            [0, 1], arrays, [FIXED_COLOR, POLYMORPHIC_COLOR]
        ):
            sampled = rng.choice(values, size=min(500, len(values)), replace=False)
            ax.scatter(
                rng.normal(position, 0.045, len(sampled)),
                sampled,
                s=7,
                alpha=0.12,
                linewidth=0,
                color=color,
            )
        ax.set_title(strain)
        ax.set_xticks([0, 1], [f"Fixed only\n(n={len(arrays[0])})", f"Polymorphic\n(n={len(arrays[1])})"])
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("log2(mean TPM + 1)")
    fig.suptitle("Expression of genes containing fixed-only or polymorphic introners")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_results(
    output: Path,
    summaries: pd.DataFrame,
    continuous: pd.DataFrame,
    binary: pd.DataFrame,
    gee: pd.DataFrame,
    libraries: pd.DataFrame,
) -> None:
    primary = continuous[continuous["endpoint"].eq("log2_mean_tpm_plus1")]
    lines = [
        "# RNA-seq count-matrix expression comparison",
        "",
        "The analysis uses all eligible genes in `merged_counts_matrix.csv`, performs",
        "each comparison separately by strain, and does not pool expression values",
        "across strains. TPM corrects for both library size and featureCounts gene",
        "length. The negative-binomial GEE uses replicate raw counts, library-size",
        "and gene-length offsets, replicate fixed effects, and gene-clustered errors.",
        "",
        "## Primary gene-level results",
        "",
        "| Strain | Fixed n | Polymorphic n | Poly - fixed log2(TPM+1) | 95% bootstrap CI | Permutation P | BH FDR |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in primary.itertuples():
        lines.append(
            f"| {row.strain} | {int(row.fixed_n)} | {int(row.polymorphic_n)} | "
            f"{row.polymorphic_minus_fixed:.4f} | {row.bootstrap_ci_low:.4f} to "
            f"{row.bootstrap_ci_high:.4f} | {row.permutation_p_two_sided:.4g} | "
            f"{row.permutation_fdr_bh:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Replicate-level negative-binomial GEE sensitivity analysis",
            "",
            "| Strain | Polymorphic/fixed expression-rate ratio | 95% CI | P | BH FDR |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in gee.itertuples():
        lines.append(
            f"| {row.strain} | {row.rate_ratio_polymorphic:.4f} | "
            f"{row.rate_ratio_ci_low:.4f} to {row.rate_ratio_ci_high:.4f} | "
            f"{row.p_value:.4g} | {row.fdr_bh:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Secondary low-expression threshold tests",
            "",
            "| Strain | Threshold | Fixed fraction | Polymorphic fraction | Odds ratio | Fisher P | BH FDR across six tests |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in binary.itertuples():
        lines.append(
            f"| {row.strain} | {row.endpoint.replace('_', ' ')} | "
            f"{row.fixed_low_fraction:.4f} | {row.polymorphic_low_fraction:.4f} | "
            f"{row.odds_ratio_polymorphic:.4f} | {row.fisher_p_two_sided:.4g} | "
            f"{row.fisher_fdr_bh:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `fixed_only` is strict: every present, gene-assigned introner must be confidently fixed-present in the focal population.",
            "- `polymorphic` means that the sampled strain carries at least one fully callable polymorphic introner; these genes may also carry fixed introners.",
            "- Mating-type-region orthologs are excluded using the CCMP1545 interval.",
            "- Expression percentiles and low-expression categories are defined against all 10,648 count-matrix genes within each strain.",
            "- Group 1 and Group 2 use their respective population-specific introner classifications.",
        ]
    )
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    if args.permutations < 1 or args.bootstraps < 1:
        raise ValueError("--permutations and --bootstraps must be positive.")
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrix = load_genotypes(args.matrix)
    classes, class_audit = build_gene_classes(matrix)
    normalized_counts, libraries = normalize_counts(pd.read_csv(args.counts))
    genes, all_expression = build_expression_table(classes, normalized_counts)

    rng = np.random.default_rng(args.seed)
    summaries = summarize(genes)
    continuous = continuous_tests(
        genes, rng, args.permutations, args.bootstraps
    )
    binary = binary_tests(genes)
    gee = fit_count_gee(genes, normalized_counts, libraries)

    genes.to_csv(args.outdir / "gene_class_expression.tsv", sep="\t", index=False)
    all_expression.to_csv(
        args.outdir / "all_count_matrix_gene_expression.tsv", sep="\t", index=False
    )
    summaries.to_csv(args.outdir / "sample_gene_class_summary.tsv", sep="\t", index=False)
    continuous.to_csv(args.outdir / "continuous_expression_tests.tsv", sep="\t", index=False)
    binary.to_csv(args.outdir / "low_expression_tests.tsv", sep="\t", index=False)
    gee.to_csv(args.outdir / "negative_binomial_gee.tsv", sep="\t", index=False)
    libraries.to_csv(args.outdir / "library_sizes.tsv", sep="\t", index=False)
    class_audit.to_csv(args.outdir / "class_assignment_audit.tsv", sep="\t", index=False)
    plot_expression(genes, args.outdir / "gene_expression_distributions.png")
    write_results(
        args.outdir / "RESULTS.md",
        summaries,
        continuous,
        binary,
        gee,
        libraries,
    )

    print(f"Eligible gene-strain observations: {len(genes)}")
    print(summaries[["strain", "gene_class", "n_genes"]].to_string(index=False))
    print(f"Wrote results to {args.outdir}")


if __name__ == "__main__":
    main()

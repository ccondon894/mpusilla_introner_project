#!/usr/bin/env python3
"""Fit reference-independent within-gene introner-count/isoform models."""

from __future__ import annotations

import argparse
import os
import warnings
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalPoisson


STRAINS = ("CCMP1545", "RCC1614", "RCC1749")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-long-read-count", type=int, default=10)
    return parser.parse_args()


def supported_rows(data: pd.DataFrame, threshold: int) -> pd.DataFrame:
    required = [
        "common_gene_id",
        "strain",
        "n_extra_isoforms",
        "current_introner_count",
        "current_missing_locus_count",
        "total_long_read_count",
        "log1p_long_read_count",
    ]
    frame = data.dropna(subset=required).copy()
    frame = frame[
        (frame["strain"].isin(STRAINS))
        & (frame["total_long_read_count"] >= threshold)
        & (frame["current_missing_locus_count"] == 0)
    ].copy()
    strain_counts = frame.groupby("common_gene_id")["strain"].nunique()
    return frame[frame["common_gene_id"].isin(strain_counts[strain_counts >= 2].index)].copy()


def fit_conditional_poisson(frame: pd.DataFrame) -> tuple[pd.DataFrame, object]:
    frame = frame.copy()
    frame["log_reads_within_gene"] = frame["log1p_long_read_count"] - frame.groupby(
        "common_gene_id"
    )["log1p_long_read_count"].transform("mean")
    # Strain indicators are nuisance adjustments only; the introner-count
    # coefficient is based on within-gene changes and has no reference-locus label.
    frame["strain_RCC1614"] = (frame["strain"] == "RCC1614").astype(float)
    frame["strain_RCC1749"] = (frame["strain"] == "RCC1749").astype(float)
    terms = [
        "current_introner_count",
        "strain_RCC1614",
        "strain_RCC1749",
        "log_reads_within_gene",
    ]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = ConditionalPoisson(
            frame["n_extra_isoforms"].astype(float),
            frame[terms].astype(float),
            groups=frame["common_gene_id"],
        ).fit(method="BFGS", maxiter=500, disp=False)
    params = pd.Series(np.asarray(result.params), index=terms)
    errors = pd.Series(np.asarray(result.bse), index=terms)
    pvalues = pd.Series(np.asarray(result.pvalues), index=terms)
    rows = []
    for term in terms:
        estimate = float(params[term])
        se = float(errors[term])
        rows.append(
            {
                "term": term,
                "coefficient": estimate,
                "std_error": se,
                "p_value": float(pvalues[term]),
                "rate_ratio": float(np.exp(estimate)),
                "rate_ratio_ci_lower": float(np.exp(estimate - 1.96 * se)),
                "rate_ratio_ci_upper": float(np.exp(estimate + 1.96 * se)),
                "n_observations": int(result.model.nobs),
                "n_genes": len(result.model._endog_grp),
                "fit_warning": "; ".join(sorted({str(x.message) for x in caught})),
            }
        )
    return pd.DataFrame(rows), result


def make_pairwise(frame: pd.DataFrame) -> pd.DataFrame:
    indexed = {
        (row.common_gene_id, row.strain): row
        for row in frame.itertuples(index=False)
    }
    rows = []
    genes = sorted(frame["common_gene_id"].unique())
    for gene in genes:
        for strain_a, strain_b in combinations(STRAINS, 2):
            a = indexed.get((gene, strain_a))
            b = indexed.get((gene, strain_b))
            if a is None or b is None:
                continue
            delta_introner = float(b.current_introner_count - a.current_introner_count)
            delta_isoform = float(b.n_extra_isoforms - a.n_extra_isoforms)
            delta_reads = float(b.log1p_long_read_count - a.log1p_long_read_count)
            if delta_introner == 0:
                oriented_isoform = np.nan
                oriented_reads = np.nan
                absolute_change = 0
            else:
                direction = np.sign(delta_introner)
                oriented_isoform = delta_isoform * direction
                oriented_reads = delta_reads * direction
                absolute_change = int(abs(delta_introner))
            rows.append(
                {
                    "common_gene_id": gene,
                    "strain_pair": f"{strain_a}_vs_{strain_b}",
                    "delta_introner_count": delta_introner,
                    "delta_extra_isoforms": delta_isoform,
                    "delta_log1p_long_read_count": delta_reads,
                    "absolute_introner_change": absolute_change,
                    "extra_isoform_change_higher_minus_lower_introner": oriented_isoform,
                    "read_change_higher_minus_lower_introner": oriented_reads,
                }
            )
    return pd.DataFrame(rows)


def fit_pairwise(pairwise: pd.DataFrame) -> pd.DataFrame:
    result = smf.ols(
        "delta_extra_isoforms ~ delta_introner_count + "
        "delta_log1p_long_read_count + C(strain_pair)",
        data=pairwise,
    ).fit(cov_type="cluster", cov_kwds={"groups": pairwise["common_gene_id"]})
    rows = []
    for term in result.params.index:
        rows.append(
            {
                "term": term,
                "coefficient": float(result.params[term]),
                "std_error": float(result.bse[term]),
                "p_value": float(result.pvalues[term]),
                "ci_lower": float(result.conf_int().loc[term, 0]),
                "ci_upper": float(result.conf_int().loc[term, 1]),
                "n_pairwise_comparisons": int(result.nobs),
                "n_genes": int(pairwise["common_gene_id"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def fit_within_between_negative_binomial(frame: pd.DataFrame) -> pd.DataFrame:
    """Separate within-gene change from average cross-gene burden in one GEE."""
    model_data = frame.dropna(subset=["log_gene_span"]).copy()
    for variable in ("current_introner_count", "log1p_long_read_count"):
        between = model_data.groupby("common_gene_id")[variable].transform("mean")
        model_data[f"{variable}_between"] = between
        model_data[f"{variable}_within"] = model_data[variable] - between
    formula = (
        "n_extra_isoforms ~ current_introner_count_within + "
        "current_introner_count_between + C(strain) + "
        "log1p_long_read_count_within + log1p_long_read_count_between + "
        "log_gene_span"
    )
    result = smf.gee(
        formula,
        groups="common_gene_id",
        data=model_data,
        family=sm.families.NegativeBinomial(alpha=1.0),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit()
    rows = []
    for term in result.params.index:
        estimate = float(result.params[term])
        se = float(result.bse[term])
        rows.append(
            {
                "term": term,
                "coefficient": estimate,
                "std_error": se,
                "p_value": float(result.pvalues[term]),
                "rate_ratio": float(np.exp(estimate)),
                "rate_ratio_ci_lower": float(np.exp(estimate - 1.96 * se)),
                "rate_ratio_ci_upper": float(np.exp(estimate + 1.96 * se)),
                "n_observations": len(model_data),
                "n_genes": model_data["common_gene_id"].nunique(),
            }
        )
    return pd.DataFrame(rows)


def plot_changes(pairwise: pd.DataFrame, output: Path) -> None:
    changed = pairwise[pairwise["absolute_introner_change"] > 0].copy()
    changed = changed[changed["absolute_introner_change"] <= 4]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    groups = sorted(changed["absolute_introner_change"].unique())
    values = [
        changed.loc[
            changed["absolute_introner_change"] == group,
            "extra_isoform_change_higher_minus_lower_introner",
        ].values
        for group in groups
    ]
    if values:
        ax.boxplot(values, tick_labels=groups, showfliers=False)
        rng = np.random.default_rng(20260716)
        for index, vals in enumerate(values, start=1):
            sample = vals if len(vals) <= 500 else rng.choice(vals, 500, replace=False)
            ax.scatter(
                rng.normal(index, 0.055, len(sample)),
                sample,
                s=8,
                alpha=0.22,
                color="#336699",
                linewidth=0,
            )
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Absolute within-gene change in introner count")
    ax.set_ylabel("Change in extra isoforms\n(higher- minus lower-introner strain)")
    ax.set_title("Reference-independent pairwise changes")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input, sep="\t")

    sensitivity = []
    negative_binomial_sensitivity = []
    primary_frame = None
    primary_coefficients = None
    primary_nb_coefficients = None
    for threshold in (5, 10, 20):
        frame = supported_rows(data, threshold)
        coefficients, _ = fit_conditional_poisson(frame)
        focal = coefficients[coefficients["term"] == "current_introner_count"].copy()
        focal.insert(0, "minimum_long_read_count", threshold)
        sensitivity.append(focal)
        nb_coefficients = fit_within_between_negative_binomial(frame)
        nb_focal = nb_coefficients[
            nb_coefficients["term"].isin(
                ["current_introner_count_within", "current_introner_count_between"]
            )
        ].copy()
        nb_focal.insert(0, "minimum_long_read_count", threshold)
        negative_binomial_sensitivity.append(nb_focal)
        if threshold == args.minimum_long_read_count:
            primary_frame = frame
            primary_coefficients = coefficients
            primary_nb_coefficients = nb_coefficients

    if (
        primary_frame is None
        or primary_coefficients is None
        or primary_nb_coefficients is None
    ):
        raise SystemExit("Primary support threshold must be one of 5, 10, or 20")

    pairwise = make_pairwise(primary_frame)
    pairwise_coefficients = fit_pairwise(pairwise)
    changed = pairwise[pairwise["absolute_introner_change"] > 0].copy()
    pairwise_summary = (
        changed.groupby("absolute_introner_change")
        .agg(
            n_comparisons=("common_gene_id", "size"),
            n_genes=("common_gene_id", "nunique"),
            mean_extra_isoform_change=(
                "extra_isoform_change_higher_minus_lower_introner",
                "mean",
            ),
            median_extra_isoform_change=(
                "extra_isoform_change_higher_minus_lower_introner",
                "median",
            ),
            fraction_higher_isoform_count=(
                "extra_isoform_change_higher_minus_lower_introner",
                lambda x: float((x > 0).mean()),
            ),
            fraction_lower_isoform_count=(
                "extra_isoform_change_higher_minus_lower_introner",
                lambda x: float((x < 0).mean()),
            ),
        )
        .reset_index()
    )

    primary_coefficients.to_csv(
        args.output_dir / "conditional_poisson_coefficients.tsv", sep="\t", index=False
    )
    pd.concat(sensitivity, ignore_index=True).to_csv(
        args.output_dir / "support_threshold_sensitivity.tsv", sep="\t", index=False
    )
    primary_nb_coefficients.to_csv(
        args.output_dir / "within_between_negative_binomial.tsv", sep="\t", index=False
    )
    pd.concat(negative_binomial_sensitivity, ignore_index=True).to_csv(
        args.output_dir / "within_between_negative_binomial_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    pairwise.to_csv(args.output_dir / "pairwise_gene_changes.tsv", sep="\t", index=False)
    pairwise_coefficients.to_csv(
        args.output_dir / "pairwise_change_model.tsv", sep="\t", index=False
    )
    pairwise_summary.to_csv(
        args.output_dir / "pairwise_change_summary.tsv", sep="\t", index=False
    )
    plot_changes(pairwise, args.output_dir / "pairwise_isoform_change.png")

    focal = primary_coefficients.loc[
        primary_coefficients["term"] == "current_introner_count"
    ].iloc[0]
    pair_focal = pairwise_coefficients.loc[
        pairwise_coefficients["term"] == "delta_introner_count"
    ].iloc[0]
    nb_within = primary_nb_coefficients.loc[
        primary_nb_coefficients["term"] == "current_introner_count_within"
    ].iloc[0]
    nb_between = primary_nb_coefficients.loc[
        primary_nb_coefficients["term"] == "current_introner_count_between"
    ].iloc[0]
    n_changed_genes = changed["common_gene_id"].nunique()
    lines = [
        "Reference-independent within-gene isoform model",
        "================================================",
        f"Minimum long-read support: {args.minimum_long_read_count}",
        f"Supported observations before conditional-model invariant strata are dropped: {len(primary_frame)}",
        f"Supported genes with at least two strains: {primary_frame['common_gene_id'].nunique()}",
        f"Genes with at least one pairwise introner-count difference: {n_changed_genes}",
        "",
        "Gene-conditioned Poisson model",
        f"Rate ratio per one-introner within-gene increase: {focal.rate_ratio:.6f}",
        f"95% CI: {focal.rate_ratio_ci_lower:.6f}-{focal.rate_ratio_ci_upper:.6f}",
        f"Wald p-value: {focal.p_value:.6g}",
        f"Conditional-model observations: {int(focal.n_observations)}",
        f"Conditional-model genes: {int(focal.n_genes)}",
        "",
        "Pairwise first-difference OLS check (gene-clustered SE)",
        f"Extra-isoform change per one-introner change: {pair_focal.coefficient:.6f}",
        f"95% CI: {pair_focal.ci_lower:.6f}-{pair_focal.ci_upper:.6f}",
        f"Wald p-value: {pair_focal.p_value:.6g}",
        "",
        "Negative-binomial GEE within/between robustness model",
        f"Within-gene rate ratio: {nb_within.rate_ratio:.6f}",
        f"Within-gene 95% CI: {nb_within.rate_ratio_ci_lower:.6f}-{nb_within.rate_ratio_ci_upper:.6f}",
        f"Within-gene Wald p-value: {nb_within.p_value:.6g}",
        f"Between-gene burden rate ratio: {nb_between.rate_ratio:.6f}",
        f"Between-gene 95% CI: {nb_between.rate_ratio_ci_lower:.6f}-{nb_between.rate_ratio_ci_upper:.6f}",
        f"Between-gene Wald p-value: {nb_between.p_value:.6g}",
        "",
        "Interpretation: the focal coefficient uses current introner count and",
        "within-gene differences; it does not label any strain as the ancestral",
        "state or interpret CCMP1545-relative differences as gains or losses.",
    ]
    (args.output_dir / "summary.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

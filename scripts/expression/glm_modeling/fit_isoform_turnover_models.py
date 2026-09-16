"""Fit reference-independent within-gene and cross-gene isoform models."""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalLogit, ConditionalPoisson


STRAINS = ("CCMP1545", "RCC1614", "RCC1749")
TERM_LABELS = {
    "current_introner_count": "Current introner count",
}


def result_rows(model_name, result, terms, scale, nobs, ngenes, note=""):
    rows = []
    params = pd.Series(np.asarray(result.params), index=terms)
    errors = pd.Series(np.asarray(result.bse), index=terms)
    pvalues = pd.Series(np.asarray(result.pvalues), index=terms)
    for term in terms:
        coefficient = float(params[term])
        standard_error = float(errors[term])
        lower = coefficient - 1.96 * standard_error
        upper = coefficient + 1.96 * standard_error
        rows.append(
            {
                "model": model_name,
                "term": term,
                "coefficient": coefficient,
                "std_error": standard_error,
                "p_value": float(pvalues[term]),
                "ci_lower": lower,
                "ci_upper": upper,
                "effect_scale": scale,
                "effect_ratio": float(np.exp(coefficient)),
                "effect_ratio_ci_lower": float(np.exp(lower)),
                "effect_ratio_ci_upper": float(np.exp(upper)),
                "n_observations": int(nobs),
                "n_genes": int(ngenes),
                "note": note,
            }
        )
    return rows


def prepare_conditional_data(data, strains):
    frame = data[data["strain"].isin(strains)].copy()
    counts = frame.groupby("common_gene_id")["strain"].nunique()
    frame = frame[frame["common_gene_id"].isin(counts[counts >= 2].index)].copy()
    frame["log_reads_within_gene"] = frame["log1p_long_read_count"] - frame.groupby(
        "common_gene_id"
    )["log1p_long_read_count"].transform("mean")
    for strain in strains:
        if strain != "CCMP1545":
            frame[f"strain_{strain}"] = (frame["strain"] == strain).astype(float)
    return frame


def fit_conditional_poisson(frame, model_name):
    terms = ["current_introner_count"] + [
        f"strain_{strain}" for strain in STRAINS if f"strain_{strain}" in frame.columns
    ] + ["log_reads_within_gene"]
    exog = frame[terms].astype(float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = ConditionalPoisson(
            frame["n_extra_isoforms"].astype(float),
            exog,
            groups=frame["common_gene_id"],
        ).fit(method="BFGS", maxiter=500, disp=False)
    note = "; ".join(sorted({str(item.message) for item in caught}))
    rows = result_rows(
        model_name,
        result,
        terms,
        "rate_ratio",
        int(result.model.nobs),
        len(result.model._endog_grp),
        note,
    )
    return rows, result


def fit_conditional_logit(frame, model_name):
    terms = ["current_introner_count"] + [
        f"strain_{strain}" for strain in STRAINS if f"strain_{strain}" in frame.columns
    ] + ["log_reads_within_gene"]
    exog = frame[terms].astype(float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = ConditionalLogit(
            frame["has_alternative_isoform"].astype(float),
            exog,
            groups=frame["common_gene_id"],
        ).fit(method="BFGS", maxiter=500, disp=False)
    note = "; ".join(sorted({str(item.message) for item in caught}))
    rows = result_rows(
        model_name,
        result,
        terms,
        "odds_ratio",
        int(result.model.nobs),
        len(result.model._endog_grp),
        note,
    )
    return rows, result


def fit_architecture_models(data):
    primary_formula = (
        "{outcome} ~ current_introner_count + C(strain) + log1p_long_read_count + log_gene_span"
    )
    exon_adjusted_formula = primary_formula + " + max_annotated_exons"
    rows = []
    fitted = {}

    count_model = smf.gee(
        primary_formula.format(outcome="n_extra_isoforms"),
        groups="common_gene_id",
        data=data,
        family=sm.families.NegativeBinomial(alpha=1.0),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit()
    terms = list(count_model.params.index)
    rows.extend(
        result_rows(
            "architecture_negative_binomial_gee",
            count_model,
            terms,
            "rate_ratio",
            len(data),
            data["common_gene_id"].nunique(),
            "Negative-binomial GEE with alpha fixed at 1; gene-clustered working correlation.",
        )
    )
    fitted["architecture_negative_binomial_gee"] = count_model

    exon_adjusted_count = smf.gee(
        exon_adjusted_formula.format(outcome="n_extra_isoforms"),
        groups="common_gene_id",
        data=data,
        family=sm.families.NegativeBinomial(alpha=1.0),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit()
    terms = list(exon_adjusted_count.params.index)
    rows.extend(
        result_rows(
            "architecture_negative_binomial_exon_adjusted_gee",
            exon_adjusted_count,
            terms,
            "rate_ratio",
            len(data),
            data["common_gene_id"].nunique(),
            "Sensitivity model conditioning on maximum annotated exon count.",
        )
    )
    fitted["architecture_negative_binomial_exon_adjusted_gee"] = exon_adjusted_count

    hurdle_model = smf.gee(
        primary_formula.format(outcome="has_alternative_isoform"),
        groups="common_gene_id",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit()
    terms = list(hurdle_model.params.index)
    rows.extend(
        result_rows(
            "architecture_any_alternative_isoform_gee",
            hurdle_model,
            terms,
            "odds_ratio",
            len(data),
            data["common_gene_id"].nunique(),
            "Gene-clustered logistic GEE.",
        )
    )
    fitted["architecture_any_alternative_isoform_gee"] = hurdle_model

    exon_adjusted_hurdle = smf.gee(
        exon_adjusted_formula.format(outcome="has_alternative_isoform"),
        groups="common_gene_id",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit()
    terms = list(exon_adjusted_hurdle.params.index)
    rows.extend(
        result_rows(
            "architecture_any_alternative_isoform_exon_adjusted_gee",
            exon_adjusted_hurdle,
            terms,
            "odds_ratio",
            len(data),
            data["common_gene_id"].nunique(),
            "Sensitivity model conditioning on maximum annotated exon count.",
        )
    )
    fitted["architecture_any_alternative_isoform_exon_adjusted_gee"] = exon_adjusted_hurdle

    positive = data[data["n_extra_isoforms"] > 0].copy()
    positive["additional_isoforms_after_first"] = positive["n_extra_isoforms"] - 1
    positive_model = smf.glm(
        primary_formula.format(outcome="additional_isoforms_after_first"),
        data=positive,
        family=sm.families.Poisson(),
    ).fit(cov_type="cluster", cov_kwds={"groups": positive["common_gene_id"]})
    terms = list(positive_model.params.index)
    rows.extend(
        result_rows(
            "architecture_positive_count_clustered_poisson",
            positive_model,
            terms,
            "rate_ratio",
            len(positive),
            positive["common_gene_id"].nunique(),
            "Cluster-robust Poisson GLM for additional isoforms beyond the first alternative isoform; descriptive hurdle component.",
        )
    )
    fitted["architecture_positive_count_clustered_poisson"] = positive_model
    return rows, fitted


def coefficient_plot(coefficients, output):
    selected_models = [
        "within_gene_conditional_poisson_all",
        "within_gene_conditional_logit_all",
        "architecture_negative_binomial_gee",
        "architecture_negative_binomial_exon_adjusted_gee",
        "architecture_any_alternative_isoform_gee",
    ]
    selected = coefficients[
        coefficients["model"].isin(selected_models)
        & (coefficients["term"] == "current_introner_count")
    ].copy()
    if selected.empty:
        return
    selected["label"] = selected["term"].map(TERM_LABELS)
    selected["model_label"] = selected["model"].map(
        {
            "within_gene_conditional_poisson_all": "Within gene: extra-isoform count",
            "within_gene_conditional_logit_all": "Within gene: any alternative isoform",
            "architecture_negative_binomial_gee": "Architecture: extra-isoform count",
            "architecture_negative_binomial_exon_adjusted_gee": "Architecture: count, exon-adjusted",
            "architecture_any_alternative_isoform_gee": "Architecture: any alternative isoform",
        }
    )
    fig, ax = plt.subplots(figsize=(9, 5.5))
    y_positions = np.arange(len(selected))
    colors = selected["model"].map(
        {
            "within_gene_conditional_poisson_all": "#0072B2",
            "within_gene_conditional_logit_all": "#56B4E9",
            "architecture_negative_binomial_gee": "#009E73",
            "architecture_negative_binomial_exon_adjusted_gee": "#D55E00",
            "architecture_any_alternative_isoform_gee": "#CC79A7",
        }
    )
    for index, (_, row) in enumerate(selected.iterrows()):
        color = colors.iloc[index]
        ax.errorbar(
            row["effect_ratio"],
            y_positions[index],
            xerr=[[row["effect_ratio"] - row["effect_ratio_ci_lower"]],
                  [row["effect_ratio_ci_upper"] - row["effect_ratio"]]],
            fmt="o",
            color=color,
            capsize=3,
            linewidth=1.5,
            markersize=6,
        )
    ax.axvline(1, color="black", linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(selected["model_label"] + " — " + selected["label"])
    ax.set_xlabel("Rate ratio or odds ratio (95% CI)")
    ax.set_title("Introner number and isoform richness")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def overview_plot(data, coefficients, output):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for strain, group in data.groupby("strain"):
        axes[0, 0].hist(
            group["n_isoforms"], bins=np.arange(0.5, group["n_isoforms"].max() + 1.5),
            histtype="step", linewidth=1.7, label=strain
        )
    axes[0, 0].set_xlabel("Detected qualifying isoforms")
    axes[0, 0].set_ylabel("Gene-strain observations")
    axes[0, 0].legend(frameon=False)

    introner_counts = [
        data.loc[data["strain"] == strain, "current_introner_count"].values
        for strain in STRAINS
    ]
    axes[0, 1].boxplot(introner_counts, tick_labels=STRAINS, showfliers=False)
    axes[0, 1].set_ylabel("Current introner count per gene")
    axes[0, 1].set_xlabel("")
    axes[0, 1].tick_params(axis="x", rotation=25)

    axes[1, 0].scatter(
        data["total_long_read_count"], data["n_isoforms"], s=6, alpha=0.2, color="#444444"
    )
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xlabel("Long-read support per gene (log scale)")
    axes[1, 0].set_ylabel("Detected qualifying isoforms")

    architecture = coefficients[
        (coefficients["model"] == "architecture_negative_binomial_gee")
        & (coefficients["term"] == "current_introner_count")
    ]
    turnover = coefficients[
        (coefficients["model"] == "within_gene_conditional_poisson_all")
        & (coefficients["term"] == "current_introner_count")
    ]
    display = pd.concat([architecture, turnover], ignore_index=True)
    if not display.empty:
        labels = display["model"].map(
            {
                "within_gene_conditional_poisson_all": "Within-gene count change",
                "architecture_negative_binomial_gee": "Cross-gene burden",
            }
        )
        y = np.arange(len(display))
        axes[1, 1].errorbar(
            display["effect_ratio"], y,
            xerr=np.vstack([
                display["effect_ratio"] - display["effect_ratio_ci_lower"],
                display["effect_ratio_ci_upper"] - display["effect_ratio"],
            ]),
            fmt="o", color="#333333", capsize=3,
        )
        axes[1, 1].axvline(1, linestyle="--", color="black", linewidth=1)
        axes[1, 1].set_yticks(y)
        axes[1, 1].set_yticklabels(labels)
        axes[1, 1].set_xlabel("Count rate ratio (95% CI)")
        axes[1, 1].set_title("Primary count-model estimates")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "results" / "expression" / "isoform_analysis" / "turnover"
        / "isoform_introner_model_data.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "results" / "expression" / "glm_modeling" / "isoform_turnover",
    )
    parser.add_argument("--minimum-long-read-count", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw_data = pd.read_csv(args.input, sep="\t")
    architecture_required = [
        "n_extra_isoforms", "has_alternative_isoform", "log1p_long_read_count",
        "gene_span", "max_annotated_exons",
    ]
    within_gene_required = [
        "n_extra_isoforms", "has_alternative_isoform", "log1p_long_read_count",
        "current_introner_count", "total_long_read_count",
        "current_missing_locus_count", "common_gene_id", "strain",
    ]
    architecture_data = raw_data.dropna(subset=architecture_required).copy()
    within_gene_data = raw_data.dropna(subset=within_gene_required).copy()
    supported = architecture_data[
        (architecture_data["total_long_read_count"] >= args.minimum_long_read_count)
        & (architecture_data["current_missing_locus_count"] == 0)
    ].copy()
    within_gene_supported = within_gene_data[
        (within_gene_data["total_long_read_count"] >= args.minimum_long_read_count)
        & (within_gene_data["current_missing_locus_count"] == 0)
    ].copy()
    # Reference-independent within-gene analysis: current introner count is
    # compared across strains without assigning CCMP1545 as an ancestral state.
    within_gene = prepare_conditional_data(within_gene_supported, STRAINS)

    coefficient_rows = []
    model_notes = []
    fitted_models = {}
    for fitter, label in [
        (fit_conditional_poisson, "within_gene_conditional_poisson_all"),
        (fit_conditional_logit, "within_gene_conditional_logit_all"),
    ]:
        try:
            rows, result = fitter(within_gene, label)
            coefficient_rows.extend(rows)
            fitted_models[label] = result
        except Exception as error:
            model_notes.append(f"{label}: FAILED: {error}")

    try:
        rows, architecture_models = fit_architecture_models(supported)
        coefficient_rows.extend(rows)
        fitted_models.update(architecture_models)
    except Exception as error:
        model_notes.append(f"architecture models: FAILED: {error}")

    coefficients = pd.DataFrame(coefficient_rows)
    coefficients.to_csv(args.output_dir / "model_coefficients.tsv", sep="\t", index=False)

    support = (
        supported.groupby("strain")
        .agg(
            observations=("common_gene_id", "size"),
            genes=("common_gene_id", "nunique"),
            genes_with_introners=("current_introner_count", lambda values: int((values > 0).sum())),
            mean_current_introner_count=("current_introner_count", "mean"),
            genes_with_gain=("reference_relative_gain_count", lambda values: int((values > 0).sum())),
            genes_with_loss=("reference_relative_loss_count", lambda values: int((values > 0).sum())),
            total_gain_events=("reference_relative_gain_count", "sum"),
            total_loss_events=("reference_relative_loss_count", "sum"),
            mean_isoforms=("n_isoforms", "mean"),
            median_long_read_count=("total_long_read_count", "median"),
        )
        .reset_index()
    )
    support.to_csv(args.output_dir / "model_data_support.tsv", sep="\t", index=False)

    sensitivity_rows = []
    for threshold in (5, 10, 20):
        threshold_within_gene_data = within_gene_data[
            (within_gene_data["total_long_read_count"] >= threshold)
            & (within_gene_data["current_missing_locus_count"] == 0)
        ].copy()
        threshold_architecture_data = architecture_data[
            (architecture_data["total_long_read_count"] >= threshold)
            & (architecture_data["current_missing_locus_count"] == 0)
        ].copy()
        try:
            conditional = prepare_conditional_data(threshold_within_gene_data, STRAINS)
            rows, _ = fit_conditional_poisson(
                conditional, f"within_gene_conditional_poisson_support_{threshold}"
            )
            for row in rows:
                if row["term"] == "current_introner_count":
                    row["minimum_long_read_count"] = threshold
                    sensitivity_rows.append(row)
        except Exception as error:
            model_notes.append(f"turnover support threshold {threshold}: FAILED: {error}")
        try:
            architecture = smf.gee(
                "n_extra_isoforms ~ current_introner_count + C(strain) + "
                "log1p_long_read_count + log_gene_span",
                groups="common_gene_id",
                data=threshold_architecture_data,
                family=sm.families.NegativeBinomial(alpha=1.0),
                cov_struct=sm.cov_struct.Exchangeable(),
            ).fit()
            rows = result_rows(
                f"architecture_negative_binomial_support_{threshold}",
                architecture,
                list(architecture.params.index),
                "rate_ratio",
                len(threshold_architecture_data),
                threshold_architecture_data["common_gene_id"].nunique(),
            )
            for row in rows:
                if row["term"] == "current_introner_count":
                    row["minimum_long_read_count"] = threshold
                    sensitivity_rows.append(row)
        except Exception as error:
            model_notes.append(f"architecture support threshold {threshold}: FAILED: {error}")
    pd.DataFrame(sensitivity_rows).to_csv(
        args.output_dir / "read_support_threshold_sensitivity.tsv", sep="\t", index=False
    )

    coefficient_png = args.output_dir / "isoform_turnover_model_coefficients.png"
    overview_png = args.output_dir / "isoform_turnover_overview.png"
    coefficient_plot(coefficients, coefficient_png)
    coefficient_plot(coefficients, coefficient_png.with_suffix(".pdf"))
    overview_plot(supported, coefficients, overview_png)
    overview_plot(supported, coefficients, overview_png.with_suffix(".pdf"))

    key = coefficients[coefficients["term"] == "current_introner_count"].copy()
    lines = [
        "Isoform introner-number model summary",
        "====================================",
        f"Minimum long-read support: {args.minimum_long_read_count}",
        f"Architecture observations: {len(supported)} ({supported['common_gene_id'].nunique()} genes)",
        f"Within-gene eligible observations: {len(within_gene)} ({within_gene['common_gene_id'].nunique()} genes)",
        "",
        "Selected coefficients:",
    ]
    for row in key.itertuples(index=False):
        if row.model not in {
            "within_gene_conditional_poisson_all",
            "within_gene_conditional_logit_all",
            "architecture_negative_binomial_gee",
            "architecture_negative_binomial_exon_adjusted_gee",
            "architecture_any_alternative_isoform_gee",
        }:
            continue
        lines.append(
            f"  {row.model} | {row.term}: {row.effect_ratio:.4f} "
            f"({row.effect_ratio_ci_lower:.4f}-{row.effect_ratio_ci_upper:.4f}), p={row.p_value:.4g}"
        )
    lines.extend(["", "Fit notes:"] + [f"  {note}" for note in model_notes])
    lines.extend(
        [
            "",
            "Interpretation guardrails:",
            "  - The conditional models compare current introner number within the same common gene across strains.",
            "  - The within-gene coefficient is reference-independent and does not assign gain/loss direction.",
            "  - Architecture models summarize current introner burden across genes and strains.",
            "  - Exact splice-chain richness awaits the matching R2C2 transcript GTFs.",
        ]
    )
    (args.output_dir / "model_summary.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

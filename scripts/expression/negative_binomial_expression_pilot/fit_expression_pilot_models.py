#!/usr/bin/env python3
"""Fit paired turnover and within-between NB-GEE expression pilot models."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import minimize
from scipy.special import ndtr


REFERENCE = "CCMP1545"
TURNOVER_VARS = ["reference_relative_gain_count", "reference_relative_loss_count"]
STATE_VARS = ["current_introner_count"]


def add_within_between(data: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    result = data.copy()
    for variable in variables:
        between = result.groupby("gene_id")[variable].transform("mean")
        result[f"{variable}_between"] = between
        result[f"{variable}_within"] = result[variable] - between
    return result


def estimate_within_replicate_alpha(data: pd.DataFrame, offset_column: str) -> float:
    work = data[["gene_id", "strain", "raw_count", offset_column]].copy()
    work["exposure"] = np.exp(work[offset_column])
    totals = work.groupby(["gene_id", "strain"]).agg(
        sum_count=("raw_count", "sum"), sum_exposure=("exposure", "sum")
    )
    totals["rate"] = totals["sum_count"] / totals["sum_exposure"]
    work = work.join(totals["rate"], on=["gene_id", "strain"])
    mu = work["exposure"] * work["rate"]
    alpha = (((work["raw_count"] - mu) ** 2 - work["raw_count"]).sum() / (mu**2).sum())
    return float(max(alpha, 1e-6))


def coefficient_frame(result, model_name: str, family: str) -> pd.DataFrame:
    confidence = result.conf_int()
    rows = []
    for term in result.params.index:
        coefficient = float(result.params[term])
        rows.append(
            {
                "model": model_name,
                "family": family,
                "term": term,
                "coefficient": coefficient,
                "std_error": float(result.bse[term]),
                "statistic": float(result.tvalues[term]),
                "p_value": float(result.pvalues[term]),
                "ci_lower": float(confidence.loc[term, 0]),
                "ci_upper": float(confidence.loc[term, 1]),
                "rate_ratio": float(np.exp(coefficient)),
                "rr_ci_lower": float(np.exp(confidence.loc[term, 0])),
                "rr_ci_upper": float(np.exp(confidence.loc[term, 1])),
            }
        )
    return pd.DataFrame(rows)


class FixedEffectPPMLResult:
    def __init__(self, terms, params, covariance, fittedvalues, converged, iterations):
        self.params = pd.Series(params, index=terms)
        self.bse = pd.Series(np.sqrt(np.diag(covariance)), index=terms)
        self.tvalues = self.params / self.bse
        self.pvalues = pd.Series(2 * (1 - ndtr(np.abs(self.tvalues))), index=terms)
        self._confidence = pd.DataFrame(
            {0: self.params - 1.96 * self.bse, 1: self.params + 1.96 * self.bse}
        )
        self.fittedvalues = pd.Series(fittedvalues)
        self.converged = converged
        self.iterations = iterations

    def conf_int(self):
        return self._confidence


def fit_gene_fe_ppml(
    data: pd.DataFrame,
    variables: list[str],
    offset_column: str,
    model_name: str,
) -> FixedEffectPPMLResult:
    """Poisson pseudo-likelihood with gene fixed effects and gene-clustered SEs."""
    work = data.copy()
    work = work[work.groupby("gene_id")["raw_count"].transform("sum").gt(0)].copy()
    design = pd.DataFrame(
        {
            f"strain_{strain}": work["strain"].eq(strain).astype(float)
            for strain in ["RCC1614", "RCC1749"]
            if work["strain"].eq(strain).any()
        },
        index=work.index,
    )
    for variable in variables:
        design[variable] = work[variable].astype(float)
    design["GC_content"] = work["GC_content"].astype(float)

    x = design.to_numpy(dtype=float)
    y = work["raw_count"].to_numpy(dtype=float)
    offset = work[offset_column].to_numpy(dtype=float)
    group_codes, _ = pd.factorize(work["gene_id"], sort=False)
    n_groups = int(group_codes.max() + 1)
    group_totals = np.bincount(group_codes, weights=y, minlength=n_groups)
    total_count = float(y.sum())

    def probabilities(eta):
        maxima = np.full(n_groups, -np.inf)
        np.maximum.at(maxima, group_codes, eta)
        exponentiated = np.exp(eta - maxima[group_codes])
        denominators = np.bincount(
            group_codes, weights=exponentiated, minlength=n_groups
        )
        probabilities = exponentiated / denominators[group_codes]
        log_denominators = maxima + np.log(denominators)
        return probabilities, log_denominators

    def objective_and_gradient(beta):
        eta = x @ beta + offset
        probabilities_, log_denominators = probabilities(eta)
        log_probabilities = eta - log_denominators[group_codes]
        residual = y - group_totals[group_codes] * probabilities_
        objective = -(y @ log_probabilities) / total_count
        gradient = -(x.T @ residual) / total_count
        return objective, gradient

    optimization = minimize(
        lambda beta: objective_and_gradient(beta),
        np.zeros(x.shape[1]),
        jac=True,
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 500},
    )
    gradient_max = float(np.max(np.abs(optimization.jac)))
    numerically_converged = bool(optimization.success or gradient_max < 1e-7)
    if not numerically_converged:
        raise RuntimeError(f"{model_name} did not converge: {optimization.message}")

    beta = optimization.x
    probabilities_, _ = probabilities(x @ beta + offset)
    fitted = group_totals[group_codes] * probabilities_
    residual = y - fitted

    group_means = np.zeros((n_groups, x.shape[1]))
    np.add.at(group_means, group_codes, probabilities_[:, None] * x)
    centered = x - group_means[group_codes]
    weights = group_totals[group_codes] * probabilities_
    bread = centered.T @ (weights[:, None] * centered)

    group_scores = np.zeros((n_groups, x.shape[1]))
    np.add.at(group_scores, group_codes, residual[:, None] * x)
    meat = group_scores.T @ group_scores
    bread_inverse = np.linalg.pinv(bread)
    covariance = bread_inverse @ meat @ bread_inverse
    correction = (n_groups / (n_groups - 1)) * ((len(work) - 1) / (len(work) - x.shape[1]))
    covariance *= correction

    result = FixedEffectPPMLResult(
        list(design.columns),
        beta,
        covariance,
        pd.Series(fitted, index=work.index),
        numerically_converged,
        optimization.nit,
    )
    result.model_data = work
    return result


def fit_gee(
    data: pd.DataFrame,
    formula: str,
    offset_column: str,
    alpha: float,
    model_name: str,
):
    model = smf.gee(
        formula,
        groups="gene_id",
        data=data,
        family=sm.families.NegativeBinomial(alpha=alpha),
        cov_struct=sm.cov_struct.Exchangeable(),
        offset=data[offset_column],
    )
    result = model.fit(maxiter=100, cov_type="robust")
    if not np.isfinite(result.params).all():
        raise RuntimeError(f"{model_name} produced non-finite coefficients")
    return result


def fit_within_log_sensitivity(
    data: pd.DataFrame,
    introner_variables: list[str],
    model_name: str,
) -> tuple[object, pd.DataFrame]:
    aggregate = (
        data.groupby(["gene_id", "strain"], as_index=False)
        .agg(
            log_expression=("log1p_normalized_count_per_kb", "mean"),
            GC_content=("GC_content", "first"),
            **{variable: (variable, "first") for variable in introner_variables},
        )
    )
    strain = pd.get_dummies(aggregate["strain"], prefix="strain", drop_first=True, dtype=float)
    aggregate = pd.concat([aggregate, strain], axis=1)
    predictors = list(strain.columns) + introner_variables + ["GC_content"]
    values = aggregate[["log_expression", *predictors]]
    demeaned = values - aggregate.groupby("gene_id")[["log_expression", *predictors]].transform(
        "mean"
    )
    nonconstant = demeaned[predictors].std().gt(1e-12)
    predictors = list(nonconstant[nonconstant].index)
    result = sm.OLS(demeaned["log_expression"], demeaned[predictors]).fit(
        cov_type="cluster", cov_kwds={"groups": aggregate["gene_id"]}
    )
    return result, coefficient_frame(result, model_name, "gene-demeaned log expression")


def support_row(name: str, data: pd.DataFrame, alpha: float | None = None, dependence=None):
    unique = data.drop_duplicates(["gene_id", "strain"])
    return {
        "model": name,
        "observations": len(data),
        "genes": data["gene_id"].nunique(),
        "gene_strain_rows": len(unique),
        "gain_positive_gene_strain_rows": int(
            (unique["reference_relative_gain_count"] > 0).sum()
        ),
        "loss_positive_gene_strain_rows": int(
            (unique["reference_relative_loss_count"] > 0).sum()
        ),
        "current_positive_gene_strain_rows": int((unique["current_introner_count"] > 0).sum()),
        "alpha": alpha,
        "exchangeable_dependence": dependence,
    }


def plot_outputs(
    coefficients: pd.DataFrame,
    normalization: pd.DataFrame,
    primary_result,
    output_dir: Path,
) -> None:
    selections = [
        ("paired_turnover_gene_fe_ppml", "reference_relative_gain_count", "Gain, within gene"),
        ("paired_turnover_gene_fe_ppml", "reference_relative_loss_count", "Loss, within gene"),
        ("current_state_gene_fe_ppml", "current_introner_count", "Current burden, within gene"),
        ("current_state_nb_gee", "current_introner_count_between", "Current burden, between genes"),
    ]
    forest_rows = []
    for model, term, label in selections:
        row = coefficients[(coefficients["model"] == model) & (coefficients["term"] == term)].iloc[0].copy()
        row["label"] = label
        forest_rows.append(row)
    forest = pd.DataFrame(forest_rows)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), facecolor="white")
    for axis in axes.flat:
        axis.set_facecolor("white")
    ax = axes[0, 0]
    y = np.arange(len(forest))[::-1]
    colors = ["#2C7FB8", "#2C7FB8", "#F28E2B", "#F28E2B"]
    for position, (_, row), color in zip(y, forest.iterrows(), colors):
        ax.errorbar(
            row["rate_ratio"], position,
            xerr=[[row["rate_ratio"] - row["rr_ci_lower"]], [row["rr_ci_upper"] - row["rate_ratio"]]],
            fmt="o", color=color, capsize=3,
        )
    ax.set_yticks(y, forest["label"])
    ax.axvline(1, color="black", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel("Expression rate ratio")
    ax.set_title("Primary introner associations")

    ax = axes[0, 1]
    ax.scatter(
        normalization["library_size_factor"],
        normalization["median_ratio_size_factor"],
        c=pd.Categorical(normalization["strain"]).codes,
        s=55,
    )
    low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([low, high], [low, high], color="black", linestyle="--", linewidth=1)
    for row in normalization.itertuples(index=False):
        ax.annotate(row.replicate, (row.library_size_factor, row.median_ratio_size_factor), fontsize=7)
    ax.set_xlabel("Total-count library factor")
    ax.set_ylabel("Median-ratio size factor")
    ax.set_title("Normalization comparison")

    primary_data = primary_result.model_data
    fitted = primary_result.fittedvalues
    ax = axes[1, 0]
    ax.hexbin(
        np.log1p(primary_data["raw_count"]),
        np.log1p(fitted),
        gridsize=55,
        bins="log",
        mincnt=1,
        cmap="viridis",
    )
    ax.plot([0, np.log1p(primary_data["raw_count"].max())], [0, np.log1p(primary_data["raw_count"].max())], "--", color="white")
    ax.set_xlabel("log(1 + observed count)")
    ax.set_ylabel("log(1 + fitted count)")
    ax.set_title("Paired-turnover fitted values")

    residual = (primary_data["raw_count"] - fitted) / np.sqrt(fitted)
    ax = axes[1, 1]
    lower, upper = residual.quantile([0.01, 0.99])
    central = residual[(residual >= lower) & (residual <= upper)]
    ax.hist(central, bins=70, color="#5B8DB8", edgecolor="white")
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Pearson residual (central 98%)")
    ax.set_ylabel("Observations")
    ax.set_title("Gene-FE PPML residual distribution")
    ax.text(
        0.98,
        0.95,
        f"1st-99th percentile: {lower:.1f} to {upper:.1f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(
        output_dir / "expression_pilot_diagnostics.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).resolve().parent
    parser.add_argument("--input", type=Path, default=base / "data" / "expression_model_data.tsv")
    parser.add_argument("--normalization", type=Path, default=base / "data" / "normalization_factors.tsv")
    parser.add_argument("--output-dir", type=Path, default=base / "models")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input, sep="\t")
    normalization = pd.read_csv(args.normalization, sep="\t")

    data = add_within_between(data, TURNOVER_VARS + STATE_VARS + ["GC_content"])
    paired = data[data.groupby("gene_id")["strain"].transform("nunique").ge(2)].copy()
    paired_complete = paired[paired["complete_event_callability"]].copy()

    turnover_formula = (
        "raw_count ~ C(strain, Treatment(reference='CCMP1545'))"
        " + reference_relative_gain_count_within"
        " + reference_relative_loss_count_within"
        " + reference_relative_gain_count_between"
        " + reference_relative_loss_count_between"
        " + GC_content_within + GC_content_between"
    )
    state_formula = (
        "raw_count ~ C(strain, Treatment(reference='CCMP1545'))"
        " + current_introner_count_within + current_introner_count_between"
        " + GC_content_within + GC_content_between"
    )

    alpha_median = estimate_within_replicate_alpha(paired_complete, "offset_median_ratio_length")
    alpha_library = estimate_within_replicate_alpha(paired_complete, "offset_library_length")

    turnover_ppml = fit_gene_fe_ppml(
        paired_complete,
        TURNOVER_VARS,
        "offset_median_ratio_length",
        "paired_turnover_gene_fe_ppml",
    )
    current_ppml = fit_gene_fe_ppml(
        paired,
        STATE_VARS,
        "offset_median_ratio_length",
        "current_state_gene_fe_ppml",
    )
    turnover_all_calls_ppml = fit_gene_fe_ppml(
        paired,
        TURNOVER_VARS,
        "offset_median_ratio_length",
        "paired_turnover_all_callability_ppml",
    )
    turnover_library_ppml = fit_gene_fe_ppml(
        paired_complete,
        TURNOVER_VARS,
        "offset_library_length",
        "paired_turnover_library_size_ppml",
    )

    turnover_gee = fit_gee(
        paired_complete,
        turnover_formula,
        "offset_median_ratio_length",
        alpha_median,
        "paired_turnover_nb_gee_sensitivity",
    )
    state = fit_gee(
        paired,
        state_formula,
        "offset_median_ratio_length",
        alpha_median,
        "current_state_nb_gee",
    )
    dispersion_rows = []
    for sensitivity_alpha in [0.02, alpha_median, 0.25, 1.0, 2.65]:
        sensitivity_result = (
            state
            if np.isclose(sensitivity_alpha, alpha_median)
            else fit_gee(
                paired,
                state_formula,
                "offset_median_ratio_length",
                sensitivity_alpha,
                f"current_state_alpha_{sensitivity_alpha:g}",
            )
        )
        confidence = sensitivity_result.conf_int()
        for term in ["current_introner_count_within", "current_introner_count_between"]:
            dispersion_rows.append(
                {
                    "alpha": sensitivity_alpha,
                    "term": term,
                    "coefficient": sensitivity_result.params[term],
                    "std_error": sensitivity_result.bse[term],
                    "p_value": sensitivity_result.pvalues[term],
                    "rate_ratio": np.exp(sensitivity_result.params[term]),
                    "rr_ci_lower": np.exp(confidence.loc[term, 0]),
                    "rr_ci_upper": np.exp(confidence.loc[term, 1]),
                    "exchangeable_dependence": float(sensitivity_result.cov_struct.dep_params),
                }
            )
    dispersion_sensitivity = pd.DataFrame(dispersion_rows)
    dispersion_sensitivity.to_csv(
        output_dir / "architecture_dispersion_sensitivity.tsv", sep="\t", index=False
    )
    _, log_turnover_coef = fit_within_log_sensitivity(
        paired_complete, TURNOVER_VARS, "paired_turnover_within_log_sensitivity"
    )
    _, log_state_coef = fit_within_log_sensitivity(
        paired, STATE_VARS, "current_state_within_log_sensitivity"
    )

    coefficients = pd.concat(
        [
            coefficient_frame(
                turnover_ppml,
                "paired_turnover_gene_fe_ppml",
                "gene-fixed-effect Poisson pseudo-likelihood",
            ),
            coefficient_frame(
                current_ppml,
                "current_state_gene_fe_ppml",
                "gene-fixed-effect Poisson pseudo-likelihood",
            ),
            coefficient_frame(
                turnover_all_calls_ppml,
                "paired_turnover_all_callability_ppml",
                "gene-fixed-effect Poisson pseudo-likelihood",
            ),
            coefficient_frame(
                turnover_library_ppml,
                "paired_turnover_library_size_ppml",
                "gene-fixed-effect Poisson pseudo-likelihood",
            ),
            coefficient_frame(
                turnover_gee,
                "paired_turnover_nb_gee_sensitivity",
                "negative binomial GEE",
            ),
            coefficient_frame(state, "current_state_nb_gee", "negative binomial GEE"),
            log_turnover_coef,
            log_state_coef,
        ],
        ignore_index=True,
    )
    coefficients.to_csv(output_dir / "model_coefficients.tsv", sep="\t", index=False)

    support = pd.DataFrame(
        [
            support_row(
                "paired_turnover_gene_fe_ppml",
                paired_complete,
                None,
                None,
            ),
            support_row(
                "current_state_gene_fe_ppml",
                paired,
                None,
                None,
            ),
            support_row(
                "current_state_nb_gee",
                paired,
                alpha_median,
                float(state.cov_struct.dep_params),
            ),
            support_row(
                "paired_turnover_nb_gee_sensitivity",
                paired_complete,
                alpha_median,
                float(turnover_gee.cov_struct.dep_params),
            ),
            support_row(
                "paired_turnover_all_callability_ppml",
                paired,
                None,
                None,
            ),
            support_row(
                "paired_turnover_library_size_ppml",
                paired_complete,
                None,
                None,
            ),
        ]
    )
    support.to_csv(output_dir / "model_support.tsv", sep="\t", index=False)

    selected = coefficients[
        coefficients["term"].isin(
            [
                "reference_relative_gain_count_within",
                "reference_relative_loss_count_within",
                "current_introner_count_within",
                "current_introner_count_between",
                "reference_relative_gain_count",
                "reference_relative_loss_count",
                "current_introner_count",
            ]
        )
    ].copy()
    selected.to_csv(output_dir / "introner_effect_sensitivity.tsv", sep="\t", index=False)

    plot_outputs(coefficients, normalization, turnover_ppml, output_dir)

    def line(model: str, term: str) -> str:
        row = coefficients[(coefficients["model"] == model) & (coefficients["term"] == term)].iloc[0]
        return (
            f"{model} | {term}: RR {row.rate_ratio:.3f} "
            f"({row.rr_ci_lower:.3f}-{row.rr_ci_upper:.3f}), p={row.p_value:.4g}"
        )

    summary = [
        "Negative-binomial expression pilot summary",
        "============================================",
        f"Within-replicate NB2 alpha (median-ratio offset): {alpha_median:.5f}",
        f"Within-replicate NB2 alpha (library-size offset): {alpha_library:.5f}",
        "",
        "Primary paired turnover effects (gene-FE PPML, gene-clustered SE):",
        "  " + line("paired_turnover_gene_fe_ppml", "reference_relative_gain_count"),
        "  " + line("paired_turnover_gene_fe_ppml", "reference_relative_loss_count"),
        "",
        "Current-state decomposition:",
        "  " + line("current_state_gene_fe_ppml", "current_introner_count"),
        "  " + line("current_state_nb_gee", "current_introner_count_between"),
        "  Architecture alpha sensitivity: between-gene RR range "
        f"{dispersion_sensitivity.loc[dispersion_sensitivity['term'].eq('current_introner_count_between'), 'rate_ratio'].min():.3f}-"
        f"{dispersion_sensitivity.loc[dispersion_sensitivity['term'].eq('current_introner_count_between'), 'rate_ratio'].max():.3f}",
        "",
        "Interpretation guardrails:",
        "  - The paired PPML conditions out each gene's baseline expression.",
        "  - Gene-clustered sandwich errors allow overdispersion in paired PPML.",
        "  - Between effects describe gene architecture and are not turnover effects.",
        "  - Gain/loss are callable CCMP1545-relative differences, not ancestral events.",
        "  - Isoform counts are excluded to avoid expression-dependent detection bias.",
        "  - Effective exonic length is an offset, not an estimated biological effect.",
    ]
    (output_dir / "model_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()

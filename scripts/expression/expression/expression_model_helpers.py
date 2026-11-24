#!/usr/bin/env python3
"""Lightweight count-model helpers for gene-FE PPML and NB-GEE expression models."""

from __future__ import annotations

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
    alpha = ((work["raw_count"] - mu) ** 2 - work["raw_count"]).sum() / (mu**2).sum()
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
        probabilities_ = exponentiated / denominators[group_codes]
        log_denominators = maxima + np.log(denominators)
        return probabilities_, log_denominators

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

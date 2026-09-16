#!/usr/bin/env python3
"""Estimate expression-adjusted retention rates with beta-binomial regression.

The response is boundary-level RNA-seq evidence. For the primary combined
endpoint, retained events are U5 + U3 and spliced events are 2 * S. This
preserves the existing retention PSI exactly:

    (U5 + U3) / (U5 + U3 + 2*S)

Models include intron class, a natural cubic spline for log2 gene CPM, and
replicate fixed effects. Beta-binomial overdispersion handles residual
heterogeneity, while sandwich covariance clustered by host gene accounts for
repeated replicates and multiple introns per gene. Adjusted class rates are
average standardized predictions over the common observed expression and
replicate distribution.

"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
from scipy import optimize, special, stats
from statsmodels.tools.numdiff import approx_hess


CLASSES = ["fixed_present", "polymorphic", "conventional_intron"]
CLASS_LABELS = {
    "fixed_present": "Fixed-present introners",
    "polymorphic": "Polymorphic introners",
    "conventional_intron": "Canonical introns",
}
COLORS = {
    "fixed_present": "#326b77",
    "polymorphic": "#80ae9a",
    "conventional_intron": "#bc272d",
}
PRIMARY_ENDPOINT = "combined_boundary_expression_interaction"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts",
        type=Path,
        default=Path("results/expression/counts/merged_counts_matrix.csv"),
    )
    parser.add_argument(
        "--psi-counts",
        type=Path,
        default=Path(
            "results/expression/functional/"
            "corrected_retention_psi/per_locus_replicate_U5_U3_S.tsv"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/expression/functional/introner_retention_adjusted_rates"),
    )
    parser.add_argument("--min-pooled-junction", type=int, default=5)
    parser.add_argument("--expression-spline-df", type=int, default=3)
    parser.add_argument("--prediction-draws", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=28411)
    return parser.parse_args()


def as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "t", "1", "yes"})


def load_expression(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(path)
    if counts.Geneid.duplicated().any():
        raise ValueError("Geneid is not unique in the merged count matrix")
    replicate_columns = [
        column for column in counts.columns if column not in {"Geneid", "Length"}
    ]
    numeric = counts[replicate_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    library_sizes = numeric.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("At least one featureCounts library has zero total counts")
    long = (
        numeric.assign(Geneid=counts.Geneid)
        .melt(id_vars="Geneid", var_name="replicate", value_name="gene_count")
        .rename(columns={"Geneid": "gene_id"})
    )
    long["library_size"] = long.replicate.map(library_sizes)
    long["gene_cpm"] = 1e6 * long.gene_count / long.library_size
    long["log2_cpm"] = np.log2(long.gene_cpm + 0.5)
    audit = pd.DataFrame(
        {
            "replicate": replicate_columns,
            "library_size": [library_sizes[x] for x in replicate_columns],
        }
    )
    return long, audit


def prepare_data(
    psi_path: Path,
    expression: pd.DataFrame,
    min_pooled_junction: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(psi_path, sep="\t", low_memory=False)
    required = {
        "sample",
        "replicate",
        "gene_id",
        "locus_id",
        "analysis_class",
        "primary_boundary_eligible",
        "U5_donor_crossing",
        "U3_acceptor_crossing",
        "S_junction_signal_exact",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing PSI columns: {sorted(missing)}")
    data = data[data.analysis_class.isin(CLASSES)].copy()
    data["primary_boundary_eligible"] = as_bool(data.primary_boundary_eligible)
    count_columns = [
        "U5_donor_crossing",
        "U3_acceptor_crossing",
        "S_junction_signal_exact",
    ]
    for column in count_columns:
        data[column] = (
            pd.to_numeric(data[column], errors="coerce").fillna(0).astype(np.int64)
        )
        if (data[column] < 0).any():
            raise ValueError(f"Negative count in {column}")

    locus_keys = ["sample", "locus_id", "analysis_class"]
    pooled = (
        data.groupby(locus_keys, observed=True, as_index=False)
        .agg(
            pooled_junction_support=("S_junction_signal_exact", "sum"),
            pooled_U5=("U5_donor_crossing", "sum"),
            pooled_U3=("U3_acceptor_crossing", "sum"),
            boundary_eligible=("primary_boundary_eligible", "all"),
        )
    )
    pooled["pooled_total_boundary_events"] = (
        pooled.pooled_U5 + pooled.pooled_U3 + 2 * pooled.pooled_junction_support
    )
    data = data.merge(
        pooled[
            locus_keys
            + [
                "pooled_junction_support",
                "pooled_total_boundary_events",
                "boundary_eligible",
            ]
        ],
        on=locus_keys,
        validate="many_to_one",
    )
    data["eligible_locus"] = (
        data.boundary_eligible
        & (data.pooled_junction_support >= min_pooled_junction)
    )
    before_expression = len(data)
    data = data.merge(
        expression,
        on=["gene_id", "replicate"],
        how="left",
        validate="many_to_one",
    )
    expression_missing = int(data.log2_cpm.isna().sum())
    data["retained_combined"] = data.U5_donor_crossing + data.U3_acceptor_crossing
    data["spliced_combined"] = 2 * data.S_junction_signal_exact
    data["total_combined"] = data.retained_combined + data.spliced_combined
    data["raw_boundary_retention_psi"] = np.where(
        data.total_combined > 0,
        data.retained_combined / data.total_combined,
        np.nan,
    )
    data["expression_z"] = (
        data.log2_cpm - data.loc[data.eligible_locus, "log2_cpm"].mean()
    ) / data.loc[data.eligible_locus, "log2_cpm"].std(ddof=0)

    audit_rows = [
        {"metric": "input_replicate_rows", "value": before_expression},
        {"metric": "rows_missing_expression", "value": expression_missing},
        {
            "metric": "eligible_loci",
            "value": int(
                pooled[
                    pooled.boundary_eligible
                    & (pooled.pooled_junction_support >= min_pooled_junction)
                ].shape[0]
            ),
        },
    ]
    for class_name in CLASSES:
        subset = pooled[pooled.analysis_class == class_name]
        eligible = subset[
            subset.boundary_eligible
            & (subset.pooled_junction_support >= min_pooled_junction)
        ]
        audit_rows.extend(
            [
                {
                    "metric": f"{class_name}_candidate_loci",
                    "value": len(subset),
                },
                {
                    "metric": f"{class_name}_eligible_loci",
                    "value": len(eligible),
                },
            ]
        )
    return data, pd.DataFrame(audit_rows)


class BetaBinomialRegression:
    """Maximum-likelihood beta-binomial regression with analytic scores."""

    def __init__(self, design: np.ndarray, retained: np.ndarray, spliced: np.ndarray):
        self.x = np.asarray(design, dtype=float)
        self.y = np.asarray(retained, dtype=float)
        self.failure = np.asarray(spliced, dtype=float)
        self.n = self.y + self.failure
        if np.any(self.y < 0) or np.any(self.failure < 0):
            raise ValueError("Counts must be non-negative")
        if np.any(np.mod(self.y, 1) != 0) or np.any(np.mod(self.failure, 1) != 0):
            raise ValueError("Beta-binomial counts must be integers")

    def loglike_score(
        self, parameters: np.ndarray, per_observation: bool = False
    ) -> tuple[float, np.ndarray]:
        beta = parameters[:-1]
        log_kappa = float(parameters[-1])
        kappa = math.exp(log_kappa)
        eta = np.clip(self.x @ beta, -30, 30)
        probability = special.expit(eta)
        alpha = np.clip(probability * kappa, 1e-10, None)
        beta_shape = np.clip((1 - probability) * kappa, 1e-10, None)

        loglik = (
            special.gammaln(self.n + 1)
            - special.gammaln(self.y + 1)
            - special.gammaln(self.failure + 1)
            + special.betaln(self.y + alpha, self.failure + beta_shape)
            - special.betaln(alpha, beta_shape)
        )
        d_probability = kappa * (
            special.digamma(self.y + alpha)
            - special.digamma(self.failure + beta_shape)
            - special.digamma(alpha)
            + special.digamma(beta_shape)
        )
        d_eta = d_probability * probability * (1 - probability)
        d_log_kappa = kappa * (
            probability * special.digamma(self.y + alpha)
            + (1 - probability) * special.digamma(self.failure + beta_shape)
            - special.digamma(self.n + kappa)
            - probability * special.digamma(alpha)
            - (1 - probability) * special.digamma(beta_shape)
            + special.digamma(kappa)
        )
        observation_scores = np.column_stack(
            [self.x * d_eta[:, None], d_log_kappa]
        )
        if per_observation:
            return float(loglik.sum()), observation_scores
        return float(loglik.sum()), observation_scores.sum(axis=0)

    def fit(self) -> optimize.OptimizeResult:
        pooled_probability = np.clip(self.y.sum() / self.n.sum(), 1e-5, 1 - 1e-5)
        initial = np.zeros(self.x.shape[1] + 1)
        initial[0] = special.logit(pooled_probability)
        initial[-1] = math.log(20)

        def objective(par: np.ndarray) -> tuple[float, np.ndarray]:
            loglike, score = self.loglike_score(par)
            return -loglike, -score

        bounds = [(None, None)] * self.x.shape[1] + [(-8, 16)]
        result = optimize.minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": 1500, "ftol": 1e-11, "gtol": 1e-6},
        )
        if not result.success:
            raise RuntimeError(f"Beta-binomial fit failed: {result.message}")
        return result


def clustered_covariance(
    model: BetaBinomialRegression,
    parameters: np.ndarray,
    clusters: pd.Series,
) -> tuple[np.ndarray, np.ndarray]:
    def negative_loglike(par: np.ndarray) -> float:
        return -model.loglike_score(par)[0]

    hessian = approx_hess(parameters, negative_loglike)
    bread = np.linalg.pinv(hessian, rcond=1e-10)
    _, observation_scores = model.loglike_score(parameters, per_observation=True)
    codes, levels = pd.factorize(clusters, sort=False)
    cluster_scores = np.zeros((len(levels), observation_scores.shape[1]))
    np.add.at(cluster_scores, codes, observation_scores)
    meat = cluster_scores.T @ cluster_scores
    n, p, g = len(codes), observation_scores.shape[1], len(levels)
    correction = (g / (g - 1)) * ((n - 1) / (n - p)) if g > 1 and n > p else 1
    covariance = correction * bread @ meat @ bread
    covariance = (covariance + covariance.T) / 2
    return covariance, hessian


def nearest_psd(covariance: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh((covariance + covariance.T) / 2)
    floor = max(float(values.max()) * 1e-10, 1e-12)
    return vectors @ np.diag(np.maximum(values, floor)) @ vectors.T


def standardized_rates(
    beta: np.ndarray,
    beta_covariance: np.ndarray,
    design_info: patsy.DesignInfo,
    prediction_frame: pd.DataFrame,
    draws: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    draw_beta = rng.multivariate_normal(
        beta, nearest_psd(beta_covariance), size=draws, method="eigh"
    )
    rate_draws: dict[str, np.ndarray] = {}
    point_rates: dict[str, float] = {}
    for class_name in CLASSES:
        counterfactual = prediction_frame.copy()
        counterfactual["analysis_class"] = class_name
        x = np.asarray(
            patsy.build_design_matrices(
                [design_info], counterfactual, return_type="dataframe"
            )[0],
            dtype=float,
        )
        point_rates[class_name] = float(special.expit(x @ beta).mean())
        class_draws = np.empty(draws)
        chunk = 100
        for start in range(0, draws, chunk):
            stop = min(start + chunk, draws)
            class_draws[start:stop] = special.expit(
                x @ draw_beta[start:stop].T
            ).mean(axis=0)
        rate_draws[class_name] = class_draws

    rate_rows = []
    for class_name in CLASSES:
        values = rate_draws[class_name]
        rate_rows.append(
            {
                "analysis_class": class_name,
                "class_label": CLASS_LABELS[class_name],
                "adjusted_retention_rate": point_rates[class_name],
                "ci_low": np.quantile(values, 0.025),
                "ci_high": np.quantile(values, 0.975),
                "adjusted_spliced_fraction": 1 - point_rates[class_name],
            }
        )

    contrast_rows = []
    for focal, reference in [
        ("polymorphic", "fixed_present"),
        ("fixed_present", "conventional_intron"),
        ("polymorphic", "conventional_intron"),
    ]:
        difference = rate_draws[focal] - rate_draws[reference]
        ratio = rate_draws[focal] / rate_draws[reference]
        point_difference = point_rates[focal] - point_rates[reference]
        point_ratio = point_rates[focal] / point_rates[reference]
        contrast_rows.append(
            {
                "focal_class": focal,
                "reference_class": reference,
                "adjusted_rate_difference": point_difference,
                "difference_ci_low": np.quantile(difference, 0.025),
                "difference_ci_high": np.quantile(difference, 0.975),
                "adjusted_rate_ratio": point_ratio,
                "ratio_ci_low": np.quantile(ratio, 0.025),
                "ratio_ci_high": np.quantile(ratio, 0.975),
                "two_sided_tail_probability": 2
                * min(float((difference <= 0).mean()), float((difference >= 0).mean())),
            }
        )
    return pd.DataFrame(rate_rows), pd.DataFrame(contrast_rows), rate_draws


def fit_endpoint(
    data: pd.DataFrame,
    endpoint: str,
    retained_column: str,
    spliced_column: str,
    spline_df: int,
    draws: int,
    seed: int,
    class_expression_interaction: bool = False,
) -> dict[str, object]:
    model_data = data[
        data.eligible_locus
        & data.log2_cpm.notna()
        & ((data[retained_column] + data[spliced_column]) > 0)
    ].copy()
    model_data["analysis_class"] = pd.Categorical(
        model_data.analysis_class, categories=CLASSES
    )
    class_term = "C(analysis_class, Treatment(reference='fixed_present'))"
    expression_term = (
        f"cr(expression_z, df={spline_df}, constraints='center')"
    )
    if class_expression_interaction:
        formula = f"{class_term} * {expression_term} + C(replicate)"
    else:
        formula = f"{class_term} + {expression_term} + C(replicate)"
    design = patsy.dmatrix(formula, model_data, return_type="dataframe")
    model = BetaBinomialRegression(
        design,
        model_data[retained_column].to_numpy(),
        model_data[spliced_column].to_numpy(),
    )
    result = model.fit()
    covariance, hessian = clustered_covariance(
        model, result.x, model_data.gene_id
    )
    standard_errors = np.sqrt(np.clip(np.diag(covariance), 0, None))
    names = list(design.columns) + ["log_kappa"]
    z_values = result.x / standard_errors
    coefficients = pd.DataFrame(
        {
            "endpoint": endpoint,
            "term": names,
            "estimate": result.x,
            "cluster_robust_se": standard_errors,
            "z": z_values,
            "p_value": 2 * stats.norm.sf(np.abs(z_values)),
            "ci_low": result.x - 1.96 * standard_errors,
            "ci_high": result.x + 1.96 * standard_errors,
        }
    )
    rates, contrasts, rate_draws = standardized_rates(
        result.x[:-1],
        covariance[:-1, :-1],
        design.design_info,
        model_data,
        draws,
        seed,
    )
    _, final_score = model.loglike_score(result.x)
    hessian_eigenvalues = np.linalg.eigvalsh((hessian + hessian.T) / 2)
    diagnostics = pd.DataFrame(
        [
            {
                "endpoint": endpoint,
                "converged": bool(result.success),
                "optimizer_message": str(result.message),
                "n_observations": len(model_data),
                "n_host_gene_clusters": model_data.gene_id.nunique(),
                "n_parameters": len(result.x),
                "design_matrix_rank": int(np.linalg.matrix_rank(model.x)),
                "log_likelihood": -float(result.fun),
                "max_absolute_score": float(np.abs(final_score).max()),
                "minimum_hessian_eigenvalue": float(hessian_eigenvalues.min()),
                "maximum_hessian_eigenvalue": float(hessian_eigenvalues.max()),
                "hessian_condition_number": float(
                    hessian_eigenvalues.max() / hessian_eigenvalues.min()
                ),
                "kappa": math.exp(float(result.x[-1])),
                "overdispersion_rho": 1 / (math.exp(float(result.x[-1])) + 1),
            }
        ]
    )
    rates.insert(0, "endpoint", endpoint)
    contrasts.insert(0, "endpoint", endpoint)
    return {
        "data": model_data,
        "formula": formula,
        "design": design,
        "result": result,
        "covariance": covariance,
        "hessian": hessian,
        "coefficients": coefficients,
        "rates": rates,
        "contrasts": contrasts,
        "diagnostics": diagnostics,
        "rate_draws": rate_draws,
        "kappa": math.exp(float(result.x[-1])),
        "rho": 1 / (math.exp(float(result.x[-1])) + 1),
    }


def plot_rates(rates: pd.DataFrame, path: Path) -> None:
    primary = rates[rates.endpoint == PRIMARY_ENDPOINT].copy()
    primary["analysis_class"] = pd.Categorical(
        primary.analysis_class, categories=CLASSES, ordered=True
    )
    primary = primary.sort_values("analysis_class")
    figure, axis = plt.subplots(figsize=(6.5, 4.4))
    x = np.arange(len(primary))
    y = primary.adjusted_retention_rate.to_numpy()
    lower = y - primary.ci_low.to_numpy()
    upper = primary.ci_high.to_numpy() - y
    axis.errorbar(
        x,
        y,
        yerr=np.vstack([lower, upper]),
        fmt="none",
        ecolor="#333333",
        elinewidth=1.4,
        capsize=4,
        zorder=2,
    )
    for position, row in zip(x, primary.itertuples()):
        axis.scatter(
            position,
            row.adjusted_retention_rate,
            s=95,
            color=COLORS[row.analysis_class],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
    axis.set_xticks(x, primary.class_label, rotation=12, ha="right")
    axis.set_ylabel("Adjusted intron-retention rate")
    axis.set_ylim(bottom=0)
    axis.grid(axis="y", color="#d8d8d8", linewidth=0.7, alpha=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_title("Expression-adjusted retention rates")
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)


def write_summary(
    cfg: argparse.Namespace,
    fits: dict[str, dict[str, object]],
    rates: pd.DataFrame,
    contrasts: pd.DataFrame,
    path: Path,
) -> None:
    additive = fits["combined_boundary"]
    primary = fits[PRIMARY_ENDPOINT]
    likelihood_ratio = 2 * (
        -float(primary["result"].fun) + float(additive["result"].fun)
    )
    likelihood_ratio_df = (
        len(primary["result"].x) - len(additive["result"].x)
    )
    likelihood_ratio_p = stats.chi2.sf(likelihood_ratio, likelihood_ratio_df)
    lines = [
        "Expression-adjusted intron-retention beta-binomial analysis",
        "",
        "Primary response: retained boundary events = U5 + U3; "
        "spliced boundary events = 2*S.",
        f"Locus eligibility: primary_boundary_eligible and pooled exact "
        f"junction support >= {cfg.min_pooled_junction}.",
        "Covariates: intron class, a natural cubic spline of log2(CPM + 0.5), "
        "their interaction, and RNA-seq replicate fixed effects.",
        "Uncertainty: sandwich covariance clustered by host gene.",
        "Adjusted rates: average standardized predictions over the common "
        "observed expression and replicate distribution.",
        "",
        f"Primary observations: {len(primary['data'])}",
        f"Host-gene clusters: {primary['data'].gene_id.nunique()}",
        f"Log likelihood: {primary['result'].fun * -1:.6f}",
        f"Beta-binomial kappa: {primary['kappa']:.6f}",
        f"Intraclass overdispersion rho: {primary['rho']:.6f}",
        f"Optimizer convergence: {primary['result'].message}",
        f"Class-by-expression interaction likelihood-ratio test: "
        f"chi-square={likelihood_ratio:.6f}, df={likelihood_ratio_df}, "
        f"p={likelihood_ratio_p:.6g}",
        "",
        "Adjusted class rates",
        rates[rates.endpoint == PRIMARY_ENDPOINT].to_string(index=False),
        "",
        "Adjusted contrasts",
        contrasts[contrasts.endpoint == PRIMARY_ENDPOINT].to_string(index=False),
        "",
        "Sensitivity endpoints",
        rates[rates.endpoint != PRIMARY_ENDPOINT].to_string(index=False),
        "",
        "Interpretive limitation: introner boundary resolution requires exact "
        "junction evidence, so adjusted rates apply to junction-resolved introners "
        "and should not be interpreted as unconditional rates over all genomic "
        "introner loci.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    cfg = arguments()
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    expression, library_audit = load_expression(cfg.counts)
    data, selection_audit = prepare_data(
        cfg.psi_counts, expression, cfg.min_pooled_junction
    )
    endpoints = {
        "combined_boundary": (
            "retained_combined",
            "spliced_combined",
            False,
        ),
        "donor_boundary": (
            "U5_donor_crossing",
            "S_junction_signal_exact",
            False,
        ),
        "acceptor_boundary": (
            "U3_acceptor_crossing",
            "S_junction_signal_exact",
            False,
        ),
        "combined_boundary_expression_interaction": (
            "retained_combined",
            "spliced_combined",
            True,
        ),
    }
    fits = {}
    for index, (endpoint, columns) in enumerate(endpoints.items()):
        print(f"Fitting {endpoint}", flush=True)
        fits[endpoint] = fit_endpoint(
            data,
            endpoint,
            columns[0],
            columns[1],
            cfg.expression_spline_df,
            cfg.prediction_draws,
            cfg.seed + index,
            columns[2],
        )
    if cfg.min_pooled_junction < 10:
        sensitivity = data.copy()
        sensitivity["eligible_locus"] = (
            sensitivity.boundary_eligible
            & (sensitivity.pooled_junction_support >= 10)
        )
        endpoint = "combined_boundary_min_junction_10"
        print(f"Fitting {endpoint}", flush=True)
        fits[endpoint] = fit_endpoint(
            sensitivity,
            endpoint,
            "retained_combined",
            "spliced_combined",
            cfg.expression_spline_df,
            cfg.prediction_draws,
            cfg.seed + len(fits),
            False,
        )

    rates = pd.concat([fit["rates"] for fit in fits.values()], ignore_index=True)
    contrasts = pd.concat(
        [fit["contrasts"] for fit in fits.values()], ignore_index=True
    )
    coefficients = pd.concat(
        [fit["coefficients"] for fit in fits.values()], ignore_index=True
    )
    diagnostics = pd.concat(
        [fit["diagnostics"] for fit in fits.values()], ignore_index=True
    )
    model_input = fits[PRIMARY_ENDPOINT]["data"].copy()
    keep = [
        "sample",
        "replicate",
        "gene_id",
        "locus_id",
        "analysis_class",
        "U5_donor_crossing",
        "U3_acceptor_crossing",
        "S_junction_signal_exact",
        "retained_combined",
        "spliced_combined",
        "total_combined",
        "raw_boundary_retention_psi",
        "pooled_junction_support",
        "gene_count",
        "library_size",
        "gene_cpm",
        "log2_cpm",
        "expression_z",
    ]
    model_input[keep].to_csv(
        cfg.outdir / "model_input.tsv", sep="\t", index=False
    )
    rates.to_csv(cfg.outdir / "adjusted_class_retention_rates.tsv", sep="\t", index=False)
    contrasts.to_csv(
        cfg.outdir / "adjusted_class_retention_contrasts.tsv", sep="\t", index=False
    )
    coefficients.to_csv(
        cfg.outdir / "model_coefficients.tsv", sep="\t", index=False
    )
    diagnostics.to_csv(
        cfg.outdir / "model_diagnostics.tsv", sep="\t", index=False
    )
    selection_audit.to_csv(
        cfg.outdir / "selection_audit.tsv", sep="\t", index=False
    )
    library_audit.to_csv(
        cfg.outdir / "expression_library_audit.tsv", sep="\t", index=False
    )
    plot_rates(rates, cfg.outdir / "adjusted_class_retention_rates.png")
    write_summary(
        cfg,
        fits,
        rates,
        contrasts,
        cfg.outdir / "model_summary.txt",
    )
    print(f"Wrote {cfg.outdir}", flush=True)


if __name__ == "__main__":
    main()

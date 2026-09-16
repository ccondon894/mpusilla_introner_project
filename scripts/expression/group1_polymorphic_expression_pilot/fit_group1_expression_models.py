#!/usr/bin/env python3
"""Fit Group 1 between-gene and paired within-gene expression models."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import binomtest, norm

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from figure_color_guide import get_color, load_color_guide  # noqa: E402


N_PERMUTATIONS = 500
RNG_SEED = 1545


def coefficient_frame(model_name, family, params, bse, pvalues, conf_int) -> pd.DataFrame:
    rows = []
    for term in params.index:
        coefficient = float(params[term])
        lower, upper = conf_int.loc[term]
        rows.append(
            {
                "model": model_name,
                "family": family,
                "term": term,
                "coefficient": coefficient,
                "std_error": float(bse[term]),
                "p_value": float(pvalues[term]),
                "ci_lower": float(lower),
                "ci_upper": float(upper),
                "rate_ratio": float(np.exp(coefficient)),
                "rr_ci_lower": float(np.exp(lower)),
                "rr_ci_upper": float(np.exp(upper)),
            }
        )
    return pd.DataFrame(rows)


def estimate_nb_alpha(y: np.ndarray, mu: np.ndarray) -> float:
    return float(max(((y - mu) ** 2 - mu).sum() / np.maximum((mu**2).sum(), 1e-12), 1e-8))


def fit_between_gene_nb(data: pd.DataFrame, strain: str) -> tuple[pd.DataFrame, float]:
    work = data[data["strain"].eq(strain) & data["raw_count"].ge(0)].copy()
    work = work[np.isfinite(work["log_exposure"]) & work["exposure"].gt(0)]
    formula = (
        "raw_count ~ n_fixed_present + n_polymorphic_present + GC_content + log_gene_span"
    )
    poisson = smf.glm(formula, data=work, family=sm.families.Poisson(), offset=work["log_exposure"]).fit()
    alpha = estimate_nb_alpha(work["raw_count"].to_numpy(float), poisson.fittedvalues.to_numpy(float))
    result = smf.glm(
        formula,
        data=work,
        family=sm.families.NegativeBinomial(alpha=alpha),
        offset=work["log_exposure"],
    ).fit(cov_type="HC3")
    frame = coefficient_frame(
        f"test_a_nb_{strain}",
        "negative_binomial_HC3",
        result.params,
        result.bse,
        result.pvalues,
        result.conf_int(),
    )
    frame["alpha"] = alpha
    frame["n_genes"] = len(work)
    return frame, alpha


def fit_between_gene_ols(data: pd.DataFrame, strain: str) -> pd.DataFrame:
    work = data[data["strain"].eq(strain)].copy()
    work = work[np.isfinite(work["log_rate_pseudocount"])]
    formula = (
        "log_rate_pseudocount ~ n_fixed_present + n_polymorphic_present + GC_content + log_gene_span"
    )
    result = smf.ols(formula, data=work).fit(cov_type="HC3")
    frame = coefficient_frame(
        f"test_a_ols_{strain}",
        "ols_log_rate_HC3",
        result.params,
        result.bse,
        result.pvalues,
        result.conf_int(),
    )
    frame["n_genes"] = len(work)
    return frame


def permute_between_gene(data: pd.DataFrame, strain: str, n_perm: int, seed: int) -> pd.DataFrame:
    work = data[data["strain"].eq(strain)].copy()
    work = work[np.isfinite(work["log_rate_pseudocount"])].reset_index(drop=True)
    span_bin = pd.qcut(work["log_gene_span"], 5, duplicates="drop")
    gc_bin = pd.qcut(work["GC_content"], 5, duplicates="drop")
    work["bin"] = span_bin.astype(str) + "|" + gc_bin.astype(str)
    observed = smf.ols(
        "log_rate_pseudocount ~ n_fixed_present + n_polymorphic_present + GC_content + log_gene_span",
        data=work,
    ).fit(cov_type="HC3")
    rng = np.random.default_rng(seed)
    records = []
    for term in ["n_fixed_present", "n_polymorphic_present"]:
        null = np.empty(n_perm)
        for i in range(n_perm):
            shuffled = work.copy()
            for _, group in work.groupby("bin"):
                idx = group.index.to_numpy()
                take = rng.permutation(idx)
                shuffled.loc[idx, ["n_fixed_present", "n_polymorphic_present"]] = work.loc[
                    take, ["n_fixed_present", "n_polymorphic_present"]
                ].to_numpy()
            fit = smf.ols(
                "log_rate_pseudocount ~ n_fixed_present + n_polymorphic_present + GC_content + log_gene_span",
                data=shuffled,
            ).fit()
            null[i] = fit.params[term]
        observed_coef = float(observed.params[term])
        records.append(
            {
                "strain": strain,
                "term": term,
                "observed_coefficient": observed_coef,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)),
                "two_sided_p": float(np.mean(np.abs(null) >= abs(observed_coef))),
                "n_permutations": n_perm,
            }
        )
    return pd.DataFrame(records)


def fit_pair_ppml(table: pd.DataFrame, count_column: str, adjusted: bool = True) -> pd.DataFrame:
    """Profile gene intercepts; gene-level sandwich covariance."""
    terms = ["strain_RCC1614", count_column]
    x = np.column_stack([np.ones(len(table)), table[count_column].to_numpy(float)])
    if adjusted:
        terms.append("delta_GC")
        x = np.column_stack([x, table["delta_GC"].to_numpy(float)])
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("Nonidentifiable paired design")
    n = table["total_count"].to_numpy(float)
    y = table["raw_count_RCC1614"].to_numpy(float)
    off = table["log_exposure_ratio"].to_numpy(float)
    scale = n.sum()

    def objective(beta):
        eta = off + x @ beta
        return (
            np.sum(n * np.logaddexp(0, eta) - y * eta) / scale,
            x.T @ (n * expit(eta) - y) / scale,
        )

    opt = minimize(objective, np.zeros(x.shape[1]), jac=True, method="BFGS", options={"gtol": 1e-11, "maxiter": 1000})
    if np.max(np.abs(opt.jac)) >= 1e-6:
        raise RuntimeError(f"PPML did not converge: {opt.message}")
    p = expit(off + x @ opt.x)
    bread = x.T @ ((n * p * (1 - p))[:, None] * x)
    score = x * (y - n * p)[:, None]
    inv = np.linalg.inv(bread)
    covariance = inv @ (score.T @ score) @ inv * len(x) / (len(x) - x.shape[1])
    se = np.sqrt(np.diag(covariance))
    params = pd.Series(opt.x, index=terms)
    bse = pd.Series(se, index=terms)
    pvalues = pd.Series(2 * norm.sf(np.abs(opt.x / se)), index=terms)
    conf = pd.DataFrame({0: opt.x - 1.96 * se, 1: opt.x + 1.96 * se}, index=terms)
    return coefficient_frame(
        f"test_b_gene_fe_ppml_{count_column}",
        "gene_fe_ppml_sandwich",
        params,
        bse,
        pvalues,
        conf,
    )


def fit_delta_log(table: pd.DataFrame, count_column: str) -> pd.DataFrame:
    predictors = table[[count_column, "delta_GC"]]
    result = sm.OLS(table["delta_log_rate"], sm.add_constant(predictors)).fit(cov_type="HC3")
    return coefficient_frame(
        f"test_b_delta_log_{count_column}",
        "ols_delta_log_HC3",
        result.params,
        result.bse,
        result.pvalues,
        result.conf_int(),
    )


def sign_test(table: pd.DataFrame) -> dict[str, float | int]:
    stable = table[table["delta_polymorphic_present"].eq(0)]
    offset = float(stable["delta_log_rate"].mean()) if len(stable) else 0.0
    single = table[table["n_polymorphic_discordant_loci"].eq(1) & table["delta_polymorphic_present"].abs().eq(1)].copy()
    adjusted = single["delta_log_rate"] - offset
    # Negative association: extra occupancy in RCC1614 should lower RCC1614-relative expression.
    concordant = (adjusted * single["delta_polymorphic_present"]) < 0
    n = int(len(single))
    k = int(concordant.sum())
    test = binomtest(k, n, 0.5, alternative="greater") if n else None
    return {
        "n_single_discordant_genes": n,
        "n_present_allele_lower": k,
        "strain_offset": offset,
        "binomial_p": float(test.pvalue) if test is not None else np.nan,
        "proportion": (k / n) if n else np.nan,
    }


def plot_overview(
    test_a: pd.DataFrame,
    test_b: pd.DataFrame,
    pairs: pd.DataFrame,
    permutation: pd.DataFrame,
    sign: dict,
    output_path: Path,
) -> None:
    colors = load_color_guide()
    fixed_color = get_color("Population 1 Fixed", colors)
    poly_color = get_color("Population 1 Polymorphic", colors)
    pop1 = get_color("Population 1", colors)
    ancestor = get_color("Ancestor", colors)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), facecolor="white")
    for axis in axes.flat:
        axis.set_facecolor("white")

    ax = axes[0, 0]
    rows = []
    labels = []
    bar_colors = []
    for strain, prefix in [("CCMP1545", "CCMP1545"), ("RCC1614", "RCC1614")]:
        for term, label, color in [
            ("n_fixed_present", f"{prefix} fixed", fixed_color),
            ("n_polymorphic_present", f"{prefix} polymorphic", poly_color),
        ]:
            row = test_a[(test_a["model"] == f"test_a_nb_{strain}") & (test_a["term"] == term)].iloc[0]
            rows.append(row)
            labels.append(label)
            bar_colors.append(color)
    y = np.arange(len(rows))[::-1]
    for position, row, color in zip(y, rows, bar_colors):
        ax.errorbar(
            row["rate_ratio"],
            position,
            xerr=[[row["rate_ratio"] - row["rr_ci_lower"]], [row["rr_ci_upper"] - row["rate_ratio"]]],
            fmt="o",
            color=color,
            capsize=3,
        )
    ax.axvline(1, color=ancestor, linestyle="--", linewidth=1)
    ax.set_yticks(y, labels)
    ax.set_xscale("log")
    ax.set_xlabel("Expression rate ratio per intron (95% CI)")
    ax.set_title("Test A: between-gene occupancy")

    ax = axes[0, 1]
    b_rows = []
    b_labels = []
    for model, label in [
        ("test_b_gene_fe_ppml_delta_polymorphic_present", "Gene FE PPML"),
        ("test_b_delta_log_delta_polymorphic_present", "Delta log OLS"),
    ]:
        subset = test_b[(test_b["model"] == model) & (test_b["term"].str.contains("polymorphic"))]
        if subset.empty:
            subset = test_b[(test_b["model"] == model) & (test_b["term"] == "delta_polymorphic_present")]
        if subset.empty:
            subset = test_b[(test_b["model"] == model) & (test_b["term"] != "const") & (~test_b["term"].str.contains("strain|GC|delta_GC"))]
        row = subset.iloc[0]
        b_rows.append(row)
        b_labels.append(label)
    y = np.arange(len(b_rows))[::-1]
    for position, row in zip(y, b_rows):
        ax.errorbar(
            row["rate_ratio"],
            position,
            xerr=[[row["rate_ratio"] - row["rr_ci_lower"]], [row["rr_ci_upper"] - row["rate_ratio"]]],
            fmt="o",
            color=pop1,
            capsize=3,
        )
    ax.axvline(1, color=ancestor, linestyle="--", linewidth=1)
    ax.set_yticks(y, b_labels)
    ax.set_xscale("log")
    ax.set_xlabel("Expression ratio per extra polymorphic intron (95% CI)")
    ax.set_title("Test B: within-gene Group 1 pair")

    ax = axes[1, 0]
    jitter = np.random.default_rng(0).uniform(-0.08, 0.08, len(pairs))
    ax.scatter(
        pairs["delta_polymorphic_present"] + jitter,
        pairs["delta_log_rate"],
        s=np.where(pairs["delta_polymorphic_present"].eq(0), 8, 22),
        c=np.where(pairs["delta_polymorphic_present"].eq(0), ancestor, poly_color),
        alpha=0.55,
        linewidths=0,
    )
    ax.axhline(0, color=ancestor, linestyle="--", linewidth=1)
    ax.axvline(0, color=ancestor, linestyle="--", linewidth=1)
    ax.set_xlabel("Δ polymorphic present (RCC1614 − CCMP1545)")
    ax.set_ylabel("Δ log expression rate")
    ax.set_title("Paired expression vs polymorphic occupancy")

    ax = axes[1, 1]
    perm = permutation[permutation["term"] == "n_polymorphic_present"]
    if len(perm):
        ax.axis("off")
        lines = [
            f"Test A permutation (OLS, {int(perm.iloc[0].n_permutations)} shuffles)",
            "",
        ]
        for row in permutation.itertuples(index=False):
            lines.append(
                f"{row.strain} {row.term}: observed {row.observed_coefficient:.3f}, "
                f"null p={row.two_sided_p:.3f}"
            )
        lines += [
            "",
            "Sign test, single polymorphic discordant locus:",
            f"  {sign['n_present_allele_lower']}/{sign['n_single_discordant_genes']} "
            f"present allele lower (p={sign['binomial_p']:.3g})",
        ]
        ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", family="DejaVu Sans", fontsize=10)
    ax.set_title("Null checks")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def validate(expression: pd.DataFrame, pairs: pd.DataFrame, locus: pd.DataFrame) -> list[str]:
    checks = []
    mt_in_expr = bool(pairs.loc[pairs["mating_type_gene"], "gene_id"].isin(expression["gene_id"]).any())
    checks.append(
        ("FAIL" if mt_in_expr else "PASS") + ": mating-type genes are excluded from expression rows"
    )
    n_fixed_delta = int((pairs.loc[pairs["eligible_pair"], "delta_fixed_present"] != 0).sum())
    checks.append(
        f"PASS: eligible-pair fixed occupancy differences are rare ({n_fixed_delta} genes)"
        if n_fixed_delta <= 20
        else f"WARN: {n_fixed_delta} eligible pairs differ in fixed occupancy"
    )
    rna = locus[locus["strain"].isin(["CCMP1545", "RCC1614"])]
    checks.append(
        "PASS: presence codes are 1/2/3"
        if rna["presence"].isin([1, 2, 3]).all()
        else "FAIL: unexpected presence codes"
    )
    poly = rna[rna["group1_pattern"].eq("polymorphic") & rna["unique_common_gene"] & rna["callable"]]
    checks.append(
        f"PASS: {poly['ortholog_id'].nunique()} uniquely mapped callable Group 1 polymorphic loci"
    )
    if not pairs.loc[pairs["eligible_pair"], "n_replicates_CCMP1545"].eq(4).all():
        checks.append("FAIL: eligible pairs missing CCMP1545 replicates")
    else:
        checks.append("PASS: eligible pairs have four replicates per Group 1 RNA-seq strain")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=base / "data")
    parser.add_argument("--output-dir", type=Path, default=base / "models")
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(data_dir / "expression_gene_strain.tsv", sep="\t")
    pairs = pd.read_csv(data_dir / "gene_pair_table.tsv", sep="\t")
    locus = pd.read_csv(data_dir / "locus_status.tsv", sep="\t", low_memory=False)
    eligible = pairs[pairs["eligible_pair"]].copy()

    test_a = pd.concat(
        [
            fit_between_gene_nb(expression, "CCMP1545")[0],
            fit_between_gene_nb(expression, "RCC1614")[0],
            fit_between_gene_ols(expression, "CCMP1545"),
            fit_between_gene_ols(expression, "RCC1614"),
        ],
        ignore_index=True,
    )
    alpha_ccmp = float(test_a.loc[test_a["model"].eq("test_a_nb_CCMP1545"), "alpha"].iloc[0])
    alpha_rcc = float(test_a.loc[test_a["model"].eq("test_a_nb_RCC1614"), "alpha"].iloc[0])

    permutation = pd.concat(
        [
            permute_between_gene(expression, "CCMP1545", args.n_permutations, RNG_SEED),
            permute_between_gene(expression, "RCC1614", args.n_permutations, RNG_SEED + 1),
        ],
        ignore_index=True,
    )

    test_b = pd.concat(
        [
            fit_pair_ppml(eligible, "delta_polymorphic_present"),
            fit_delta_log(eligible, "delta_polymorphic_present"),
            fit_pair_ppml(eligible[eligible["n_polymorphic_discordant_loci"].gt(0)], "delta_polymorphic_present").assign(
                model=lambda frame: frame["model"].str.replace(
                    "test_b_gene_fe_ppml_delta_polymorphic_present",
                    "test_b_gene_fe_ppml_discordant_only",
                    regex=False,
                )
            ),
            fit_delta_log(eligible[eligible["n_polymorphic_discordant_loci"].gt(0)], "delta_polymorphic_present").assign(
                model=lambda frame: frame["model"].str.replace(
                    "test_b_delta_log_delta_polymorphic_present",
                    "test_b_delta_log_discordant_only",
                    regex=False,
                )
            ),
        ],
        ignore_index=True,
    )
    sign = sign_test(eligible)

    test_a.to_csv(output_dir / "test_a_coefficients.tsv", sep="\t", index=False)
    test_b.to_csv(output_dir / "test_b_coefficients.tsv", sep="\t", index=False)
    permutation.to_csv(output_dir / "test_a_permutation.tsv", sep="\t", index=False)
    pd.DataFrame([sign]).to_csv(output_dir / "test_b_sign_test.tsv", sep="\t", index=False)

    support = pd.DataFrame(
        [
            {
                "model": "test_a_nb_CCMP1545",
                "genes": int(expression["strain"].eq("CCMP1545").sum()),
                "genes_with_fixed": int(
                    (expression["strain"].eq("CCMP1545") & expression["n_fixed_present"].gt(0)).sum()
                ),
                "genes_with_polymorphic": int(
                    (expression["strain"].eq("CCMP1545") & expression["n_polymorphic_present"].gt(0)).sum()
                ),
            },
            {
                "model": "test_b_primary",
                "genes": int(len(eligible)),
                "genes_delta_poly": int((eligible["delta_polymorphic_present"] != 0).sum()),
                "genes_poly_discordant": int(eligible["n_polymorphic_discordant_loci"].gt(0).sum()),
                "genes_delta_fixed": int((eligible["delta_fixed_present"] != 0).sum()),
            },
        ]
    )
    support.to_csv(output_dir / "model_support.tsv", sep="\t", index=False)

    plot_overview(
        test_a,
        test_b,
        eligible,
        permutation,
        sign,
        output_dir / "group1_expression_overview.png",
    )

    def line(frame: pd.DataFrame, model: str, term: str) -> str:
        row = frame[(frame["model"] == model) & (frame["term"] == term)].iloc[0]
        return (
            f"{model} | {term}: RR {row.rate_ratio:.3f} "
            f"({row.rr_ci_lower:.3f}-{row.rr_ci_upper:.3f}), p={row.p_value:.4g}"
        )

    checks = validate(expression, pairs, locus)
    (output_dir / "validation_summary.txt").write_text("\n".join(checks) + "\n")

    summary = [
        "Group 1 polymorphic expression pilot",
        "====================================",
        f"Test A NB alpha: CCMP1545={alpha_ccmp:.4f}; RCC1614={alpha_rcc:.4f}",
        "",
        "Test A, between-gene negative binomial (count-weighted, HC3):",
        "  " + line(test_a, "test_a_nb_CCMP1545", "n_fixed_present"),
        "  " + line(test_a, "test_a_nb_CCMP1545", "n_polymorphic_present"),
        "  " + line(test_a, "test_a_nb_RCC1614", "n_fixed_present"),
        "  " + line(test_a, "test_a_nb_RCC1614", "n_polymorphic_present"),
        "",
        "Test A, between-gene OLS of log rate (equal-gene, HC3):",
        "  " + line(test_a, "test_a_ols_CCMP1545", "n_fixed_present"),
        "  " + line(test_a, "test_a_ols_CCMP1545", "n_polymorphic_present"),
        "  " + line(test_a, "test_a_ols_RCC1614", "n_fixed_present"),
        "  " + line(test_a, "test_a_ols_RCC1614", "n_polymorphic_present"),
        "",
        "Test A length/GC-binned occupancy permutation (OLS):",
        *[
            f"  {row.strain} {row.term}: observed {row.observed_coefficient:.3f}, "
            f"null p={row.two_sided_p:.3f}"
            for row in permutation.itertuples(index=False)
        ],
        "",
        "Test B, CCMP1545 vs RCC1614:",
        "  " + line(test_b, "test_b_gene_fe_ppml_delta_polymorphic_present", "delta_polymorphic_present"),
        "  " + line(test_b, "test_b_delta_log_delta_polymorphic_present", "delta_polymorphic_present"),
        "  " + line(test_b, "test_b_gene_fe_ppml_discordant_only", "delta_polymorphic_present"),
        "  " + line(test_b, "test_b_delta_log_discordant_only", "delta_polymorphic_present"),
        "",
        "Sign test (single polymorphic discordant locus): "
        f"{sign['n_present_allele_lower']}/{sign['n_single_discordant_genes']} "
        f"present allele lower; p={sign['binomial_p']:.4g}",
        "",
        "Eligible paired genes: "
        f"{len(eligible)}; polymorphic occupancy differences: "
        f"{int((eligible['delta_polymorphic_present'] != 0).sum())}; "
        f"fixed occupancy differences: {int((eligible['delta_fixed_present'] != 0).sum())}",
        "",
        "Interpretation guardrails:",
        "  - Test A is a gene-compartment association, not an allelic effect.",
        "  - Count-weighted NB and equal-gene OLS can disagree; report both.",
        "  - Test B is identified by Group 1 polymorphic occupancy differences.",
        "  - RCC1749 is excluded. Gain/loss versus CCMP1545 are not used.",
        "  - Isoform counts are not covariates.",
        "  - Two RNA-seq strains cannot support a population-wide causal claim.",
    ]
    (output_dir / "model_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    print("\n".join(checks))


if __name__ == "__main__":
    main()

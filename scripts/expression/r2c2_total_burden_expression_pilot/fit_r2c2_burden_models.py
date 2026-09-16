#!/usr/bin/env python3
"""Fit R2C2 between-gene burden and CCMP1545-relative paired occupancy models."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")
os.environ.setdefault("TMPDIR", "/scratch1/chris/tmp")

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
STRAINS = ("CCMP1545", "RCC1614", "RCC1749")
OUTCOMES = ("all", "fsm", "major")


def coefficient_frame(model_name, family, params, bse, pvalues, conf_int, **extra) -> pd.DataFrame:
    rows = []
    for term in params.index:
        coefficient = float(params[term])
        lower, upper = conf_int.loc[term]
        row = {
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
        row.update(extra)
        rows.append(row)
    return pd.DataFrame(rows)


def estimate_nb_alpha(y: np.ndarray, mu: np.ndarray) -> float:
    return float(max(((y - mu) ** 2 - mu).sum() / np.maximum((mu**2).sum(), 1e-12), 1e-8))


def fit_between_gene(data: pd.DataFrame, strain: str, outcome: str, predictors: str, model_prefix: str) -> pd.DataFrame:
    work = data[data["strain"].eq(strain)].copy()
    work["raw_count"] = work[f"count_{outcome}"]
    work = work[np.isfinite(work["log_exposure"]) & work["exposure"].gt(0)]
    formula = f"raw_count ~ {predictors} + GC_content + log_gene_span"
    poisson = smf.glm(formula, data=work, family=sm.families.Poisson(), offset=work["log_exposure"]).fit()
    alpha = estimate_nb_alpha(work["raw_count"].to_numpy(float), poisson.fittedvalues.to_numpy(float))
    nb = smf.glm(
        formula,
        data=work,
        family=sm.families.NegativeBinomial(alpha=alpha),
        offset=work["log_exposure"],
    ).fit(cov_type="HC3")
    nb_frame = coefficient_frame(
        f"{model_prefix}_nb_{strain}_{outcome}",
        "negative_binomial_HC3",
        nb.params,
        nb.bse,
        nb.pvalues,
        nb.conf_int(),
        alpha=alpha,
        n_genes=len(work),
        outcome=outcome,
        strain=strain,
    )
    work["log_rate"] = np.log((work["raw_count"] + 0.5) / work["exposure"])
    ols = smf.ols(f"log_rate ~ {predictors} + GC_content + log_gene_span", data=work).fit(cov_type="HC3")
    ols_frame = coefficient_frame(
        f"{model_prefix}_ols_{strain}_{outcome}",
        "ols_log_rate_HC3",
        ols.params,
        ols.bse,
        ols.pvalues,
        ols.conf_int(),
        n_genes=len(work),
        outcome=outcome,
        strain=strain,
    )
    return pd.concat([nb_frame, ols_frame], ignore_index=True)


def permute_between_gene(data: pd.DataFrame, strain: str, n_perm: int, seed: int) -> pd.DataFrame:
    work = data[data["strain"].eq(strain)].copy()
    work["log_rate"] = np.log((work["count_all"] + 0.5) / work["exposure"])
    work = work[np.isfinite(work["log_rate"])].reset_index(drop=True)
    span_bin = pd.qcut(work["log_gene_span"], 5, duplicates="drop")
    gc_bin = pd.qcut(work["GC_content"], 5, duplicates="drop")
    work["bin"] = span_bin.astype(str) + "|" + gc_bin.astype(str)
    observed = smf.ols("log_rate ~ n_present + GC_content + log_gene_span", data=work).fit(cov_type="HC3")
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        shuffled = work.copy()
        for _, group in work.groupby("bin"):
            idx = group.index.to_numpy()
            take = rng.permutation(idx)
            shuffled.loc[idx, "n_present"] = work.loc[take, "n_present"].to_numpy()
        fit = smf.ols("log_rate ~ n_present + GC_content + log_gene_span", data=shuffled).fit()
        null[i] = fit.params["n_present"]
    observed_coef = float(observed.params["n_present"])
    return pd.DataFrame(
        [
            {
                "strain": strain,
                "term": "n_present",
                "observed_coefficient": observed_coef,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)),
                "two_sided_p": float(np.mean(np.abs(null) >= abs(observed_coef))),
                "n_permutations": n_perm,
            }
        ]
    )


def fit_pair_ppml(table: pd.DataFrame, count_columns: list[str], model_name: str) -> pd.DataFrame:
    terms = ["strain_RCC1614", *count_columns, "delta_GC"]
    x = np.column_stack(
        [np.ones(len(table)), *[table[column].to_numpy(float) for column in count_columns], table["delta_GC"].to_numpy(float)]
    )
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError(f"Nonidentifiable paired design for {model_name}")
    n = table["total_count"].to_numpy(float)
    y = table["count_rcc"].to_numpy(float)
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
    return coefficient_frame(model_name, "gene_fe_ppml_sandwich", params, bse, pvalues, conf, n_genes=len(table))


def fit_delta_log(table: pd.DataFrame, count_columns: list[str], model_name: str, outcome: str) -> pd.DataFrame:
    predictors = table[[*count_columns, "delta_GC"]]
    result = sm.OLS(table[f"delta_log_rate_{outcome}"], sm.add_constant(predictors)).fit(cov_type="HC3")
    return coefficient_frame(
        model_name,
        "ols_delta_log_HC3",
        result.params,
        result.bse,
        result.pvalues,
        result.conf_int(),
        n_genes=len(table),
        outcome=outcome,
    )


def sign_test(table: pd.DataFrame, outcome: str) -> dict:
    rate = table[f"delta_log_rate_{outcome}"]
    stable = table[table["delta_n"].eq(0)]
    offset = float(rate.loc[stable.index].mean()) if len(stable) else 0.0
    single = table[table["n_discordant"].eq(1) & table["delta_n"].abs().eq(1)]
    adjusted = rate.loc[single.index] - offset
    concordant = (adjusted * single["delta_n"]) < 0
    n = int(len(single))
    k = int(concordant.sum())
    test = binomtest(k, n, 0.5, alternative="greater") if n else None
    return {
        "outcome": outcome,
        "n_single_discordant_genes": n,
        "n_present_allele_lower": k,
        "strain_offset": offset,
        "binomial_p": float(test.pvalue) if test is not None else np.nan,
        "proportion": (k / n) if n else np.nan,
    }


def pair_for_outcome(pairs: pd.DataFrame, outcome: str) -> pd.DataFrame:
    table = pairs.copy()
    table["count_rcc"] = table[f"count_{outcome}_RCC1614"]
    table["total_count"] = table[f"count_{outcome}_CCMP1545"] + table[f"count_{outcome}_RCC1614"]
    return table[table["total_count"].gt(0)].copy()


def plot_overview(test_a, test_b, pairs, permutation, signs, output_path: Path) -> None:
    colors = load_color_guide()
    pop1 = get_color("Population 1", colors)
    pop2 = get_color("Population 2", colors)
    ancestor = get_color("Ancestor", colors)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), facecolor="white")
    for axis in axes.flat:
        axis.set_facecolor("white")

    ax = axes[0, 0]
    rows, labels, bar_colors = [], [], []
    palette = {"CCMP1545": pop1, "RCC1614": pop1, "RCC1749": pop2}
    for strain in STRAINS:
        row = test_a[
            test_a["model"].eq(f"test_a_nb_{strain}_all") & test_a["term"].eq("n_present")
        ].iloc[0]
        rows.append(row)
        labels.append(strain)
        bar_colors.append(palette[strain])
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
    ax.axvline(1, color=ancestor, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xscale("log")
    ax.set_xlabel("Expression rate ratio per present intron (95% CI)")
    ax.set_title("Test A: between-gene total burden (R2C2 sum)")

    ax = axes[0, 1]
    specs = [
        ("test_b_gene_fe_ppml_delta_n_all", "delta_n", "Gene FE PPML"),
        ("test_b_delta_log_delta_n_all", "delta_n", "Delta log OLS"),
        ("test_b_gene_fe_ppml_extra_missing_all", "n_extra", "FE extra vs CCMP1545"),
        ("test_b_gene_fe_ppml_extra_missing_all", "n_missing", "FE missing vs CCMP1545"),
        ("test_b_gene_fe_ppml_delta_n_family2_all", "delta_n_family2", "FE family 2 Δn"),
    ]
    b_rows, b_labels = [], []
    for model, term, label in specs:
        subset = test_b[test_b["model"].eq(model) & test_b["term"].eq(term)]
        if subset.empty:
            continue
        b_rows.append(subset.iloc[0])
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
    ax.axvline(1, color=ancestor, linestyle="--")
    ax.set_yticks(y, b_labels)
    ax.set_xscale("log")
    ax.set_xlabel("Expression ratio (95% CI)")
    ax.set_title("Test B: CCMP1545-relative occupancy")

    ax = axes[1, 0]
    jitter = np.random.default_rng(0).uniform(-0.08, 0.08, len(pairs))
    ax.scatter(
        pairs["delta_n"] + jitter,
        pairs["delta_log_rate_all"],
        s=np.where(pairs["delta_n"].eq(0), 8, 22),
        c=np.where(pairs["delta_n"].eq(0), ancestor, pop1),
        alpha=0.55,
        linewidths=0,
    )
    ax.axhline(0, color=ancestor, linestyle="--")
    ax.axvline(0, color=ancestor, linestyle="--")
    ax.set_xlabel("Δn (RCC1614 − CCMP1545 present count)")
    ax.set_ylabel("Δ log R2C2 rate")
    ax.set_title("Paired long-read expression vs occupancy")

    ax = axes[1, 1]
    ax.axis("off")
    lines = ["Test A permutation (OLS of log rate, length/GC bins)", ""]
    for row in permutation.itertuples(index=False):
        lines.append(f"{row.strain} n_present: observed {row.observed_coefficient:.3f}, null p={row.two_sided_p:.3f}")
    lines.append("")
    for sign in signs:
        lines.append(
            f"Sign test ({sign['outcome']}): {sign['n_present_allele_lower']}/"
            f"{sign['n_single_discordant_genes']} present allele lower (p={sign['binomial_p']:.3g})"
        )
    ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=10)
    ax.set_title("Null checks")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=base / "data")
    parser.add_argument("--output-dir", type=Path, default=base / "models")
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(args.data_dir / "expression_gene_strain.tsv", sep="\t")
    pairs = pd.read_csv(args.data_dir / "gene_pair_table.tsv", sep="\t")
    eligible = pairs[pairs["eligible_pair"]].copy()

    test_a_frames = []
    for strain in STRAINS:
        for outcome in OUTCOMES:
            test_a_frames.append(fit_between_gene(expression, strain, outcome, "n_present", "test_a"))
            test_a_frames.append(
                fit_between_gene(
                    expression,
                    strain,
                    outcome,
                    "n_present_family1 + n_present_family2 + n_present_family3 + n_present_other",
                    "test_a_family",
                )
            )
    test_a = pd.concat(test_a_frames, ignore_index=True)
    permutation = pd.concat(
        [permute_between_gene(expression, strain, args.n_permutations, RNG_SEED + i) for i, strain in enumerate(STRAINS)],
        ignore_index=True,
    )

    test_b_frames = []
    signs = []
    for outcome in OUTCOMES:
        table = pair_for_outcome(eligible, outcome)
        test_b_frames.append(fit_pair_ppml(table, ["delta_n"], f"test_b_gene_fe_ppml_delta_n_{outcome}"))
        test_b_frames.append(fit_delta_log(table, ["delta_n"], f"test_b_delta_log_delta_n_{outcome}", outcome))
        test_b_frames.append(
            fit_pair_ppml(table, ["n_extra", "n_missing"], f"test_b_gene_fe_ppml_extra_missing_{outcome}")
        )
        test_b_frames.append(
            fit_delta_log(table, ["n_extra", "n_missing"], f"test_b_delta_log_extra_missing_{outcome}", outcome)
        )
        test_b_frames.append(
            fit_pair_ppml(table, ["delta_n_family2"], f"test_b_gene_fe_ppml_delta_n_family2_{outcome}")
        )
        discordant = table[table["n_discordant"].gt(0)]
        test_b_frames.append(
            fit_pair_ppml(discordant, ["delta_n"], f"test_b_gene_fe_ppml_discordant_only_{outcome}")
        )
        test_b_frames.append(
            fit_delta_log(discordant, ["delta_n"], f"test_b_delta_log_discordant_only_{outcome}", outcome)
        )
        signs.append(sign_test(table, outcome))
        family_cols = ["delta_n_family1", "delta_n_family2", "delta_n_family3", "delta_n_familyother"]
        if np.linalg.matrix_rank(np.column_stack([np.ones(len(table)), table[family_cols], table["delta_GC"]])) == 6:
            test_b_frames.append(fit_pair_ppml(table, family_cols, f"test_b_gene_fe_ppml_family_joint_{outcome}"))
    test_b = pd.concat(test_b_frames, ignore_index=True)

    test_a.to_csv(output_dir / "test_a_coefficients.tsv", sep="\t", index=False)
    test_b.to_csv(output_dir / "test_b_coefficients.tsv", sep="\t", index=False)
    permutation.to_csv(output_dir / "test_a_permutation.tsv", sep="\t", index=False)
    pd.DataFrame(signs).to_csv(output_dir / "test_b_sign_test.tsv", sep="\t", index=False)
    support = pd.DataFrame(
        [
            {
                "model": "test_a",
                "genes_CCMP1545": int(expression["strain"].eq("CCMP1545").sum()),
                "genes_RCC1614": int(expression["strain"].eq("RCC1614").sum()),
                "genes_RCC1749": int(expression["strain"].eq("RCC1749").sum()),
                "genes_with_introners_CCMP1545": int(
                    (expression["strain"].eq("CCMP1545") & expression["n_present"].gt(0)).sum()
                ),
            },
            {
                "model": "test_b",
                "eligible_pairs": int(len(eligible)),
                "nonzero_delta_n": int((eligible["delta_n"] != 0).sum()),
                "n_extra_positive": int(eligible["n_extra"].gt(0).sum()),
                "n_missing_positive": int(eligible["n_missing"].gt(0).sum()),
                "family2_nonzero": int((eligible["delta_n_family2"] != 0).sum()),
            },
        ]
    )
    support.to_csv(output_dir / "model_support.tsv", sep="\t", index=False)
    plot_overview(test_a, test_b, eligible, permutation, signs, output_dir / "r2c2_burden_overview.png")

    def line(frame, model, term) -> str:
        row = frame[frame["model"].eq(model) & frame["term"].eq(term)].iloc[0]
        return (
            f"{model} | {term}: RR {row.rate_ratio:.3f} "
            f"({row.rr_ci_lower:.3f}-{row.rr_ci_upper:.3f}), p={row.p_value:.4g}"
        )

    checks = [
        "PASS: Test A uses total present count, not fixed/polymorphic class",
        "PASS: Test B reference is CCMP1545; RCC1749 is between-gene only",
        "PASS: R2C2 exposure is median-ratio size factor with no exon-length offset",
        f"PASS: eligible pairs have four replicates ({eligible.n_replicates_CCMP1545.eq(4).all() and eligible.n_replicates_RCC1614.eq(4).all()})",
        f"PASS: Δn equals n_extra - n_missing ({eligible.delta_n.eq(eligible.n_extra - eligible.n_missing).all()})",
    ]
    (output_dir / "validation_summary.txt").write_text("\n".join(checks) + "\n")

    summary = [
        "R2C2 total-burden expression pilot",
        "==================================",
        "",
        "Test A, between-gene NB of summed R2C2 counts:",
        "  " + line(test_a, "test_a_nb_CCMP1545_all", "n_present"),
        "  " + line(test_a, "test_a_nb_RCC1614_all", "n_present"),
        "  " + line(test_a, "test_a_nb_RCC1749_all", "n_present"),
        "",
        "Test A family split (NB, CCMP1545, all isoforms):",
        "  " + line(test_a, "test_a_family_nb_CCMP1545_all", "n_present_family1"),
        "  " + line(test_a, "test_a_family_nb_CCMP1545_all", "n_present_family2"),
        "  " + line(test_a, "test_a_family_nb_CCMP1545_all", "n_present_family3"),
        "  " + line(test_a, "test_a_family_nb_CCMP1545_all", "n_present_other"),
        "",
        "Test A OLS permutation of total burden:",
        *[
            f"  {row.strain}: observed {row.observed_coefficient:.3f}, null p={row.two_sided_p:.3f}"
            for row in permutation.itertuples(index=False)
        ],
        "",
        "Test B, CCMP1545 vs RCC1614, summed R2C2 counts:",
        "  " + line(test_b, "test_b_gene_fe_ppml_delta_n_all", "delta_n"),
        "  " + line(test_b, "test_b_delta_log_delta_n_all", "delta_n"),
        "  " + line(test_b, "test_b_gene_fe_ppml_extra_missing_all", "n_extra"),
        "  " + line(test_b, "test_b_gene_fe_ppml_extra_missing_all", "n_missing"),
        "  " + line(test_b, "test_b_gene_fe_ppml_delta_n_family2_all", "delta_n_family2"),
        "",
        "Major-isoform and FSM sensitivities (Δn FE):",
        "  " + line(test_b, "test_b_gene_fe_ppml_delta_n_major", "delta_n"),
        "  " + line(test_b, "test_b_gene_fe_ppml_delta_n_fsm", "delta_n"),
        "",
        f"Eligible pairs: {len(eligible)}; nonzero Δn: {int((eligible['delta_n'] != 0).sum())}; "
        f"extra: {int(eligible['n_extra'].gt(0).sum())}; missing: {int(eligible['n_missing'].gt(0).sum())}",
        "",
        f"Sign test (all isoforms): {signs[0]['n_present_allele_lower']}/{signs[0]['n_single_discordant_genes']} "
        f"present allele lower; p={signs[0]['binomial_p']:.4g}",
        "",
        "Interpretation guardrails:",
        "  - Test A is a between-gene occupancy association on collapsed long-read counts.",
        "  - Test B is a signed CCMP1545-relative occupancy difference within orthologs.",
        "  - RCC1749 is a between-gene check only, not an allelic contrast.",
        "  - Isoform number is not a covariate. Primary outcome sums all retained isoforms.",
    ]
    (output_dir / "model_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    print("\n".join(checks))


if __name__ == "__main__":
    main()

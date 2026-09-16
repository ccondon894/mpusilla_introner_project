#!/usr/bin/env python3
"""Create exploratory four-panel functional-association figure variants.

The three variants use the same source analyses but emphasize different layers:

1. data_forward: modeled isoform response plus observed distributions/fractions.
2. model_forward: effect-size and interval summaries throughout.
3. hybrid: modeled effects paired with compact observed summaries.

All outputs are 6.5-inch-wide PNG files intended for full-page placement in a
standard Google Doc. These are exploratory figures and are not workflow targets.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import expit
from scipy.stats import norm


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[3]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import get_color, load_color_guide  # noqa: E402


FIG_WIDTH = 6.5
FIG_HEIGHT = 6.4
AXIS_FONT = 8.5
TICK_FONT = 7.5
TITLE_FONT = 9.5
ANNOT_FONT = 7.2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_PATH.parent,
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=20260722)
    return parser.parse_args()


def p_text(value: float) -> str:
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        mantissa = value / (10 ** exponent)
        return rf"$p={mantissa:.1f}\times10^{{{exponent}}}$"
    return rf"$p={value:.2f}$"


def significance_label(value: float) -> str:
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def configure_axis(ax: plt.Axes, title: str, panel: str) -> None:
    ax.set_title(title, loc="left", fontsize=TITLE_FONT, fontweight="bold", pad=6)
    ax.text(
        -0.14,
        1.04,
        panel,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )
    ax.tick_params(
        axis="both",
        which="major",
        bottom=True,
        left=True,
        top=False,
        right=False,
        direction="out",
        labelsize=TICK_FONT,
        width=0.8,
        length=3,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.grid(False)


def load_isoform_effects(root: Path):
    directory = root / 'results/expression/glm_modeling/reference_independent_isoform_change'
    primary = pd.read_csv(directory/'conditional_poisson_coefficients.tsv',sep='\t')
    r = primary[primary.term.eq('current_introner_count')].iloc[0]
    rows = [dict(label='Within-gene change', model='conditional_poisson', estimate=r.coefficient,
                 se=r.std_error,p_value=r.p_value,ratio=r.rate_ratio,ci_low=r.rate_ratio_ci_lower,ci_high=r.rate_ratio_ci_upper)]
    robustness = pd.read_csv(directory/'within_between_negative_binomial.tsv',sep='\t')
    r = robustness[robustness.term.eq('current_introner_count_between')].iloc[0]
    rows.append(dict(label='Cross-gene burden',model='within_between_negative_binomial',estimate=r.coefficient,
                     se=r.std_error,p_value=r.p_value,ratio=r.rate_ratio,ci_low=r.rate_ratio_ci_lower,ci_high=r.rate_ratio_ci_upper))
    effects = pd.DataFrame(rows)
    return effects, effects.iloc[0].to_dict()


def fit_expression_contrasts(root: Path, replicates: int, rng):
    table = pd.read_csv(root/'results/expression/functional/short_read_models/model_coefficients.tsv',sep='\t')
    rows, draws = [], []
    for label, model, term in [
        ('Within-gene gain','paired_turnover_gene_fe_ppml','reference_relative_gain_count'),
        ('Within-gene loss','paired_turnover_gene_fe_ppml','reference_relative_loss_count'),
        ('Between-gene burden','current_state_nb_gee','current_introner_count_between')]:
        r = table[table.model.eq(model) & table.term.eq(term)].iloc[0]
        rows.append(dict(label=label,estimate=r.coefficient,se=r.std_error,p_value=r.p_value,
            ratio=r.rate_ratio,ci_low=r.rr_ci_lower,ci_high=r.rr_ci_upper,
            percent=100*(r.rate_ratio-1),percent_low=100*(r.rr_ci_lower-1),percent_high=100*(r.rr_ci_upper-1)))
        draws.extend(dict(label=label,percent=100*(np.exp(x)-1)) for x in rng.normal(r.coefficient,r.std_error,size=replicates))
    return pd.DataFrame(rows), pd.DataFrame(draws)


def bootstrap_raw_nmd(
    data: pd.DataFrame, replicates: int, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = (
        data.groupby(["common_gene_id", "has_introner"], sort=False)[
            ["n_nmd_positive", "n_nmd_evaluable"]
        ]
        .sum()
        .reset_index()
    )
    genes = grouped["common_gene_id"].unique()
    gene_index = {gene: index for index, gene in enumerate(genes)}
    positives = np.zeros((len(genes), 2), dtype=float)
    totals = np.zeros((len(genes), 2), dtype=float)
    for row in grouped.itertuples(index=False):
        index = gene_index[row.common_gene_id]
        category = int(row.has_introner)
        positives[index, category] = row.n_nmd_positive
        totals[index, category] = row.n_nmd_evaluable

    boot = np.full((replicates, 2), np.nan)
    for replicate in range(replicates):
        selected = rng.integers(0, len(genes), size=len(genes))
        p_sum = positives[selected].sum(axis=0)
        n_sum = totals[selected].sum(axis=0)
        boot[replicate] = p_sum / n_sum

    pooled = data.groupby("has_introner")[["n_nmd_positive", "n_nmd_evaluable"]].sum()
    rows = []
    draw_rows = []
    for category, label in [(0, "No introner"), (1, "≥1 introner")]:
        fraction = pooled.loc[category, "n_nmd_positive"] / pooled.loc[
            category, "n_nmd_evaluable"
        ]
        rows.append(
            {
                "label": label,
                "has_introner": category,
                "fraction": fraction,
                "ci_low": np.nanquantile(boot[:, category], 0.025),
                "ci_high": np.nanquantile(boot[:, category], 0.975),
                "n_positive": int(pooled.loc[category, "n_nmd_positive"]),
                "n_total": int(pooled.loc[category, "n_nmd_evaluable"]),
            }
        )
        draw_rows.extend(
            {
                "label": label,
                "fraction_percent": 100 * value,
            }
            for value in boot[:, category]
            if np.isfinite(value)
        )
    return pd.DataFrame(rows), pd.DataFrame(draw_rows)


def fit_adjusted_nmd(
    root: Path, replicates: int, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(
        root / "results/expression/functional/adjusted_nmd/nmd_gene_strain_data.tsv",
        sep="\t",
    )
    suffix = (
        "C(strain) + log1p_long_read_count + log_gene_span + "
        "background_exon_count + log1p_n_isoforms"
    )
    formula = "nmd_isoform_fraction ~ has_introner + " + suffix
    model = sm.GEE.from_formula(
        formula,
        groups="common_gene_id",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
        weights=data["n_nmd_evaluable"],
    )
    result = model.fit(cov_type="robust", maxiter=200)
    design_info = result.model.data.design_info
    weights = data["n_nmd_evaluable"].to_numpy(dtype=float)
    covariance = result.cov_params().to_numpy()
    covariance = (covariance + covariance.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    covariance_psd = eigenvectors @ np.diag(np.clip(eigenvalues, 0, None)) @ eigenvectors.T
    simulated_params = rng.multivariate_normal(
        result.params.to_numpy(), covariance_psd, size=replicates
    )

    adjusted_rows = []
    for category, label in [(0, "No introner"), (1, "≥1 introner")]:
        counterfactual = data.copy()
        counterfactual["has_introner"] = category
        design = np.asarray(
            patsy.build_design_matrices([design_info], counterfactual)[0]
        )
        estimate = np.average(expit(design @ result.params.to_numpy()), weights=weights)
        simulation_values = []
        for start in range(0, replicates, 250):
            chunk = simulated_params[start : start + 250]
            predictions = expit(design @ chunk.T)
            simulation_values.extend(np.average(predictions, axis=0, weights=weights))
        adjusted_rows.append(
            {
                "label": label,
                "has_introner": category,
                "fraction": estimate,
                "ci_low": np.quantile(simulation_values, 0.025),
                "ci_high": np.quantile(simulation_values, 0.975),
            }
        )

    coefficient_table = pd.read_csv(
        root / "results/expression/functional/adjusted_nmd/model_coefficients.tsv",
        sep="\t",
    )
    focal_rows = []
    selections = [
        (
            "Any current introner",
            "binary_architecture_adjusted",
            "has_introner",
        ),
        (
            "Per current introner",
            "burden_architecture_adjusted",
            "direct_current_introner_count",
        ),
    ]
    for label, model_name, term in selections:
        row = coefficient_table[
            coefficient_table["model"].eq(model_name)
            & coefficient_table["term"].eq(term)
        ].iloc[0]
        focal_rows.append(
            {
                "label": label,
                "ratio": row["odds_ratio"],
                "ci_low": row["or_ci_lower"],
                "ci_high": row["or_ci_upper"],
                "p_value": row["p_value"],
            }
        )
    return data, pd.DataFrame(adjusted_rows), pd.DataFrame(focal_rows)


def bootstrap_splicing_medians(
    data: pd.DataFrame, replicates: int, rng: np.random.Generator
) -> pd.DataFrame:
    labels = {
        "conventional_intron": "Conventional",
        "fixed_present": "Fixed introner",
        "polymorphic": "Polymorphic introner",
    }
    rows = []
    for category, label in labels.items():
        values = data.loc[data["analysis_class"].eq(category), "retention_psi"].to_numpy()
        bootstrap = np.empty(replicates)
        for replicate in range(replicates):
            bootstrap[replicate] = np.median(
                values[rng.integers(0, len(values), size=len(values))]
            )
        rows.append(
            {
                "analysis_class": category,
                "label": label,
                "n": len(values),
                "median": np.median(values),
                "q25": np.quantile(values, 0.25),
                "q75": np.quantile(values, 0.75),
                "ci_low": np.quantile(bootstrap, 0.025),
                "ci_high": np.quantile(bootstrap, 0.975),
            }
        )
    return pd.DataFrame(rows)


def load_splicing(
    root: Path, replicates: int, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = root / "results/expression/functional/corrected_retention_psi"
    loci = pd.read_csv(base / "per_locus_scores.tsv", sep="\t", low_memory=False)
    eligible = loci["primary_boundary_eligible"].astype(bool)
    confirmed = loci[eligible & (loci["junction_signal_exact"] >= 5)].copy()
    keep = ["conventional_intron", "fixed_present", "polymorphic"]
    confirmed = confirmed[confirmed["analysis_class"].isin(keep)].copy()
    summaries = bootstrap_splicing_medians(confirmed, replicates, rng)
    tests = pd.read_csv(base / "evidence_guarded_class_tests.tsv", sep="\t")
    return confirmed, summaries, tests


def plot_isoform_curve(ax: plt.Axes, within: dict[str, float], color: str,
                       hybrid: bool = False) -> None:
    changes = np.arange(-2, 3)
    beta = within["estimate"]
    lower_beta = beta - 1.96 * within["se"]
    upper_beta = beta + 1.96 * within["se"]
    estimate = np.exp(changes * beta)
    bounds = np.vstack([np.exp(changes * lower_beta), np.exp(changes * upper_beta)])
    lower = bounds.min(axis=0)
    upper = bounds.max(axis=0)
    ax.fill_between(changes, lower, upper, color=color, alpha=0.18, linewidth=0)
    ax.plot(changes, estimate, color=color, marker="o", linewidth=1.5, markersize=4)
    ax.axhline(1, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_xticks(changes)
    ax.set_xlabel("Within-gene change in introner number", fontsize=AXIS_FONT)
    ax.set_ylabel("Relative alternative-isoform count", fontsize=AXIS_FONT)
    annotation = (
        f"+1 introner: RR={within['ratio']:.2f}\n"
        f"95% CI {within['ci_low']:.2f}–{within['ci_high']:.2f}\n"
        f"{p_text(within['p_value'])}"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=ANNOT_FONT,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2},
    )
    if hybrid:
        ax.set_ylim(0.55, 1.75)


def plot_isoform_forest(ax: plt.Axes, effects: pd.DataFrame, color: str) -> None:
    positions = np.arange(len(effects))[::-1]
    ax.errorbar(
        effects["ratio"],
        positions,
        xerr=[effects["ratio"] - effects["ci_low"], effects["ci_high"] - effects["ratio"]],
        fmt="o",
        color=color,
        ecolor=color,
        capsize=2.5,
        markersize=4.5,
    )
    ax.axvline(1, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_yticks(positions, effects["label"])
    ax.set_xlabel("Rate ratio per current introner", fontsize=AXIS_FONT)
    ax.set_xlim(0.82, 1.43)


def plot_expression_percent(ax: plt.Axes, effects: pd.DataFrame, colors: list[str]) -> None:
    positions = np.arange(len(effects))[::-1]
    values = effects["percent"].to_numpy()
    ax.barh(positions, values, color=colors, edgecolor="black", linewidth=0.7, height=0.58)
    ax.errorbar(
        values,
        positions,
        xerr=[values - effects["percent_low"], effects["percent_high"] - values],
        fmt="none",
        ecolor="black",
        capsize=2.5,
        linewidth=0.9,
    )
    ax.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    labels = [
        "Current burden\n(+1)",
        "Ref-relative\ngain",
        "Ref-relative\nloss",
    ]
    ax.set_yticks(positions, labels)
    plt.setp(
        ax.get_yticklabels(),
        rotation=90,
        ha="center",
        va="center",
        rotation_mode="anchor",
        linespacing=1.65,
    )
    ax.tick_params(axis="y", pad=12)
    ax.set_xlabel("Expected expression change (%)", fontsize=AXIS_FONT)
    ax.set_xlim(-60, 38)


def plot_expression_violin(
    ax: plt.Axes,
    draws: pd.DataFrame,
    effects: pd.DataFrame,
    colors: list[str],
) -> None:
    order = [
        "Between-gene burden",
        "Within-gene gain",
        "Within-gene loss",
    ]
    palette = dict(zip(order, colors))
    sns.violinplot(
        data=draws,
        x="label",
        y="percent",
        hue="label",
        order=order,
        hue_order=order,
        palette=palette,
        inner="quartile",
        cut=0,
        density_norm="width",
        linewidth=0.8,
        legend=False,
        ax=ax,
    )
    estimates = effects.set_index("label").loc[order, "percent"].to_numpy()
    ax.scatter(
        np.arange(len(order)),
        estimates,
        s=18,
        facecolor="white",
        edgecolor="black",
        linewidth=0.7,
        zorder=4,
    )
    ax.axhline(0, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_xticks(
        np.arange(len(order)),
        ["Current burden\n(+1)", "Ref-relative\ngain", "Ref-relative\nloss"],
    )
    ax.set_xlabel("")
    ax.set_ylabel("Expected expression change (%)", fontsize=AXIS_FONT)
    lower = np.nanquantile(draws["percent"], 0.002)
    upper = np.nanquantile(draws["percent"], 0.998)
    padding = 0.08 * (upper - lower)
    ax.set_ylim(lower - padding, upper + padding)


def plot_expression_forest(ax: plt.Axes, effects: pd.DataFrame, colors: list[str]) -> None:
    positions = np.arange(len(effects))[::-1]
    for y, (_, row), color in zip(positions, effects.iterrows(), colors):
        ax.errorbar(
            row["ratio"],
            y,
            xerr=[[row["ratio"] - row["ci_low"]], [row["ci_high"] - row["ratio"]]],
            fmt="o",
            color=color,
            ecolor=color,
            capsize=2.5,
            markersize=4.5,
        )
    ax.axvline(1, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_yticks(positions, effects["label"])
    ax.set_xlabel("Expected expression ratio", fontsize=AXIS_FONT)
    ax.set_xlim(0.45, 1.38)


def plot_nmd_fractions(
    ax: plt.Axes,
    values: pd.DataFrame,
    colors: list[str],
    odds: pd.Series,
    adjusted: pd.DataFrame | None = None,
) -> None:
    positions = np.arange(len(values))
    y = 100 * values["fraction"].to_numpy()
    low = 100 * values["ci_low"].to_numpy()
    high = 100 * values["ci_high"].to_numpy()
    if adjusted is None:
        ax.bar(positions, y, color=colors, edgecolor="black", linewidth=0.7, width=0.62)
        ax.errorbar(
            positions,
            y,
            yerr=[y - low, high - y],
            fmt="none",
            ecolor="black",
            capsize=3,
            linewidth=0.9,
        )
    else:
        adjusted_y = 100 * adjusted["fraction"].to_numpy()
        adjusted_low = 100 * adjusted["ci_low"].to_numpy()
        adjusted_high = 100 * adjusted["ci_high"].to_numpy()
        for index in positions:
            ax.errorbar(
                index - 0.09,
                y[index],
                yerr=[[y[index] - low[index]], [high[index] - y[index]]],
                fmt="o",
                markerfacecolor="white",
                markeredgecolor=colors[index],
                ecolor=colors[index],
                capsize=2,
                label="Observed" if index == 0 else None,
            )
            ax.errorbar(
                index + 0.09,
                adjusted_y[index],
                yerr=[
                    [adjusted_y[index] - adjusted_low[index]],
                    [adjusted_high[index] - adjusted_y[index]],
                ],
                fmt="o",
                color=colors[index],
                ecolor=colors[index],
                capsize=2,
                label="Adjusted" if index == 0 else None,
            )
        ax.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax.set_xticks(positions, values["label"])
    ax.set_ylabel("NMD-positive isoforms (%)", fontsize=AXIS_FONT)
    ax.set_ylim(0, max(7.5, high.max() * 1.55))
    ax.text(
        0.97,
        0.96,
        f"Adjusted OR={odds['ratio']:.2f}\n"
        f"95% CI {odds['ci_low']:.2f}–{odds['ci_high']:.2f}\n"
        f"{p_text(odds['p_value'])}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=ANNOT_FONT,
    )


def plot_nmd_violin(
    ax: plt.Axes,
    draws: pd.DataFrame,
    values: pd.DataFrame,
    colors: list[str],
    odds: pd.Series,
) -> None:
    order = ["No introner", "≥1 introner"]
    palette = dict(zip(order, colors))
    sns.violinplot(
        data=draws,
        x="label",
        y="fraction_percent",
        hue="label",
        order=order,
        hue_order=order,
        palette=palette,
        inner="quartile",
        cut=0,
        density_norm="width",
        linewidth=0.8,
        legend=False,
        ax=ax,
    )
    pooled = values.set_index("label").loc[order, "fraction"].to_numpy() * 100
    ax.scatter(
        np.arange(len(order)),
        pooled,
        s=18,
        facecolor="white",
        edgecolor="black",
        linewidth=0.7,
        zorder=4,
    )
    ax.set_xticks(np.arange(len(order)), order)
    ax.set_xlabel("")
    ax.set_ylabel("NMD-positive isoforms (%)", fontsize=AXIS_FONT)
    upper = np.nanquantile(draws["fraction_percent"], 0.998)
    ax.set_ylim(0, max(7.5, upper * 1.35))
    ax.text(
        0.97,
        0.96,
        f"Adjusted OR={odds['ratio']:.2f}\n"
        f"95% CI {odds['ci_low']:.2f}–{odds['ci_high']:.2f}\n"
        f"{p_text(odds['p_value'])}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=ANNOT_FONT,
    )


def plot_nmd_forest(ax: plt.Axes, effects: pd.DataFrame, color: str) -> None:
    positions = np.arange(len(effects))[::-1]
    ax.errorbar(
        effects["ratio"],
        positions,
        xerr=[effects["ratio"] - effects["ci_low"], effects["ci_high"] - effects["ratio"]],
        fmt="o",
        color=color,
        ecolor=color,
        capsize=2.5,
        markersize=4.5,
    )
    ax.axvline(1, color="#555555", linestyle="--", linewidth=0.9)
    ax.set_yticks(positions, effects["label"])
    ax.set_xlabel("Odds ratio for NMD-positive status", fontsize=AXIS_FONT)
    ax.set_xlim(0.75, 4.05)


def add_bracket(ax: plt.Axes, x1: float, x2: float, y: float, height: float,
                label: str) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y],
            color="black", linewidth=0.8)
    ax.text((x1 + x2) / 2, y + height, label, ha="center", va="bottom",
            fontsize=ANNOT_FONT, fontweight="bold")


def plot_splicing_violin(
    ax: plt.Axes,
    data: pd.DataFrame,
    summaries: pd.DataFrame,
    tests: pd.DataFrame,
    palette: dict[str, str],
    conventional_label: str = "Conventional",
    show_medians: bool = True,
) -> None:
    labels = {
        "conventional_intron": conventional_label,
        "fixed_present": "Fixed introner",
        "polymorphic": "Polymorphic introner",
    }
    order = list(labels)
    plot_data = data.copy()
    plot_data["display"] = plot_data["analysis_class"].map(labels)
    plot_data["log_retention"] = np.log10(plot_data["retention_psi"] + 0.001)
    display_order = [labels[item] for item in order]
    display_palette = {labels[key]: value for key, value in palette.items()}
    sns.violinplot(
        data=plot_data,
        x="display",
        y="log_retention",
        hue="display",
        order=display_order,
        hue_order=display_order,
        palette=display_palette,
        inner="quartile",
        cut=0,
        density_norm="width",
        linewidth=0.8,
        legend=False,
        ax=ax,
    )
    raw_ticks = np.array([0, 0.01, 0.1, 1.0])
    ax.set_yticks(np.log10(raw_ticks + 0.001), ["0", "0.01", "0.1", "1.0"])
    ax.set_ylabel("Retention PSI (log scale)", fontsize=AXIS_FONT)
    ax.set_xlabel("")
    ax.set_xticks(
        np.arange(3),
        [conventional_label, "Fixed\nintroner", "Polymorphic\nintroner"],
    )
    ax.set_ylim(-3.12, 0.43)

    fixed_vs_conv = tests[
        tests["focal_class"].eq("fixed_present")
        & tests["reference_class"].eq("conventional_intron")
    ].iloc[0]
    poly_vs_fixed = tests[
        tests["focal_class"].eq("polymorphic")
        & tests["reference_class"].eq("fixed_present")
    ].iloc[0]
    add_bracket(
        ax,
        0,
        1,
        0.08,
        0.07,
        significance_label(fixed_vs_conv["mann_whitney_p"]),
    )
    add_bracket(
        ax,
        1,
        2,
        0.25,
        0.07,
        significance_label(poly_vs_fixed["mann_whitney_p"]),
    )
    if show_medians:
        ax.text(
            1,
            -2.96,
            f"Medians: {summaries.iloc[0]['median']:.3f} / "
            f"{summaries.iloc[1]['median']:.3f} / {summaries.iloc[2]['median']:.3f}",
            ha="center",
            va="bottom",
            fontsize=6.4,
        )


def plot_splicing_summary(
    ax: plt.Axes,
    summaries: pd.DataFrame,
    tests: pd.DataFrame,
    colors: list[str],
) -> None:
    positions = np.arange(len(summaries))
    for x, (_, row), color in zip(positions, summaries.iterrows(), colors):
        ax.errorbar(
            x,
            row["median"],
            yerr=[[row["median"] - row["ci_low"]], [row["ci_high"] - row["median"]]],
            fmt="o",
            color=color,
            ecolor=color,
            markeredgecolor="black",
            markeredgewidth=0.5,
            capsize=3,
            markersize=5,
        )
        ax.vlines(x, row["q25"], row["q75"], color=color, linewidth=4, alpha=0.35)
    ax.set_xticks(positions, ["Conventional", "Fixed\nintroner", "Polymorphic\nintroner"])
    ax.set_ylabel("Median retention PSI", fontsize=AXIS_FONT)
    ax.set_ylim(0, 0.082)
    poly_vs_fixed = tests[
        tests["focal_class"].eq("polymorphic")
        & tests["reference_class"].eq("fixed_present")
    ].iloc[0]
    ax.text(
        0.98,
        0.96,
        f"Polymorphic vs fixed\n{p_text(poly_vs_fixed['mann_whitney_p'])}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=ANNOT_FONT,
    )


def save_variant(fig: plt.Figure, output: Path) -> None:
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_data_forward(
    output: Path,
    colors: dict[str, str],
    iso_within: dict[str, float],
    expression: pd.DataFrame,
    expression_draws: pd.DataFrame,
    nmd_raw: pd.DataFrame,
    nmd_draws: pd.DataFrame,
    nmd_effects: pd.DataFrame,
    splicing: pd.DataFrame,
    splicing_summary: pd.DataFrame,
    splicing_tests: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH, FIG_HEIGHT),
                             constrained_layout=True)
    fig.set_constrained_layout_pads(
        w_pad=0.015,
        h_pad=0.025,
        wspace=0.025,
        hspace=0.055,
    )
    plot_isoform_curve(axes[0, 0], iso_within, colors["introner"])
    plot_expression_violin(
        axes[0, 1],
        expression_draws,
        expression,
        [colors["introner_dark"], colors["introner"], colors["neutral"]],
    )
    plot_nmd_violin(
        axes[1, 0],
        nmd_draws,
        nmd_raw,
        [colors["neutral"], colors["introner"]],
        nmd_effects.iloc[0],
    )
    plot_splicing_violin(
        axes[1, 1],
        splicing,
        splicing_summary,
        splicing_tests,
        {
            "conventional_intron": colors["intron"],
            "fixed_present": colors["fixed"],
            "polymorphic": colors["polymorphic"],
        },
        conventional_label="Intron",
        show_medians=False,
    )
    for ax, title, panel in zip(
        axes.flat,
        ["", "", "", ""],
        list("ABCD"),
    ):
        configure_axis(ax, title, panel)
    save_variant(fig, output)


def make_model_forward(
    output: Path,
    colors: dict[str, str],
    iso_effects: pd.DataFrame,
    expression: pd.DataFrame,
    nmd_effects: pd.DataFrame,
    splicing_summary: pd.DataFrame,
    splicing_tests: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH, FIG_HEIGHT),
                             constrained_layout=True)
    plot_isoform_forest(axes[0, 0], iso_effects, colors["introner"])
    plot_expression_forest(
        axes[0, 1],
        expression,
        [colors["introner_dark"], colors["introner"], colors["neutral"]],
    )
    plot_nmd_forest(axes[1, 0], nmd_effects, colors["introner"])
    plot_splicing_summary(
        axes[1, 1],
        splicing_summary,
        splicing_tests,
        [colors["intron"], colors["fixed"], colors["polymorphic"]],
    )
    for ax, title, panel in zip(
        axes.flat,
        ["Isoform-model estimates", "Expression-model estimates",
         "NMD-model estimates", "Retention PSI summaries"],
        list("ABCD"),
    ):
        configure_axis(ax, title, panel)
    save_variant(fig, output)


def make_hybrid(
    output: Path,
    colors: dict[str, str],
    iso_effects: pd.DataFrame,
    iso_within: dict[str, float],
    expression: pd.DataFrame,
    nmd_raw: pd.DataFrame,
    nmd_adjusted: pd.DataFrame,
    nmd_effects: pd.DataFrame,
    splicing: pd.DataFrame,
    splicing_summary: pd.DataFrame,
    splicing_tests: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH, FIG_HEIGHT),
                             constrained_layout=True)
    plot_isoform_curve(axes[0, 0], iso_within, colors["introner"], hybrid=True)
    cross = iso_effects.set_index("label")
    axes[0, 0].text(
        0.97,
        0.05,
        f"Cross-gene RR={cross.loc['Cross-gene burden', 'ratio']:.2f}\n"
        "",
        transform=axes[0, 0].transAxes,
        ha="right",
        va="bottom",
        fontsize=ANNOT_FONT,
        color="#333333",
    )
    plot_expression_forest(
        axes[0, 1],
        expression,
        [colors["introner_dark"], colors["introner"], colors["neutral"]],
    )
    for y, (_, row) in zip(np.arange(len(expression))[::-1], expression.iterrows()):
        axes[0, 1].text(
            row["ci_high"] + 0.025,
            y,
            f"{row['percent']:+.0f}%",
            va="center",
            fontsize=ANNOT_FONT,
        )
    plot_nmd_fractions(
        axes[1, 0],
        nmd_raw,
        [colors["neutral"], colors["introner"]],
        nmd_effects.iloc[0],
        adjusted=nmd_adjusted,
    )
    plot_splicing_violin(
        axes[1, 1],
        splicing,
        splicing_summary,
        splicing_tests,
        {
            "conventional_intron": colors["intron"],
            "fixed_present": colors["fixed"],
            "polymorphic": colors["polymorphic"],
        },
    )
    for ax, title, panel in zip(
        axes.flat,
        ["Within-gene isoform change", "Expression association",
         "Observed and adjusted NMD", "Retention PSI distributions"],
        list("ABCD"),
    ):
        configure_axis(ax, title, panel)
    save_variant(fig, output)


def write_statistics(
    output: Path,
    isoform: pd.DataFrame,
    expression: pd.DataFrame,
    nmd_raw: pd.DataFrame,
    nmd_adjusted: pd.DataFrame,
    nmd_effects: pd.DataFrame,
    splicing: pd.DataFrame,
    splicing_tests: pd.DataFrame,
) -> None:
    rows = []
    for table_name, table in [
        ("isoform_model", isoform),
        ("expression_contrast", expression),
        ("nmd_raw_fraction", nmd_raw),
        ("nmd_adjusted_fraction", nmd_adjusted),
        ("nmd_model", nmd_effects),
        ("splicing_summary", splicing),
    ]:
        for record in table.to_dict("records"):
            record = {"analysis": table_name, **record}
            rows.append(record)
    for record in splicing_tests.to_dict("records"):
        rows.append({"analysis": "splicing_test", **record})
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    guide = load_color_guide(root / "master_figure_color_guide.tsv")
    colors = {
        "introner": get_color("Introner", guide),
        "introner_dark": get_color("Population 1 dark", guide),
        "neutral": get_color("Ancestor", guide),
        "intron": get_color("Intron", guide),
        "polymorphic": get_color("Population 1 Polymorphic", guide),
        "fixed": get_color("Population 1 Fixed", guide),
    }
    sns.set_theme(style="white", context="paper")

    isoform_effects, isoform_within = load_isoform_effects(root)
    expression_effects, expression_draws = fit_expression_contrasts(
        root,
        args.bootstrap_replicates,
        rng,
    )
    nmd_data, nmd_adjusted, nmd_effects = fit_adjusted_nmd(
        root, args.bootstrap_replicates, rng
    )
    nmd_raw, nmd_draws = bootstrap_raw_nmd(
        nmd_data,
        args.bootstrap_replicates,
        rng,
    )
    splicing_data, splicing_summary, splicing_tests = load_splicing(
        root, args.bootstrap_replicates, rng
    )

    make_data_forward(
        output_dir / "functional_associations.data_forward.violin_bc.png",
        colors,
        isoform_within,
        expression_effects,
        expression_draws,
        nmd_raw,
        nmd_draws,
        nmd_effects,
        splicing_data,
        splicing_summary,
        splicing_tests,
    )
    make_model_forward(
        output_dir / "functional_associations.model_forward.png",
        colors,
        isoform_effects,
        expression_effects,
        nmd_effects,
        splicing_summary,
        splicing_tests,
    )
    make_hybrid(
        output_dir / "functional_associations.hybrid.png",
        colors,
        isoform_effects,
        isoform_within,
        expression_effects,
        nmd_raw,
        nmd_adjusted,
        nmd_effects,
        splicing_data,
        splicing_summary,
        splicing_tests,
    )
    write_statistics(
        output_dir / "functional_associations.plotted_statistics.tsv",
        isoform_effects,
        expression_effects,
        nmd_raw,
        nmd_adjusted,
        nmd_effects,
        splicing_summary,
        splicing_tests,
    )
    print(f"Wrote three figure variants and plotted statistics to {output_dir}")


if __name__ == "__main__":
    main()

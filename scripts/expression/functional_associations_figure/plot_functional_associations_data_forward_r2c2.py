#!/usr/bin/env python3
"""Copy of the data-forward functional-associations plot with R2C2 panel B.

This alternate preserves panels A, C, and D and the 2x2 data-forward layout
from ``plot_functional_associations_variants.py``. Only panel B changes: its
three original short-read expression contrasts are replaced by the requested
R2C2 between-gene burden and reference-relative gain estimates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import plot_functional_associations_variants as original


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            SCRIPT_PATH.parent
            / "functional_associations.data_forward.violin_bc.r2c2.png"
        ),
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=20260722)
    return parser.parse_args()


def select_one(table: pd.DataFrame, description: str, **filters: str) -> pd.Series:
    selected = table.copy()
    for column, value in filters.items():
        selected = selected[selected[column].eq(value)]
    if len(selected) != 1:
        criteria = ", ".join(f"{key}={value!r}" for key, value in filters.items())
        raise ValueError(
            f"Expected one {description} row ({criteria}); found {len(selected)}"
        )
    return selected.iloc[0]


def load_r2c2_panel_b(
    root: Path,
    replicates: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = root / "results/expression/functional/r2c2_between_gene"

    tests = pd.read_csv(result_dir / "model_coefficients.tsv", sep="\t")
    between = select_one(
        tests,
        "between-gene NB-GEE",
        model="current_state_nb_gee",
        term="current_introner_count_between",
    )

    comparison = pd.read_csv(root / "results/expression/glm_modeling/within_gene_gain/gain_effect.tsv", sep="\t")
    gain = select_one(
        comparison,
        "reference-relative gain gene-FE",
        model="paired_turnover_gene_fe_ppml",
        term="reference_relative_gain_count",
    )

    specifications = [
        (
            "Between-gene burden (NB-GEE)",
            between,
            "r2c2_between_gene/model_coefficients.tsv",
        ),
        (
            "Reference-relative gain (gene FE)",
            gain,
            "within_gene_gain/gain_effect.tsv",
        ),
    ]
    rows = []
    draw_rows = []
    for label, source, source_table in specifications:
        coefficient = float(source["coefficient"])
        standard_error = float(source["std_error"])
        ratio = float(source["rate_ratio"])
        ci_low = float(source["rr_ci_lower"])
        ci_high = float(source["rr_ci_upper"])
        rows.append(
            {
                "label": label,
                "estimate": coefficient,
                "se": standard_error,
                "p_value": float(source["p_value"]),
                "ratio": ratio,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "percent": 100 * (ratio - 1),
                "percent_low": 100 * (ci_low - 1),
                "percent_high": 100 * (ci_high - 1),
                "source_table": source_table,
                "model": str(source["model"]),
                "term": str(source["term"]),
            }
        )
        log_effects = rng.normal(coefficient, standard_error, size=replicates)
        draw_rows.extend(
            {
                "label": label,
                "percent": 100 * (np.exp(log_effect) - 1),
            }
            for log_effect in log_effects
        )
    return pd.DataFrame(rows), pd.DataFrame(draw_rows)


def plot_r2c2_expression_violin(
    ax,
    draws: pd.DataFrame,
    effects: pd.DataFrame,
    colors: list[str],
) -> None:
    order = effects["label"].tolist()
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
        ["Between-gene burden\n(NB-GEE)", "Reference-relative gain\n(gene FE)"],
    )
    ax.set_xlabel("")
    ax.set_ylabel("Expected expression change (%)", fontsize=original.AXIS_FONT)
    lower = np.nanquantile(draws["percent"], 0.002)
    upper = np.nanquantile(draws["percent"], 0.998)
    padding = 0.08 * (upper - lower)
    ax.set_ylim(lower - padding, upper + padding)


def make_data_forward_r2c2(
    output: Path,
    colors: dict[str, str],
    isoform_within: dict[str, float],
    expression: pd.DataFrame,
    expression_draws: pd.DataFrame,
    nmd_raw: pd.DataFrame,
    nmd_draws: pd.DataFrame,
    nmd_effects: pd.DataFrame,
    splicing: pd.DataFrame,
    splicing_summary: pd.DataFrame,
    splicing_tests: pd.DataFrame,
) -> None:
    # Copied from make_data_forward; panel B calls the R2C2-specific renderer.
    fig, axes = original.plt.subplots(
        2,
        2,
        figsize=(original.FIG_WIDTH, original.FIG_HEIGHT),
        constrained_layout=True,
    )
    fig.set_constrained_layout_pads(
        w_pad=0.015,
        h_pad=0.025,
        wspace=0.025,
        hspace=0.055,
    )
    original.plot_isoform_curve(axes[0, 0], isoform_within, colors["introner"])
    plot_r2c2_expression_violin(
        axes[0, 1],
        expression_draws,
        expression,
        [colors["introner_dark"], colors["introner"]],
    )
    original.plot_nmd_violin(
        axes[1, 0],
        nmd_draws,
        nmd_raw,
        [colors["neutral"], colors["introner"]],
        nmd_effects.iloc[0],
    )
    original.plot_splicing_violin(
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
    for ax, panel in zip(axes.flat, "ABCD"):
        original.configure_axis(ax, "", panel)
    original.save_variant(fig, output)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    panel_b_rng = np.random.default_rng(args.seed + 1)

    guide = original.load_color_guide(root / "master_figure_color_guide.tsv")
    colors = {
        "introner": original.get_color("Introner", guide),
        "introner_dark": original.get_color("Population 1 dark", guide),
        "neutral": original.get_color("Ancestor", guide),
        "intron": original.get_color("Intron", guide),
        "polymorphic": original.get_color("Population 1 Polymorphic", guide),
        "fixed": original.get_color("Population 1 Fixed", guide),
    }
    sns.set_theme(style="white", context="paper")

    _, isoform_within = original.load_isoform_effects(root)
    # Advance the shared RNG exactly as the source figure does before its NMD
    # and retention bootstraps, keeping panels C and D reproducible.
    rng.normal(size=3 * args.bootstrap_replicates)
    expression, expression_draws = load_r2c2_panel_b(
        root, args.bootstrap_replicates, panel_b_rng
    )
    nmd_data, _, nmd_effects = original.fit_adjusted_nmd(
        root, args.bootstrap_replicates, rng
    )
    nmd_raw, nmd_draws = original.bootstrap_raw_nmd(
        nmd_data, args.bootstrap_replicates, rng
    )
    splicing, splicing_summary, splicing_tests = original.load_splicing(
        root, args.bootstrap_replicates, rng
    )

    make_data_forward_r2c2(
        output,
        colors,
        isoform_within,
        expression,
        expression_draws,
        nmd_raw,
        nmd_draws,
        nmd_effects,
        splicing,
        splicing_summary,
        splicing_tests,
    )
    statistics_output = output.with_suffix(".panel_b.tsv")
    expression.to_csv(statistics_output, sep="\t", index=False)
    print(f"Wrote {output}")
    print(f"Wrote {statistics_output}")


if __name__ == "__main__":
    main()

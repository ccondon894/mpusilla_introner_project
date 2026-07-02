#!/usr/bin/env python3
"""Plot selected pyrho recombination categories as violin distributions."""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide


CATEGORY_ORDER = [
    "Exonic background\nPopulation 1",
    "All introners\nPopulation 1",
    "Polymorphic introners\nPopulation 1",
    "Exonic background\nPopulation 2",
    "All introners\nPopulation 2",
]

PLANNED_COMPARISONS = [
    ("Exonic background\nPopulation 1", "All introners\nPopulation 1"),
    ("All introners\nPopulation 1", "Polymorphic introners\nPopulation 1"),
    ("Exonic background\nPopulation 2", "All introners\nPopulation 2"),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group1-all-windows", required=True)
    parser.add_argument("--group1-polymorphic-windows", required=True)
    parser.add_argument("--group2-windows", required=True)
    parser.add_argument("--output-pdf", required=True)
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--color-guide", default=str(DEFAULT_GUIDE_PATH))
    parser.add_argument(
        "--min-rate",
        type=float,
        default=1e-14,
        help="Minimum positive recombination rate retained for plotting.",
    )
    parser.add_argument(
        "--max-rate",
        type=float,
        default=5e-11,
        help="Maximum recombination rate retained for plotting.",
    )
    return parser.parse_args()


def load_category(path, window_type, label, population, min_rate, max_rate):
    df = pd.read_csv(path, sep="\t")
    required = {"type", "rate"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")

    subset = df[df["type"].eq(window_type)].copy()
    subset["rate"] = pd.to_numeric(subset["rate"], errors="coerce")
    subset = subset.dropna(subset=["rate"])
    subset = subset[(subset["rate"] > 0) & (subset["rate"] >= min_rate)]
    subset = subset[subset["rate"] <= max_rate]
    subset["category"] = label
    subset["population"] = population
    subset["log10_rate"] = np.log10(subset["rate"])
    return subset


def build_plot_data(args):
    frames = [
        load_category(
            args.group1_all_windows,
            "non_introner_containing",
            "Exonic background\nPopulation 1",
            "Population 1",
            args.min_rate,
            args.max_rate,
        ),
        load_category(
            args.group1_all_windows,
            "introner_containing",
            "All introners\nPopulation 1",
            "Population 1",
            args.min_rate,
            args.max_rate,
        ),
        load_category(
            args.group1_polymorphic_windows,
            "introner_containing",
            "Polymorphic introners\nPopulation 1",
            "Population 1",
            args.min_rate,
            args.max_rate,
        ),
        load_category(
            args.group2_windows,
            "non_introner_containing",
            "Exonic background\nPopulation 2",
            "Population 2",
            args.min_rate,
            args.max_rate,
        ),
        load_category(
            args.group2_windows,
            "introner_containing",
            "All introners\nPopulation 2",
            "Population 2",
            args.min_rate,
            args.max_rate,
        ),
    ]
    return pd.concat(frames, ignore_index=True)


def category_colors(color_guide):
    return {
        "Exonic background\nPopulation 1": get_color("Population 1 light", color_guide),
        "All introners\nPopulation 1": get_color("Introner", color_guide),
        "Polymorphic introners\nPopulation 1": get_color(
            "Population 1 Polymorphic", color_guide
        ),
        "Exonic background\nPopulation 2": get_color("Population 2 light", color_guide),
        "All introners\nPopulation 2": get_color("Population 2", color_guide),
    }


def print_summary(plot_df):
    print("Recombination violin categories:")
    for category in CATEGORY_ORDER:
        values = plot_df.loc[plot_df["category"].eq(category), "rate"]
        if values.empty:
            print(f"  {category}: n=0")
            continue
        print(
            f"  {category}: n={len(values)}, "
            f"median={values.median():.3e}, mean={values.mean():.3e}"
        )


def significance_label(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def planned_tests(plot_df):
    results = []
    correction = len(PLANNED_COMPARISONS)
    for left, right in PLANNED_COMPARISONS:
        left_values = plot_df.loc[plot_df["category"].eq(left), "log10_rate"].dropna()
        right_values = plot_df.loc[plot_df["category"].eq(right), "log10_rate"].dropna()
        if left_values.empty or right_values.empty:
            results.append((left, right, np.nan, np.nan, len(left_values), len(right_values)))
            continue
        stat, p_value = mannwhitneyu(left_values, right_values, alternative="two-sided")
        p_corrected = min(p_value * correction, 1.0)
        results.append((left, right, stat, p_corrected, len(left_values), len(right_values)))
    return results


def add_significance_bars(ax, plot_df):
    tests = planned_tests(plot_df)
    y_min = plot_df["log10_rate"].min()
    y_max = plot_df["log10_rate"].max()
    y_range = max(y_max - y_min, 0.1)
    base = y_max + y_range * 0.08
    spacing = y_range * 0.12
    height = y_range * 0.035

    category_positions = {category: idx for idx, category in enumerate(CATEGORY_ORDER)}
    for idx, (left, right, _stat, p_corrected, n_left, n_right) in enumerate(tests):
        if np.isnan(p_corrected):
            continue
        x1 = category_positions[left]
        x2 = category_positions[right]
        y = base + idx * spacing
        ax.plot(
            [x1, x1, x2, x2],
            [y, y + height, y + height, y],
            color="black",
            linewidth=1.1,
        )
        ax.text(
            (x1 + x2) / 2,
            y + height,
            significance_label(p_corrected),
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )
        print(
            f"  MWU {left} vs {right}: corrected p={p_corrected:.3e}, "
            f"n={n_left}/{n_right}"
        )

    if tests:
        ax.set_ylim(top=base + len(tests) * spacing + y_range * 0.08)


def plot_violin(plot_df, color_guide, output_pdf, output_png):
    sns.set_style("white")
    palette = category_colors(color_guide)

    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    sns.violinplot(
        data=plot_df,
        x="category",
        y="log10_rate",
        hue="category",
        order=CATEGORY_ORDER,
        hue_order=CATEGORY_ORDER,
        palette=palette,
        inner="quartile",
        cut=0,
        density_norm="width",
        linewidth=1.2,
        ax=ax,
        legend=False,
    )

    ax.set_xlabel("")
    ax.set_ylabel("log10 weighted mean recombination rate", fontsize=13)
    ax.set_xticks(np.arange(len(CATEGORY_ORDER)))
    ax.set_xticklabels(CATEGORY_ORDER, rotation=35, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    add_significance_bars(ax, plot_df)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    color_guide = load_color_guide(args.color_guide)
    plot_df = build_plot_data(args)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    plot_df.to_csv(args.output_tsv, sep="\t", index=False)
    print_summary(plot_df)
    plot_violin(plot_df, color_guide, args.output_pdf, args.output_png)
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()

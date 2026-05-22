#!/usr/bin/env python3
"""Plot isoform count and expression by introner count.

This is a parameterized port of the older
visualize_introner_isoform_expression_raincloud.py analysis. It uses
ptitprince raincloud plots when available and falls back to a comparable
violin/box/strip layout otherwise.
"""

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde, mannwhitneyu

warnings.filterwarnings("ignore")

try:
    import ptitprince as pt
except Exception:
    pt = None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True,
                   help="isoform_introner_data_filtered.csv")
    p.add_argument("--output-pdf", required=True)
    p.add_argument("--output-png", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--max-isoforms", type=int, default=50)
    p.add_argument("--expression-zero-frac", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_data(path, max_isoforms):
    data = pd.read_csv(path)
    required = {"gene_id", "strain", "n_isoforms", "introner_count",
                "mean_expression"}
    missing = required - set(data.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    data = data.dropna(subset=list(required)).copy()
    data["introner_count"] = data["introner_count"].astype(int)
    data = data[data["n_isoforms"] <= max_isoforms].copy()
    data["introner_count_binned"] = data["introner_count"].apply(
        lambda x: "3+" if x >= 3 else str(int(x))
    )
    return data


def significance_stars(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def mann_whitney_tests(data, y_col, order):
    rows = []
    for i in range(len(order) - 1):
        left = order[i]
        right = order[i + 1]
        a = data.loc[data["introner_count_binned"] == left, y_col].dropna()
        b = data.loc[data["introner_count_binned"] == right, y_col].dropna()
        if len(a) == 0 or len(b) == 0:
            rows.append((left, right, len(a), len(b), float("nan"), float("nan"), "NA"))
            continue
        u_stat, p_value = mannwhitneyu(a, b, alternative="two-sided")
        rows.append((left, right, len(a), len(b), u_stat, p_value,
                     significance_stars(p_value)))
    return rows


def add_significance_bars(ax, tests, order, base_y, step, log_scale=False):
    for idx, (left, right, _n1, _n2, _u, p_value, stars) in enumerate(tests):
        if pd.isna(p_value):
            continue
        x1 = order.index(left)
        x2 = order.index(right)
        y = base_y * (1.0 + step * idx) if log_scale else base_y + step * idx
        tick = y * 0.03 if log_scale else step * 0.2
        text_y = y * 1.05 if log_scale else y + step * 0.12
        ax.plot([x1, x2], [y, y], color="black", linewidth=1.2, clip_on=False)
        ax.plot([x1, x1], [y, y - tick], color="black", linewidth=1.2,
                clip_on=False)
        ax.plot([x2, x2], [y, y - tick], color="black", linewidth=1.2,
                clip_on=False)
        ax.text((x1 + x2) / 2, text_y, stars, ha="center", va="bottom",
                fontsize=11, clip_on=False)


def subsample_zero_bin(data, fraction, seed):
    zero = data[data["introner_count_binned"] == "0"]
    other = data[data["introner_count_binned"] != "0"]
    if len(zero) == 0 or fraction >= 1:
        return data
    zero = zero.sample(frac=fraction, random_state=seed)
    return pd.concat([zero, other], ignore_index=True)


def kde_profile(values, yscale):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if yscale == "log":
        values = values[values > 0]
        kde_values = np.log10(values)
    else:
        kde_values = values

    if len(kde_values) < 3 or np.nanstd(kde_values) == 0:
        return None, None

    try:
        kde = gaussian_kde(kde_values, bw_method=0.2)
    except Exception:
        return None, None

    grid = np.linspace(kde_values.min(), kde_values.max(), 200)
    density = kde(grid)
    if density.max() == 0:
        return None, None
    y_grid = np.power(10, grid) if yscale == "log" else grid
    return y_grid, density / density.max()


def draw_manual_raincloud(ax, data, y_col, order, hue_order, yscale="linear",
                          seed=42):
    palette = dict(zip(hue_order, sns.color_palette("Set2", len(hue_order))))
    box_offsets = np.linspace(-0.14, 0.14, len(hue_order))
    cloud_offsets = np.linspace(-0.035, 0.035, len(hue_order))
    rng = np.random.default_rng(seed)

    violin_width = 0.28
    box_width = 0.11
    point_jitter = 0.035

    for group_idx, group in enumerate(order):
        for hue_idx, strain in enumerate(hue_order):
            subset = data[
                (data["introner_count_binned"] == group)
                & (data["strain"] == strain)
            ][y_col].dropna()
            if subset.empty:
                continue

            color = palette[strain]
            cloud_edge = group_idx - 0.28 + cloud_offsets[hue_idx]
            box_center = group_idx + box_offsets[hue_idx]

            y_grid, density = kde_profile(subset, yscale)
            if y_grid is not None:
                x_cloud = cloud_edge - density * violin_width
                ax.fill_betweenx(
                    y_grid,
                    x_cloud,
                    cloud_edge,
                    facecolor=color,
                    edgecolor="0.35",
                    linewidth=0.45,
                    alpha=0.42,
                    zorder=1,
                )

            ax.boxplot(
                subset,
                positions=[box_center],
                widths=box_width,
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
                boxprops={
                    "facecolor": color,
                    "edgecolor": "0.35",
                    "linewidth": 0.6,
                    "alpha": 0.55,
                },
                medianprops={"color": "0.35", "linewidth": 0.8},
                whiskerprops={"color": "0.55", "linewidth": 0.7},
                capprops={"color": "0.55", "linewidth": 0.7},
                zorder=3,
            )

            x_points = box_center + rng.normal(0, point_jitter, size=len(subset))
            ax.scatter(
                x_points,
                subset,
                s=4,
                color=color,
                edgecolors="none",
                alpha=0.85,
                zorder=2,
            )

    if yscale == "log":
        ax.set_yscale("log")
    ax.set_xlim(left=-0.7, right=len(order) - 0.7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.grid(False)
    return [Patch(facecolor=palette[strain], edgecolor="0.35", alpha=0.55)
            for strain in hue_order]


def draw_raincloud(ax, data, y_col, order, hue_order, yscale="linear", seed=42):
    if pt is not None:
        pt.RainCloud(
            x="introner_count_binned",
            y=y_col,
            hue="strain",
            data=data,
            palette="Set2",
            bw=0.2,
            width_viol=0.8,
            width_box=0.4,
            orient="v",
            alpha=0.5,
            dodge=True,
            pointplot=False,
            ax=ax,
            order=order,
            hue_order=hue_order,
            box_showfliers=False,
            move=0,
        )
        legend_handles = None
    else:
        legend_handles = draw_manual_raincloud(
            ax, data, y_col, order, hue_order, yscale=yscale, seed=seed
        )

    if yscale == "log":
        ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.grid(False)
    return legend_handles


def clean_legend(ax, hue_order, keep=True, handles_override=None):
    handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    if not keep:
        return
    if handles_override is not None:
        ax.legend(handles_override, hue_order, title="Strain",
                  bbox_to_anchor=(1.01, 1), loc="upper left",
                  fontsize=14, title_fontsize=15)
        return
    selected = []
    selected_labels = []
    for handle, label in zip(handles, labels):
        if label in hue_order and label not in selected_labels:
            selected.append(handle)
            selected_labels.append(label)
    ax.legend(selected, selected_labels, title="Strain", frameon=False,
              bbox_to_anchor=(1.01, 1), loc="upper left")


def make_plot(data, args, stats):
    sns.set_style("white")
    plt.rcParams["figure.dpi"] = 300
    order = ["0", "1", "2", "3+"]
    hue_order = sorted(data["strain"].unique())

    fig, axes = plt.subplots(1, 2, figsize=(24, 8))

    draw_raincloud(axes[0], data, "n_isoforms", order, hue_order,
                   seed=args.seed)
    axes[0].set_xlabel("Introner Count", fontsize=17)
    axes[0].set_ylabel("Number of Isoforms per Gene", fontsize=17)
    axes[0].set_ylim(1, 14)
    axes[0].set_yticks(range(2, 16, 2))
    axes[0].tick_params(labelsize=17)
    add_significance_bars(axes[0], stats["isoform_tests"], order,
                          base_y=10.5, step=1.0)
    clean_legend(axes[0], hue_order, keep=False)

    expression_data = subsample_zero_bin(
        data, args.expression_zero_frac, args.seed
    )
    expression_data = expression_data[expression_data["mean_expression"] > 0].copy()
    legend_handles = draw_raincloud(axes[1], expression_data, "mean_expression",
                                    order, hue_order, yscale="log",
                                    seed=args.seed)
    axes[1].set_xlabel("Introner Count", fontsize=17)
    axes[1].set_ylabel("Log Mean Normalized Expression", fontsize=17)
    axes[1].set_ylim(1, 4000)
    axes[1].tick_params(labelsize=17)
    add_significance_bars(axes[1], stats["expression_tests"], order,
                          base_y=1500, step=0.333, log_scale=True)
    clean_legend(axes[1], hue_order, keep=True,
                 handles_override=legend_handles)

    fig.tight_layout()
    fig.savefig(args.output_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(args.output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(data, stats, args):
    with open(args.report, "w") as fh:
        fh.write("Introner Count vs Isoform Diversity and Gene Expression\n")
        fh.write("=" * 72 + "\n\n")
        fh.write(f"Input: {args.input}\n")
        implementation = "ptitprince" if pt is not None else "manual half-violin fallback"
        fh.write(f"Plot implementation: {implementation}\n")
        fh.write(f"Total observations: {len(data):,}\n")
        fh.write(f"Unique genes: {data['gene_id'].nunique():,}\n")
        fh.write(f"Strains: {', '.join(sorted(data['strain'].unique()))}\n\n")

        fh.write("Introner count distribution\n")
        fh.write("-" * 72 + "\n")
        for group, n in stats["sample_sizes"].items():
            fh.write(f"{group:>3} introners: {n:6d} ({n / len(data) * 100:5.1f}%)\n")
        fh.write("\n")

        fh.write("Summary statistics by introner count\n")
        fh.write("-" * 72 + "\n")
        fh.write(stats["summary"].to_string())
        fh.write("\n\n")

        for title, tests in [
            ("Number of isoforms per gene", stats["isoform_tests"]),
            ("Mean normalized expression", stats["expression_tests"]),
        ]:
            fh.write(title + "\n")
            fh.write("-" * 72 + "\n")
            fh.write(f"{'Comparison':>12} {'N1':>6} {'N2':>6} {'U':>12} {'p':>12} {'Sig':>5}\n")
            for left, right, n1, n2, u_stat, p_value, stars in tests:
                comparison = f"{left} vs {right}"
                fh.write(f"{comparison:>12} {n1:6d} {n2:6d} {u_stat:12.1f} {p_value:12.4e} {stars:>5}\n")
            fh.write("\n")


def main():
    args = parse_args()
    data = load_data(args.input, args.max_isoforms)
    order = ["0", "1", "2", "3+"]
    stats = {
        "sample_sizes": data.groupby("introner_count_binned").size().reindex(order, fill_value=0),
        "summary": data.groupby("introner_count_binned").agg({
            "n_isoforms": ["mean", "median", "std", "min", "max"],
            "mean_expression": ["mean", "median", "std", "min", "max"],
        }).round(3).reindex(order),
        "isoform_tests": mann_whitney_tests(data, "n_isoforms", order),
        "expression_tests": mann_whitney_tests(data, "mean_expression", order),
    }
    make_plot(data, args, stats)
    write_report(data, stats, args)


if __name__ == "__main__":
    main()

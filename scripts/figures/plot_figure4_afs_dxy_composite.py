#!/usr/bin/env python3
"""Render Figure 4 and its separated Population 1 pi supplement.

Main figure:
  A. Population 1 allele-frequency spectra by variant/intron class
  B. Dxy distributions for population-specific and ancestral loci

Supplement:
  Population 1 polymorphic versus fixed introner-body pi
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
EVOLUTION_SCRIPTS = SCRIPTS_DIR / "evolution"
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(EVOLUTION_SCRIPTS))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide
from plot_all_samples_analysis_boxplot import make_violin


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the Figure 4 AFS/Dxy composite and pi supplement."
    )
    parser.add_argument(
        "--afs",
        type=Path,
        default=PROJECT_ROOT
        / "results/snp_popgen/selection/sfs/afs_by_class.tsv",
    )
    parser.add_argument(
        "--flanking-diversity",
        type=Path,
        default=PROJECT_ROOT
        / "results/evolution/diversity_metrics/"
        "all_samples_diversity_metrics_200bp.tsv",
    )
    parser.add_argument(
        "--body-dxy",
        type=Path,
        default=PROJECT_ROOT
        / "results/evolution/diversity_metrics/"
        "shared_introner_body_dxy_200bp.tsv",
    )
    parser.add_argument(
        "--body-decay",
        type=Path,
        default=PROJECT_ROOT
        / "results/evolution/introner_body_decay/"
        "introner_body_decay.per_locus.tsv",
    )
    parser.add_argument(
        "--color-guide",
        type=Path,
        default=DEFAULT_GUIDE_PATH,
    )
    parser.add_argument(
        "--main-output",
        type=Path,
        default=PROJECT_ROOT
        / "results/figures/snp_popgen/figure4_afs_dxy_composite.png",
    )
    parser.add_argument(
        "--pi-output",
        type=Path,
        default=PROJECT_ROOT
        / "results/figures/snp_popgen/"
        "figure4_population1_pi_supplement.png",
    )
    parser.add_argument("--width", type=float, default=6.5)
    parser.add_argument("--height", type=float, default=6.65)
    parser.add_argument("--font-size", type=float, default=9.5)
    parser.add_argument(
        "--layout",
        choices=("stacked", "horizontal"),
        default="stacked",
        help="Main-panel arrangement; horizontal uses a 3:1 AFS:Dxy ratio.",
    )
    parser.add_argument(
        "--horizontal-afs-fraction",
        type=float,
        default=0.75,
        help="Fraction of horizontal panel width assigned to the AFS.",
    )
    parser.add_argument(
        "--horizontal-gap",
        type=float,
        default=0.065,
        help="Figure-width fraction reserved between horizontal panels.",
    )
    parser.add_argument(
        "--skip-pi",
        action="store_true",
        help="Do not regenerate the separate Population 1 pi supplement.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def load_diversity_data(
    flanking_path: Path,
    body_dxy_path: Path,
    body_decay_path: Path,
) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
    flanking_df = pd.read_csv(flanking_path, sep="\t")
    body_dxy_df = pd.read_csv(body_dxy_path, sep="\t")
    body_decay_df = pd.read_csv(body_decay_path, sep="\t")

    shared_flanks = flanking_df[
        flanking_df["category"].eq("group1_fixed_group2_fixed")
    ]
    ancestral_flanks = shared_flanks[
        shared_flanks["cross_group_status"].isin(
            ["ancestral", "likely_ancestral", "ancestral_low_identity"]
        )
    ]
    ancestral_bodies = body_dxy_df[
        body_dxy_df["ancestry_class"].eq("ancestral")
    ]

    dxy_data = {
        "Population 1\nspecific flanks": flanking_df.loc[
            flanking_df["category"].eq("group1_fixed_group2_absent"),
            "dxy_group1_group2",
        ]
        .dropna()
        .tolist(),
        "Population 2\nspecific flanks": flanking_df.loc[
            flanking_df["category"].eq("group1_absent_group2_fixed"),
            "dxy_group1_group2",
        ]
        .dropna()
        .tolist(),
        "Ancestral\nflanks": ancestral_flanks[
            "dxy_group1_group2"
        ]
        .dropna()
        .tolist(),
        "Ancestral\nbody": ancestral_bodies["dxy_introner"].dropna().tolist(),
    }

    group1_bodies = body_decay_df[
        body_decay_df["analysis_scope"].eq("group1_primary")
    ]
    pi_data = {
        "Population 1\npolymorphic body": group1_bodies.loc[
            group1_bodies["analysis_class"].eq("polymorphic"),
            "pi_introner_body",
        ]
        .dropna()
        .tolist(),
        "Population 1\nfixed body": group1_bodies.loc[
            group1_bodies["analysis_class"].eq("fixed_present"),
            "pi_introner_body",
        ]
        .dropna()
        .tolist(),
    }
    return dxy_data, pi_data


def draw_afs(
    ax,
    afs_df: pd.DataFrame,
    colors: dict[str, str],
    font_size: float,
    legend_anchor_x: float = 0.62,
) -> None:
    bins = afs_df["derived_count"].to_numpy()
    width = 0.20
    series = (
        (
            "synonymous_density",
            -1.5,
            "Synonymous",
            colors["synonymous"],
        ),
        (
            "nonsynonymous_density",
            -0.5,
            "Nonsynonymous",
            colors["nonsynonymous"],
        ),
        ("introner_density", 0.5, "Introners", colors["introner"]),
        (
            "non_introner_intron_density",
            1.5,
            "Canonical Introns",
            colors["intron"],
        ),
    )
    for column, offset, label, color in series:
        ax.bar(
            bins + offset * width,
            afs_df[column],
            width,
            label=label,
            color=color,
            edgecolor="#202020",
            linewidth=0.7,
        )

    ax.set_xticks(bins)
    ax.set_xlabel(
        "Population 1 allele count (derived SNPs / present introns)",
        fontsize=font_size + 1,
    )
    ax.set_ylabel("Density", fontsize=font_size + 1)
    ax.tick_params(axis="both", labelsize=font_size, length=3.5, width=0.8)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(legend_anchor_x, 0.985),
        ncol=2,
        frameon=False,
        fontsize=font_size,
        borderaxespad=0.35,
        labelspacing=0.3,
        columnspacing=1.0,
        handlelength=1.7,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)


def draw_main_figure(
    args: argparse.Namespace,
    afs_df: pd.DataFrame,
    dxy_data: dict[str, list[float]],
    guide: dict,
) -> None:
    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")
    if args.layout == "stacked":
        grid = fig.add_gridspec(
            2,
            1,
            height_ratios=(2.15, 3.0),
            left=0.13,
            right=0.975,
            top=0.975,
            bottom=0.085,
            hspace=0.42,
        )
        ax_afs = fig.add_subplot(grid[0])
        ax_dxy = fig.add_subplot(grid[1])
    else:
        if not 0.5 < args.horizontal_afs_fraction < 0.9:
            raise ValueError("--horizontal-afs-fraction must be between 0.5 and 0.9")
        if not 0.03 <= args.horizontal_gap <= 0.18:
            raise ValueError("--horizontal-gap must be between 0.03 and 0.18")
        horizontal_left = 0.105
        horizontal_right = 0.985
        horizontal_gap = args.horizontal_gap
        available_width = horizontal_right - horizontal_left - horizontal_gap
        afs_width = available_width * args.horizontal_afs_fraction
        dxy_width = available_width - afs_width
        dxy_left = horizontal_left + afs_width + horizontal_gap
        ax_afs = fig.add_axes([horizontal_left, 0.25, afs_width, 0.68])
        ax_dxy = fig.add_axes([dxy_left, 0.25, dxy_width, 0.68])

    afs_colors = {
        "synonymous": get_color("Synonymous mutation", guide),
        "nonsynonymous": get_color("Nonsynonymous mutation", guide),
        "introner": get_color("Introner", guide),
        "intron": get_color("Intron", guide),
    }
    draw_afs(
        ax_afs,
        afs_df,
        afs_colors,
        args.font_size,
        legend_anchor_x=0.54 if args.layout == "horizontal" else 0.62,
    )

    dxy_labels = list(dxy_data)
    dxy_colors = [
        get_color("Population 1", guide),
        get_color("Population 2", guide),
        get_color("Ancestor", guide),
        get_color("Ancestor", guide),
    ]
    dxy_comparisons = [
        (
            "Population 1\nspecific flanks",
            "Population 2\nspecific flanks",
        ),
        ("Population 1\nspecific flanks", "Ancestral\nflanks"),
        ("Population 2\nspecific flanks", "Ancestral\nflanks"),
        ("Ancestral\nflanks", "Ancestral\nbody"),
    ]
    make_violin(
        ax_dxy,
        dxy_data,
        dxy_labels,
        dxy_colors,
        "Dxy",
        "",
        comparisons=dxy_comparisons,
        axis_fontsize=args.font_size + 1,
        category_fontsize=args.font_size - 0.5,
        tick_fontsize=args.font_size,
        significance_fontsize=args.font_size,
        full_border=True,
    )
    if args.layout == "horizontal":
        compact_dxy_labels = [
            "P1-specific\nflanks",
            "P2-specific\nflanks",
            "Ancestral\nflanks",
            "Ancestral\nbody",
        ]
        ax_dxy.set_xticklabels(
            compact_dxy_labels,
            rotation=90,
            ha="right",
            va="center",
            rotation_mode="anchor",
            fontsize=args.font_size,
        )
        ax_dxy.tick_params(axis="x", pad=4)
    else:
        ax_dxy.set_xticklabels(
            dxy_labels,
            rotation=0,
            ha="center",
            fontsize=args.font_size - 0.5,
        )
    ax_dxy.tick_params(axis="x", length=3.5, width=0.8)
    for spine in ax_dxy.spines.values():
        spine.set_linewidth(0.8)

    fig.canvas.draw()
    for label, ax in zip(("A", "B"), (ax_afs, ax_dxy), strict=True):
        position = ax.get_position()
        if args.layout == "stacked":
            label_x = 0.018
        elif label == "B":
            label_x = (
                ax_afs.get_position().x1 + position.x0
            ) / 2 - 0.012
        else:
            label_x = position.x0 - (args.horizontal_gap / 2 + 0.012)
        fig.text(
            label_x,
            position.y1,
            label,
            ha="left",
            va="top",
            fontsize=args.font_size + 3,
            fontweight="bold",
        )

    args.main_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.main_output,
        dpi=args.dpi,
        facecolor="white",
        bbox_inches=None,
    )
    plt.close(fig)


def draw_pi_supplement(
    args: argparse.Namespace,
    pi_data: dict[str, list[float]],
    guide: dict,
) -> None:
    fig, ax = plt.subplots(figsize=(4.35, 3.7), facecolor="white")
    labels = list(pi_data)
    colors = [
        get_color("Population 1 Polymorphic", guide),
        get_color("Population 1 Fixed", guide),
    ]
    comparisons = [
        (
            "Population 1\npolymorphic body",
            "Population 1\nfixed body",
        )
    ]
    make_violin(
        ax,
        pi_data,
        labels,
        colors,
        "Pi",
        "",
        comparisons=comparisons,
        axis_fontsize=args.font_size + 1,
        category_fontsize=args.font_size,
        tick_fontsize=args.font_size,
        significance_fontsize=args.font_size,
        full_border=True,
    )
    ax.set_xticklabels(
        labels,
        rotation=0,
        ha="center",
        fontsize=args.font_size,
    )
    ax.tick_params(axis="x", length=3.5, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    fig.subplots_adjust(left=0.18, right=0.975, top=0.965, bottom=0.18)
    args.pi_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.pi_output,
        dpi=args.dpi,
        facecolor="white",
        bbox_inches=None,
    )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": args.font_size,
            "axes.labelsize": args.font_size + 1,
            "xtick.labelsize": args.font_size,
            "ytick.labelsize": args.font_size,
        }
    )

    guide = load_color_guide(args.color_guide)
    afs_df = pd.read_csv(args.afs, sep="\t")
    dxy_data, pi_data = load_diversity_data(
        args.flanking_diversity,
        args.body_dxy,
        args.body_decay,
    )

    print("Dxy sample sizes:")
    for label, values in dxy_data.items():
        print(f"  {label.replace(chr(10), ' ')}: n={len(values)}")
    print("Pi sample sizes:")
    for label, values in pi_data.items():
        print(f"  {label.replace(chr(10), ' ')}: n={len(values)}")

    draw_main_figure(args, afs_df, dxy_data, guide)
    if not args.skip_pi:
        draw_pi_supplement(args, pi_data, guide)
    print(f"Saved main figure: {args.main_output}")
    if not args.skip_pi:
        print(f"Saved pi supplement: {args.pi_output}")


if __name__ == "__main__":
    main()

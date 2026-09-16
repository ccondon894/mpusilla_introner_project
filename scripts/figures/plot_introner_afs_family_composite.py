#!/usr/bin/env python3
"""Render the introner 2D AFS and family distribution on one shared canvas."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GENOTYPING_SCRIPTS = PROJECT_ROOT / "scripts/genotyping"
sys.path.insert(0, str(GENOTYPING_SCRIPTS))

from plot_introner_2d_afs import compute_present_afs


GROUP1_SAMPLES = [
    "CCMP1545",
    "RCC114",
    "RCC1614",
    "RCC1698",
    "RCC2482",
    "RCC373",
    "RCC465",
    "RCC629",
    "RCC692",
    "RCC693",
    "RCC833",
]
GROUP2_SAMPLES = ["RCC1749", "RCC3052"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the introner AFS/family-distribution composite."
    )
    parser.add_argument(
        "--genotype-matrix",
        type=Path,
        default=PROJECT_ROOT / "results/genotyping/genotype_matrix.final.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT
        / "results/figures/genotyping/introner_afs_family_composite.png",
    )
    parser.add_argument("--width", type=float, default=6.5)
    parser.add_argument("--height", type=float, default=4.65)
    parser.add_argument("--font-size", type=float, default=9.5)
    parser.add_argument("--legend-font-size", type=float, default=6.2)
    parser.add_argument(
        "--layout",
        choices=("stacked", "stacked_tall_histogram", "horizontal"),
        default="stacked",
        help=(
            "Panel arrangement; horizontal uses an approximately 3:1 width "
            "ratio, while stacked_tall_histogram doubles the family-panel height."
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def prepare_family_data(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[int], list[str], dict[str, np.ndarray]]:
    filtered = df[
        (df["presence"] != 3)
        & (df["family"] != -1)
        & df["family"].notna()
    ].copy()
    top_families = (
        filtered["family"].value_counts().head(10).index.astype(int).tolist()
    )
    family_df = filtered[filtered["family"].isin(top_families)]

    observed_samples = set(family_df["sample"].unique())
    samples = [sample for sample in GROUP1_SAMPLES if sample in observed_samples]
    samples.extend(
        sample for sample in GROUP2_SAMPLES if sample in observed_samples
    )

    result = {}
    for family in top_families:
        values = []
        for sample in samples:
            values.append(
                int(
                    (
                        (family_df["family"] == family)
                        & (family_df["sample"] == sample)
                        & (family_df["presence"] == 1)
                    ).sum()
                )
            )
        result[family] = values

    plot_data = pd.DataFrame(result, index=samples)

    group1_observed = [sample for sample in samples if sample in GROUP1_SAMPLES]
    group2_observed = [sample for sample in samples if sample in GROUP2_SAMPLES]
    sample_colors = {}
    for sample, color in zip(
        group1_observed,
        cm.viridis(np.linspace(0.2, 0.8, len(group1_observed))),
        strict=True,
    ):
        sample_colors[sample] = color
    group2_colors = ["#B91C1C", "#DC2626"]
    for index, sample in enumerate(group2_observed):
        sample_colors[sample] = group2_colors[index % len(group2_colors)]

    return plot_data, top_families, samples, sample_colors


def draw_afs(
    ax,
    cax,
    spectrum: np.ndarray,
    font_size: float,
    colorbar_label_side: str = "right",
) -> None:
    log_spectrum = np.log10(spectrum + 1)
    image = ax.imshow(
        log_spectrum,
        origin="lower",
        cmap="viridis",
        interpolation="nearest",
        aspect="equal",
    )

    midpoint = (log_spectrum.min() + log_spectrum.max()) / 2
    for row in range(spectrum.shape[0]):
        for column in range(spectrum.shape[1]):
            value = log_spectrum[row, column]
            ax.text(
                column,
                row,
                str(spectrum[row, column]),
                ha="center",
                va="center",
                fontsize=font_size - 0.5,
                color="#202020" if value > midpoint else "white",
            )

    ax.set_xticks(np.arange(spectrum.shape[1]))
    ax.set_yticks(np.arange(spectrum.shape[0]))
    ax.set_xlabel("Population 1\nintroner count", fontsize=font_size + 1)
    ax.set_ylabel("Population 2\nintroner count", fontsize=font_size + 1)
    ax.tick_params(axis="both", labelsize=font_size, length=3.5, width=0.8)

    ax.set_xticks(np.arange(-0.5, spectrum.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, spectrum.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    colorbar = ax.figure.colorbar(image, cax=cax)
    colorbar.set_label(
        r"$\log_{10}$(count + 1)",
        fontsize=font_size + 0.5,
    )
    if colorbar_label_side == "left":
        colorbar.ax.yaxis.set_label_position("left")
        colorbar.ax.yaxis.set_ticks_position("right")
        colorbar.ax.yaxis.labelpad = 2
    elif colorbar_label_side == "compact":
        colorbar.set_label("")
        colorbar.ax.yaxis.set_ticks_position("left")
    colorbar.ax.tick_params(labelsize=font_size - 0.5, length=3)
    colorbar.outline.set_linewidth(0.8)


def draw_family_distribution(
    ax,
    plot_data: pd.DataFrame,
    families: list[int],
    samples: list[str],
    sample_colors: dict[str, np.ndarray],
    font_size: float,
    legend_font_size: float,
    show_legend: bool = True,
) -> None:
    x = np.arange(len(families))
    bar_width = 0.8 / len(samples)

    for index, sample in enumerate(samples):
        positions = x + (index - len(samples) / 2 + 0.5) * bar_width
        ax.bar(
            positions,
            plot_data.loc[sample],
            width=bar_width,
            label=sample,
            color=sample_colors[sample],
            alpha=0.8,
            edgecolor="#202020",
            linewidth=0.55,
        )

    ax.set_yscale("log")
    ax.set_ylim(bottom=1)
    ax.set_xticks(x)
    ax.set_xticklabels(families)
    ax.set_xlabel("Introner Family", fontsize=font_size + 1)
    ax.set_ylabel("Log(Frequency)", fontsize=font_size + 1)
    ax.tick_params(axis="both", which="major", labelsize=font_size, length=3.5)
    ax.tick_params(axis="y", which="minor", length=2)
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)

    if show_legend:
        legend = ax.legend(
            title="Sample",
            fontsize=legend_font_size,
            title_fontsize=legend_font_size + 0.4,
            ncol=3,
            frameon=True,
            edgecolor="#202020",
            fancybox=True,
            loc="upper right",
            borderpad=0.35,
            labelspacing=0.28,
            columnspacing=0.8,
            handlelength=1.7,
            handletextpad=0.45,
        )
        legend.get_frame().set_linewidth(0.8)


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

    print(f"Reading {args.genotype_matrix}")
    genotype_df = pd.read_csv(args.genotype_matrix, sep="\t")
    spectrum, _, used, skipped = compute_present_afs(
        genotype_df,
        GROUP1_SAMPLES,
        GROUP2_SAMPLES,
    )
    plot_data, families, samples, sample_colors = prepare_family_data(genotype_df)
    print(f"AFS loci used: {used:,}")
    print(
        "AFS loci skipped: "
        + ", ".join(f"{reason}={count:,}" for reason, count in skipped.items())
    )

    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")

    if args.layout in {"stacked", "stacked_tall_histogram"}:
        # The AFS has 12 columns and 3 rows. The heatmap dimensions retain
        # square cells; the alternative stacked layout doubles only the family
        # distribution panel height.
        axes_height_inches = 1.24
        axes_height = axes_height_inches / args.height
        afs_width_inches = axes_height_inches * (
            spectrum.shape[1] / spectrum.shape[0]
        )
        afs_width = afs_width_inches / args.width

        left = 0.125
        if args.layout == "stacked_tall_histogram":
            afs_bottom = 0.72
            family_bottom = 0.105
            family_height = 2 * axes_height
        else:
            afs_bottom = 0.675
            family_bottom = 0.19
            family_height = axes_height
        ax_afs = fig.add_axes([left, afs_bottom, afs_width, axes_height])
        cax = fig.add_axes(
            [
                left + afs_width + 0.032,
                afs_bottom,
                0.019,
                axes_height,
            ]
        )
        ax_family = fig.add_axes([left, family_bottom, 0.84, family_height])

        draw_afs(ax_afs, cax, spectrum, args.font_size)
        draw_family_distribution(
            ax_family,
            plot_data,
            families,
            samples,
            sample_colors,
            args.font_size,
            args.legend_font_size,
        )
        panel_positions = [
            (0.018, afs_bottom + axes_height, "A"),
            (0.018, family_bottom + family_height, "B"),
        ]
    else:
        outer_left = 0.06
        outer_right = 0.985
        panel_gap = 0.035
        available_width = outer_right - outer_left - panel_gap
        afs_panel_width = available_width * 0.75
        family_panel_width = available_width * 0.25
        family_left = outer_left + afs_panel_width + panel_gap

        heatmap_left = outer_left + 0.045
        colorbar_width = 0.017
        colorbar_gap = 0.025
        heatmap_right_padding = 0.015
        afs_width = (
            outer_left
            + afs_panel_width
            - heatmap_right_padding
            - colorbar_width
            - colorbar_gap
            - heatmap_left
        )
        afs_height_inches = (
            afs_width
            * args.width
            * spectrum.shape[0]
            / spectrum.shape[1]
        )
        afs_height = afs_height_inches / args.height
        family_bottom = 0.34
        family_height = 0.55
        shared_center = family_bottom + family_height / 2
        afs_bottom = shared_center - afs_height / 2

        ax_afs = fig.add_axes([heatmap_left, afs_bottom, afs_width, afs_height])
        cax = fig.add_axes(
            [
                heatmap_left + afs_width + colorbar_gap,
                afs_bottom,
                colorbar_width,
                afs_height,
            ]
        )
        ax_family = fig.add_axes(
            [family_left, family_bottom, family_panel_width, family_height]
        )

        draw_afs(
            ax_afs,
            cax,
            spectrum,
            args.font_size,
            colorbar_label_side="compact",
        )
        fig.text(
            outer_left + afs_panel_width - 0.006,
            afs_bottom + afs_height + 0.025,
            r"$\log_{10}$(count + 1)",
            ha="right",
            va="bottom",
            fontsize=args.font_size - 0.3,
        )
        draw_family_distribution(
            ax_family,
            plot_data,
            families,
            samples,
            sample_colors,
            args.font_size,
            args.legend_font_size,
            show_legend=False,
        )
        ax_family.set_ylabel("")
        ax_family.set_title(
            "Frequency (log scale)",
            fontsize=args.font_size + 0.3,
            pad=4,
        )
        ax_family.set_xlabel("Introner Family", fontsize=args.font_size + 0.5)

        handles, labels = ax_family.get_legend_handles_labels()
        legend = fig.legend(
            handles,
            labels,
            fontsize=args.legend_font_size,
            ncol=7,
            frameon=True,
            edgecolor="#202020",
            fancybox=True,
            loc="lower center",
            bbox_to_anchor=(0.52, 0.018),
            borderpad=0.3,
            labelspacing=0.25,
            columnspacing=0.7,
            handlelength=1.5,
            handletextpad=0.35,
        )
        legend.get_frame().set_linewidth(0.8)
        panel_positions = [
            (0.018, 0.94, "A"),
            (family_left - 0.026, 0.94, "B"),
        ]

    for x, y, label in panel_positions:
        fig.text(
            x,
            y,
            label,
            ha="left",
            va="top",
            fontsize=args.font_size + 3,
            fontweight="bold",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white", bbox_inches=None)
    plt.close(fig)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

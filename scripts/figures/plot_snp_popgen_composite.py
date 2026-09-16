#!/usr/bin/env python3
"""Render the SNP population-genomics figure on one shared canvas.

Panels:
  A. Folded Population 1 x Population 2 2D SNP AFS
  B. Rooted 4-fold-degenerate-site phylogeny
  C. Best-fitting split-migration demographic model
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from Bio import Phylo
from matplotlib.patches import FancyArrowPatch


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide


GROUP2_SAMPLES = {"RCC1749", "RCC3052"}
OBSOLETE_SAMPLES = {"CCMP490", "RCC647", "RCC835"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the three-panel SNP population-genomics figure."
    )
    parser.add_argument(
        "--vcf",
        default=PROJECT_ROOT
        / "results/snp_popgen/vcf/mpusilla.snps.4d.notMT.vcf.gz",
        type=Path,
    )
    parser.add_argument(
        "--tree",
        default=PROJECT_ROOT
        / "results/snp_popgen/phylogenetics/mpusilla.snps.4d.rooted.treefile",
        type=Path,
    )
    parser.add_argument(
        "--model-fits",
        default=PROJECT_ROOT / "results/snp_popgen/demography/model_fits.4d.txt",
        type=Path,
    )
    parser.add_argument(
        "--target-bed",
        default=PROJECT_ROOT
        / "results/snp_popgen/degenotate/degeneracy-all-sites.4d.bed.gz",
        type=Path,
    )
    parser.add_argument(
        "--color-guide",
        default=DEFAULT_GUIDE_PATH,
        type=Path,
    )
    parser.add_argument(
        "--output",
        default=PROJECT_ROOT
        / "results/figures/snp_popgen/snp_popgen_composite.png",
        type=Path,
    )
    parser.add_argument("--width", type=float, default=6.5)
    parser.add_argument("--height", type=float, default=8.2)
    parser.add_argument("--font-size", type=float, default=9.5)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--mutation-rate", type=float, default=9.8e-10)
    parser.add_argument("--ancestral-time-factor", type=float, default=40.0)
    parser.add_argument("--ancestral-width-factor", type=float, default=3.5)
    return parser.parse_args()


def compute_folded_afs(vcf_path: Path) -> np.ndarray:
    """Match the folding and sample assignment used by the existing AFS panel."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        pop1_indices = [
            idx for idx, sample in enumerate(sample_names) if sample in GROUP2_SAMPLES
        ]
        pop2_indices = [
            idx
            for idx, sample in enumerate(sample_names)
            if sample not in GROUP2_SAMPLES and sample not in OBSOLETE_SAMPLES
        ]
        total_samples = len(pop1_indices) + len(pop2_indices)
        max_minor_count = total_samples // 2
        spectrum = np.zeros(
            (
                min(len(pop1_indices), max_minor_count) + 1,
                max_minor_count + 1,
            ),
            dtype=int,
        )

        for record in vcf:
            if not record.alts or len(record.alts) != 1:
                continue
            pop1_alt = sum(
                1
                for idx in pop1_indices
                if record.samples[idx]["GT"][0]
                and record.samples[idx]["GT"][0] > 0
            )
            pop2_alt = sum(
                1
                for idx in pop2_indices
                if record.samples[idx]["GT"][0]
                and record.samples[idx]["GT"][0] > 0
            )
            total_alt = pop1_alt + pop2_alt
            minor_count = min(total_alt, total_samples - total_alt)
            if minor_count == 0:
                continue
            if total_alt == minor_count:
                pop1_minor, pop2_minor = pop1_alt, pop2_alt
            else:
                pop1_minor = len(pop1_indices) - pop1_alt
                pop2_minor = len(pop2_indices) - pop2_alt
            spectrum[pop1_minor, pop2_minor] += 1

    return spectrum


def draw_afs(ax, cax, spectrum: np.ndarray, font_size: float) -> None:
    log_spectrum = np.log10(spectrum + 1)
    image = ax.imshow(
        log_spectrum,
        origin="lower",
        cmap="viridis",
        interpolation="nearest",
        aspect="equal",
    )

    for row in range(log_spectrum.shape[0]):
        for column in range(log_spectrum.shape[1]):
            value = log_spectrum[row, column]
            ax.text(
                column,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=font_size - 0.5,
                color="#202020" if value >= 4.0 else "white",
            )

    ax.set_xticks(np.arange(log_spectrum.shape[1]))
    ax.set_yticks(np.arange(log_spectrum.shape[0]))
    ax.set_xlabel("Population 1 minor allele count", fontsize=font_size + 1)
    ax.set_ylabel("Population 2\nminor allele count", fontsize=font_size + 1)
    ax.tick_params(axis="both", labelsize=font_size, length=3.5, width=0.8)

    ax.set_xticks(np.arange(-0.5, log_spectrum.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, log_spectrum.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.9)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    colorbar = ax.figure.colorbar(image, cax=cax)
    colorbar.set_label(r"$\log_{10}$(count)", fontsize=font_size + 0.5)
    colorbar.ax.tick_params(labelsize=font_size - 0.5, length=3)
    colorbar.outline.set_visible(False)


def descendant_names(clade) -> set[str]:
    return {tip.name for tip in clade.get_terminals()}


def draw_tree(ax, tree_path: Path, colors: dict[str, str], font_size: float) -> None:
    tree = Phylo.read(str(tree_path), "newick")
    depths = tree.depths()
    max_depth = max(depths[tip] for tip in tree.get_terminals())
    x_positions = {clade: max_depth - depth for clade, depth in depths.items()}

    terminals = tree.get_terminals()
    y_positions = {tip: float(idx) for idx, tip in enumerate(terminals)}

    def assign_internal_y(clade) -> float:
        if clade in y_positions:
            return y_positions[clade]
        child_positions = [assign_internal_y(child) for child in clade.clades]
        y_positions[clade] = (min(child_positions) + max(child_positions)) / 2
        return y_positions[clade]

    assign_internal_y(tree.root)

    def edge_color(clade) -> str:
        names = descendant_names(clade)
        if names and names.issubset(GROUP2_SAMPLES):
            return colors["Population 2"]
        return colors["Population 1"]

    for parent in tree.find_clades(order="preorder"):
        if not parent.clades:
            continue
        parent_x = x_positions[parent]
        parent_y = y_positions[parent]
        for child in parent.clades:
            child_x = x_positions[child]
            child_y = y_positions[child]
            color = edge_color(child)
            ax.plot(
                [parent_x, child_x],
                [child_y, child_y],
                color=color,
                linewidth=1.55,
                solid_capstyle="round",
            )
            ax.plot(
                [parent_x, parent_x],
                [parent_y, child_y],
                color=color,
                linewidth=1.55,
                solid_capstyle="round",
            )

    label_offset = max_depth * 0.035
    for tip in terminals:
        color = (
            colors["Population 2"]
            if tip.name in GROUP2_SAMPLES
            else colors["Population 1"]
        )
        ax.text(
            x_positions[tip] - label_offset,
            y_positions[tip],
            tip.name,
            va="center",
            ha="left",
            fontsize=font_size,
            color=color,
        )

    legend_x = max_depth * 0.99
    legend_y = len(terminals) - 0.55
    for offset, label in enumerate(("Population 1", "Population 2")):
        y = legend_y - offset * 0.62
        ax.scatter(
            [legend_x],
            [y],
            s=28,
            color=colors[label],
            edgecolor="none",
            zorder=3,
        )
        ax.text(
            legend_x - max_depth * 0.02,
            y,
            label,
            va="center",
            ha="left",
            fontsize=font_size,
            color="#333333",
        )

    ax.set_xlim(max_depth * 1.015, -max_depth * 0.16)
    ax.set_ylim(-0.65, len(terminals) + 0.35)
    ticks = np.arange(0, np.floor(max_depth * 10) / 10 + 0.001, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        ["0" if np.isclose(tick, 0) else f"{tick:.1f}" for tick in ticks],
        fontsize=font_size,
    )
    ax.set_xlabel(
        "Substitutions per 4-fold degenerate site",
        fontsize=font_size + 1,
    )
    ax.tick_params(axis="x", length=3.5, width=0.8)
    ax.set_yticks([])
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)


def count_bed_sites(path: Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            total += int(fields[2]) - int(fields[1])
    return total


def best_demographic_parameters(path: Path) -> pd.Series:
    columns = ["log_likelihood", "N1", "N2", "T", "M", "theta"]
    fits = pd.read_csv(path, sep=r"\s+", header=None, names=columns)
    return fits.loc[fits["log_likelihood"].idxmax()]


def scientific_mathtext(value: float) -> str:
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / (10**exponent)
    return rf"{coefficient:.2f}\times10^{{{exponent}}}"


def scaled_tube_half_height(
    size: float,
    max_size: float,
    min_height: float = 0.11,
    max_height: float = 0.32,
) -> float:
    return min_height + (max_height - min_height) * np.sqrt(size / max_size)


def draw_demography(
    ax,
    params: pd.Series,
    target_length: int,
    mutation_rate: float,
    ancestral_time_factor: float,
    ancestral_width_factor: float,
    colors: dict[str, str],
    font_size: float,
) -> None:
    ancestral_ne = params["theta"] / (4 * mutation_rate * target_length)
    population_1_ne = int(params["N1"] * ancestral_ne)
    population_2_ne = int(params["N2"] * ancestral_ne)
    ancestral_plot_size = min(
        ancestral_ne,
        max(population_1_ne, population_2_ne) * ancestral_width_factor,
    )
    split_time = int(params["T"] * 2 * ancestral_ne)
    top_time = int(max(split_time + 1, split_time * ancestral_time_factor))
    migration_rate = params["M"] / (2 * ancestral_ne)

    sizes = {
        "Ancestral": ancestral_plot_size,
        "Population 1": population_1_ne,
        "Population 2": population_2_ne,
    }
    max_size = max(sizes.values())
    half_heights = {
        name: scaled_tube_half_height(size, max_size)
        for name, size in sizes.items()
    }
    centers = {"Ancestral": 0.0, "Population 1": 0.72, "Population 2": -0.72}

    def draw_tube(
        name: str,
        x_start: float,
        x_end: float,
        display_size: float,
    ) -> None:
        center = centers[name]
        half_height = half_heights[name]
        ax.fill_between(
            [x_start, x_end],
            center - half_height,
            center + half_height,
            color=colors[name],
            alpha=0.88,
            linewidth=0.9,
            edgecolor=colors[name],
            zorder=2,
        )
        label_x = np.sqrt(max(x_start, 1) * max(x_end, 1))
        ax.text(
            label_x,
            center,
            f"{name}\n" rf"($N_e={scientific_mathtext(display_size)}$)",
            ha="center",
            va="center",
            fontsize=font_size,
            linespacing=1.05,
            color="#202020",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.2, alpha=0.78),
            zorder=3,
        )

    draw_tube("Ancestral", top_time, split_time, ancestral_ne)
    draw_tube("Population 1", split_time, 1, population_1_ne)
    draw_tube("Population 2", split_time, 1, population_2_ne)

    ancestral_top = centers["Ancestral"] + half_heights["Ancestral"]
    ancestral_bottom = centers["Ancestral"] - half_heights["Ancestral"]
    for name, start_y, end_y in (
        (
            "Population 1",
            ancestral_top,
            centers["Population 1"] - half_heights["Population 1"],
        ),
        (
            "Population 2",
            ancestral_bottom,
            centers["Population 2"] + half_heights["Population 2"],
        ),
    ):
        ax.add_patch(
            FancyArrowPatch(
                (split_time, start_y),
                (split_time, end_y),
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=1.0,
                color=colors[name],
                zorder=4,
            )
        )

    migration_time = np.sqrt(split_time)
    arrows = (
        (
            migration_time * 1.18,
            centers["Population 1"] - 0.18,
            centers["Population 2"] + 0.18,
            colors["Population 1"],
        ),
        (
            migration_time / 1.18,
            centers["Population 2"] + 0.18,
            centers["Population 1"] - 0.18,
            colors["Population 2"],
        ),
    )
    for x, start_y, end_y, color in arrows:
        ax.add_patch(
            FancyArrowPatch(
                (x, start_y),
                (x, end_y),
                arrowstyle="-|>",
                mutation_scale=7,
                linewidth=0.9,
                color=color,
                zorder=4,
            )
        )
    ax.text(
        migration_time,
        0,
        "symmetric migration\n"
        rf"($\lambda={scientific_mathtext(migration_rate)}$ generation$^{{-1}}$)",
        ha="center",
        va="center",
        fontsize=font_size - 1,
        color="#444444",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.2, alpha=0.9),
        zorder=5,
    )

    ax.set_xscale("log")
    ax.set_xlim(top_time * 1.08, 1)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel("Generations ago", fontsize=font_size + 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", which="major", labelsize=font_size, length=3.5)
    ax.tick_params(axis="x", which="minor", length=2)
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)


def main() -> None:
    args = parse_args()
    guide = load_color_guide(args.color_guide)
    colors = {
        "Ancestral": get_color("Ancestor", guide),
        "Population 1": get_color("Population 1", guide),
        "Population 2": get_color("Population 2", guide),
    }

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": args.font_size,
            "axes.labelsize": args.font_size + 1,
            "xtick.labelsize": args.font_size,
            "ytick.labelsize": args.font_size,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    print("Computing folded 2D SNP AFS...")
    spectrum = compute_folded_afs(args.vcf)
    print(f"AFS shape: {spectrum.shape}; segregating sites: {spectrum.sum():,}")

    print("Loading demographic parameters...")
    params = best_demographic_parameters(args.model_fits)
    target_length = count_bed_sites(args.target_bed)
    print(f"4D target length: {target_length:,} bp")

    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=(2.45, 3.1, 2.15),
        left=0.135,
        right=0.965,
        top=0.985,
        bottom=0.065,
        hspace=0.38,
    )
    top = outer[0].subgridspec(1, 2, width_ratios=(1, 0.035), wspace=0.08)
    ax_afs = fig.add_subplot(top[0, 0])
    cax_afs = fig.add_subplot(top[0, 1])
    ax_tree = fig.add_subplot(outer[1])
    ax_demography = fig.add_subplot(outer[2])

    draw_afs(ax_afs, cax_afs, spectrum, args.font_size)
    draw_tree(ax_tree, args.tree, colors, args.font_size)
    draw_demography(
        ax_demography,
        params,
        target_length,
        args.mutation_rate,
        args.ancestral_time_factor,
        args.ancestral_width_factor,
        colors,
        args.font_size,
    )

    fig.canvas.draw()
    for label, ax in zip(
        ("A", "B", "C"),
        (ax_afs, ax_tree, ax_demography),
        strict=True,
    ):
        position = ax.get_position()
        fig.text(
            0.018,
            position.y1,
            label,
            ha="left",
            va="top",
            fontsize=args.font_size + 3,
            fontweight="bold",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.output,
        dpi=args.dpi,
        facecolor="white",
        bbox_inches=None,
    )
    plt.close(fig)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

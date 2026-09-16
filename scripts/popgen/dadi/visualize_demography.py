#!/usr/bin/env python3

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, TextArea, VPacker
from matplotlib.patches import FancyArrowPatch
import demes

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide


def parse_best_fit(params_file):
    """Parse the best-fit parameters from model fits file"""
    # Read the results file
    results = pd.read_csv(params_file, sep='\t', header=None)
    results.columns = ['log_likelihood', 'N1', 'N2', 'T', 'M', 'theta']

    # Find the best fit (highest log-likelihood)
    best_idx = results['log_likelihood'].idxmax()
    best_fit = results.iloc[best_idx]

    print("Best fit parameters:")
    print(f"Log-likelihood: {best_fit['log_likelihood']:.2f}")
    print(f"N1 (intronerful): {best_fit['N1']:.3f}")
    print(f"N2 (intronerless): {best_fit['N2']:.3f}")
    print(f"T (split time): {best_fit['T']:.3f}")
    print(f"M (migration): {best_fit['M']:.6f}")
    print(f"Theta: {best_fit['theta']:.0f}")

    return best_fit


def count_bed_sites(path):
    """Count single-base target sites represented by a BED/BED.gz file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            total += int(fields[2]) - int(fields[1])
    return total


def create_demographic_model(
    params,
    mutation_rate,
    target_length,
    ancestral_time_factor,
    ancestral_width_factor,
):
    """Create a demes demographic model from the best-fit parameters"""
    
    # Extract parameters
    N1 = params['N1']  # Relative size of intronerful population
    N2 = params['N2']  # Relative size of intronerless population  
    T = params['T']    # Split time in 2*Ne generations
    M = params['M']    # Scaled symmetric migration rate, 2*Nanc*m
    theta = params['theta']
    
    # Convert theta to Ne estimate
    # theta = 4*Nanc*mu*L
    Ne_ancestral = theta / (4 * mutation_rate * target_length)
    
    # Convert relative population sizes to absolute sizes
    N1_abs = int(N1 * Ne_ancestral)
    N2_abs = int(N2 * Ne_ancestral)
    ancestral_plot_size = min(
        Ne_ancestral,
        max(N1_abs, N2_abs) * ancestral_width_factor,
    )
    
    # Convert split time from 2*Ne units to generations
    split_time_gens = int(T * 2 * Ne_ancestral)
    migration_per_generation = M / (2 * Ne_ancestral)
    ancestral_start_time = int(max(split_time_gens + 1, split_time_gens * ancestral_time_factor))
    
    print(f"\nDemographic model parameters:")
    print(f"Mutation rate: {mutation_rate:.3e} mutations/site/generation")
    print(f"Mutation target length: {target_length:,} sites")
    print(f"Ancestral Ne: {Ne_ancestral:,.0f}")
    print(f"Group 1 Ne: {N1_abs:,.0f}")
    print(f"Group 2 Ne: {N2_abs:,.0f}")
    print(f"Ancestral tube plotting size: {ancestral_plot_size:,.0f}")
    print(f"Split time: {split_time_gens:,.0f} generations ago")
    print(f"Plot ancestral start time: {ancestral_start_time:,.0f} generations ago")
    print(f"Scaled migration M: {M:.6f}")
    print(f"Migration rate: {migration_per_generation:.3e} per generation")
    
    # Create the demes model
    # Note: demes uses time going backwards, so we need to set times appropriately
    b = demes.Builder(
        description="M. pusilla split-migration demographic model",
        time_units="generations"
    )
    
    # Add ancestral population
    b.add_deme(
        name="ancestral",
        epochs=[
            dict(
                end_time=split_time_gens,
                start_size=ancestral_plot_size,
            )
        ]
    )
    
    # Add daughter populations after split
    b.add_deme(
        name="Population_1", 
        ancestors=["ancestral"],
        epochs=[
            dict(end_time=0, start_size=N1_abs)
        ]
    )
    
    b.add_deme(
        name="Population_2",
        ancestors=["ancestral"], 
        epochs=[
            dict(end_time=0, start_size=N2_abs)
        ]
    )
    
    # Add symmetric migration between populations
    if M > 0:
        b.add_migration(
            demes=["Population_1", "Population_2"],
            rate=migration_per_generation
        )
    
    graph = b.resolve()
    graph.metadata["plot_ancestral_start_time"] = ancestral_start_time
    graph.metadata["true_ancestral_ne"] = Ne_ancestral
    graph.metadata["ancestral_plot_size"] = ancestral_plot_size
    return graph


def deme_plot_size(graph, name):
    deme = next(deme for deme in graph.demes if deme.name == name)
    return deme.epochs[0].start_size


def scaled_tube_half_height(size, max_size, min_height=0.11, max_height=0.32):
    """Map population size to a compact, readable horizontal tube height."""
    return min_height + (max_height - min_height) * np.sqrt(size / max_size)


def scientific_mathtext(value):
    """Format a positive value compactly for a Matplotlib math-text label."""
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / (10 ** exponent)
    return rf"{coefficient:.2f}\times10^{{{exponent}}}"


def add_migration_arrows(
    ax,
    split_time,
    upper_y,
    lower_y,
    migration_rate,
    colors,
    font_size,
):
    """Show symmetric post-split migration without crowding the short panel."""
    migration_time = np.sqrt(split_time)
    upper_edge = upper_y - 0.18
    lower_edge = lower_y + 0.18

    arrows = (
        (migration_time * 1.18, upper_edge, lower_edge, colors["Population_1"]),
        (migration_time / 1.18, lower_edge, upper_edge, colors["Population_2"]),
    )
    for x, start_y, end_y, color in arrows:
        ax.add_patch(
            FancyArrowPatch(
                (x, start_y),
                (x, end_y),
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=1.0,
                color=color,
                alpha=0.9,
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
        fontsize=font_size - 1.5,
        color="#444444",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.9),
        zorder=5,
    )


def visualize_demography(
    graph,
    params,
    output_pdf,
    output_png,
    color_guide_path,
    figure_width,
    figure_height,
    font_size,
):
    """Create a compact, left-to-right demographic model schematic."""

    fig, ax = plt.subplots(1, 1, figsize=(figure_width, figure_height))

    guide = load_color_guide(color_guide_path)

    # Define custom colors for each population
    colors = {
        "ancestral": get_color("Ancestor", guide),
        "Population_1": get_color("Population 1", guide),
        "Population_2": get_color("Population 2", guide),
    }

    split_time = next(
        deme.start_time for deme in graph.demes if deme.name == "Population_1"
    )
    top_time = graph.metadata.get("plot_ancestral_start_time", split_time * 10)
    sizes = {name: deme_plot_size(graph, name) for name in colors}
    max_size = max(sizes.values())
    half_heights = {
        name: scaled_tube_half_height(size, max_size)
        for name, size in sizes.items()
    }
    centers = {"ancestral": 0.0, "Population_1": 0.72, "Population_2": -0.72}

    def draw_tube(name, x_start, x_end, label, parameter_label):
        center = centers[name]
        half_height = half_heights[name]
        ax.fill_between(
            [x_start, x_end],
            center - half_height,
            center + half_height,
            color=colors[name],
            alpha=0.88,
            linewidth=1.0,
            edgecolor=colors[name],
            zorder=2,
        )
        label_x = np.sqrt(max(x_start, 1) * max(x_end, 1))
        label_box = VPacker(
            children=[
                TextArea(
                    label,
                    textprops=dict(fontsize=font_size, color="#202020"),
                ),
                TextArea(
                    parameter_label,
                    textprops=dict(fontsize=font_size - 1, color="#202020"),
                ),
            ],
            align="center",
            pad=0,
            sep=0,
        )
        ax.add_artist(
            AnnotationBbox(
                label_box,
                (label_x, center),
                xycoords="data",
                frameon=True,
                box_alignment=(0.5, 0.5),
                bboxprops=dict(
                    facecolor="white",
                    edgecolor="none",
                    boxstyle="square,pad=0.12",
                    alpha=0.78,
                ),
                zorder=3,
            )
        )

    ancestral_ne = graph.metadata["true_ancestral_ne"]
    population_1_ne = deme_plot_size(graph, "Population_1")
    population_2_ne = deme_plot_size(graph, "Population_2")
    draw_tube(
        "ancestral",
        top_time,
        split_time,
        "Ancestral",
        rf"($N_e={scientific_mathtext(ancestral_ne)}$)",
    )
    draw_tube(
        "Population_1",
        split_time,
        1,
        "Population 1",
        rf"($N_e={scientific_mathtext(population_1_ne)}$)",
    )
    draw_tube(
        "Population_2",
        split_time,
        1,
        "Population 2",
        rf"($N_e={scientific_mathtext(population_2_ne)}$)",
    )

    # Vertical arrows at the split show ancestry without adding diagonal arms
    # or a full-height guide line.
    ancestral_top = centers["ancestral"] + half_heights["ancestral"]
    ancestral_bottom = centers["ancestral"] - half_heights["ancestral"]
    for name in ("Population_1", "Population_2"):
        if name == "Population_1":
            start_y = ancestral_top
            end_y = centers[name] - half_heights[name]
        else:
            start_y = ancestral_bottom
            end_y = centers[name] + half_heights[name]
        ax.add_patch(
            FancyArrowPatch(
                (split_time, start_y),
                (split_time, end_y),
                arrowstyle="-|>",
                mutation_scale=9,
                linewidth=1.2,
                color=colors[name],
                zorder=4,
            )
        )

    if params["M"] > 0:
        add_migration_arrows(
            ax,
            split_time,
            centers["Population_1"],
            centers["Population_2"],
            graph.migrations[0].rate,
            colors,
            font_size,
        )

    ax.set_xscale("log")
    ax.set_xlim(top_time * 1.08, 1)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel("Generations ago", fontsize=font_size + 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", which="major", labelsize=font_size - 1, length=3)
    ax.tick_params(axis="x", which="minor", length=2)
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout(pad=0.45)

    # Save the plot
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(output_png, dpi=300, bbox_inches="tight", pad_inches=0.03)


    return fig

def main():
    """Main function to create demographic visualization"""
    parser = argparse.ArgumentParser(description="Visualize split-migration demographic model")
    parser.add_argument("--params", required=True, help="Model fits TSV file")
    parser.add_argument("--output_pdf", required=True, help="Output plot file (pdf)")
    parser.add_argument("--output_png", required=True, help="Output plot file (png)")
    parser.add_argument(
        "--color-guide",
        default=str(DEFAULT_GUIDE_PATH),
        help="Master figure color guide TSV. Default: %(default)s",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=9.8e-10,
        help="Mutation rate per site per generation. Default: %(default)s",
    )
    parser.add_argument(
        "--ancestral-time-factor",
        type=float,
        default=40.0,
        help=(
            "How far above the fitted split time to start the ancestral deme "
            "for plotting. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--ancestral-width-factor",
        type=float,
        default=3.5,
        help=(
            "Plotting-only cap for ancestral tube width, as a multiple of the "
            "larger daughter population size. Set high to show the fitted "
            "ancestral Ne at full width. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--figure-width",
        type=float,
        default=6.5,
        help="Figure width in inches. Default: %(default)s",
    )
    parser.add_argument(
        "--figure-height",
        type=float,
        default=2.35,
        help="Figure height in inches. Default: %(default)s",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=10.0,
        help="Base font size in points. Default: %(default)s",
    )
    target = parser.add_mutually_exclusive_group()
    target.add_argument(
        "--target-length",
        type=int,
        help="Number of neutral target sites used for theta scaling.",
    )
    target.add_argument(
        "--target-bed",
        help="BED/BED.gz of neutral target sites used for theta scaling.",
    )
    args = parser.parse_args()

    print("Parsing best-fit parameters...")
    best_params = parse_best_fit(args.params)
    if args.target_bed:
        target_length = count_bed_sites(args.target_bed)
    elif args.target_length:
        target_length = args.target_length
    else:
        target_length = 1
        print(
            "\nWARNING: No --target-length or --target-bed supplied. "
            "Using L=1, which will inflate Ne and split time."
        )

    print("\nCreating demographic model...")
    graph = create_demographic_model(
        best_params,
        args.mutation_rate,
        target_length,
        args.ancestral_time_factor,
        args.ancestral_width_factor,
    )

    print("\nGenerating visualization...")
    fig = visualize_demography(
        graph,
        best_params,
        args.output_pdf,
        args.output_png,
        args.color_guide,
        args.figure_width,
        args.figure_height,
        args.font_size,
    )

    print(f"\nVisualization saved: {args.output_png}")


if __name__ == "__main__":
    main()

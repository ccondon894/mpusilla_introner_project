#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import demes

try:
    import demesdraw
except ImportError:
    demesdraw = None

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


def create_demographic_model(params, mutation_rate, target_length, ancestral_time_factor):
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
                start_size=Ne_ancestral,
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
    return graph

def visualize_demography(graph, params, output_pdf, output_png):
    """Create visualization of the demographic model"""

    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # Define custom colors for each population
    colors = {
        'ancestral': "#2e8b57",
        'Population_1': '#1f77b4',
        'Population_2': "#ff8c00"
    }

    split_time = next(deme.start_time for deme in graph.demes if deme.name == "Population_1")
    top_time = graph.metadata.get("plot_ancestral_start_time", split_time * 10)
    if demesdraw is not None:
        demesdraw.tubes(
            graph, 
            ax=ax, 
            log_time=True, 
            colours=colors, 
            labels="xticks-mid",
            num_lines_per_migration=2
        )
        ax.set_ylim(bottom=1, top=top_time)

        # Give demesdraw's layout extra horizontal breathing room.
        left, right = ax.get_xlim()
        pad = 0.05 * (right - left)
        ax.set_xlim(left - pad, right + pad)
    else:
        # Fallback for local environments without demesdraw. Snakemake's
        # moments env includes demesdraw and will use the renderer above.
        tube_width = 0.34
        pop1_x, anc_x, pop2_x = 0.0, 1.0, 2.0

        def draw_tube(x, y0, y1, color, label):
            ax.fill_betweenx(
                [y0, y1],
                x - tube_width,
                x + tube_width,
                color=color,
                alpha=0.88,
                linewidth=0,
            )
            y_label = 10 ** ((np.log10(max(y0, 1)) + np.log10(y1)) / 2)
            ax.text(x, y_label, label, ha="center", va="center", fontsize=13)

        draw_tube(anc_x, split_time, top_time, colors["ancestral"], "ancestral")
        draw_tube(pop1_x, 1, split_time, colors["Population_1"], "intronerful")
        draw_tube(pop2_x, 1, split_time, colors["Population_2"], "intronerless")
        ax.plot(
            [pop1_x, anc_x, pop2_x],
            [split_time, split_time, split_time],
            color="0.2",
            linewidth=1.4,
        )
        ax.annotate(
            "",
            xy=(pop2_x - tube_width, split_time / 30),
            xytext=(pop1_x + tube_width, split_time / 30),
            arrowprops=dict(arrowstyle="<->", color="0.25", lw=1.2),
        )
        ax.set_xticks([])
        ax.set_xlim(-0.7, 2.7)
        ax.set_yscale("log")
        ax.set_ylim(bottom=1, top=top_time)

    ax.set_ylabel('Generations Ago', fontsize=16)

    # Increase tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Increase population label sizes (labels are now mid-positioned on tubes)
    for text in ax.texts:
        text.set_fontsize(14)

    plt.tight_layout()

    # Save the plot
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')


    return fig

def main():
    """Main function to create demographic visualization"""
    parser = argparse.ArgumentParser(description="Visualize split-migration demographic model")
    parser.add_argument("--params", required=True, help="Model fits TSV file")
    parser.add_argument("--output_pdf", required=True, help="Output plot file (pdf)")
    parser.add_argument("--output_png", required=True, help="Output plot file (png)")
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=9.8e-10,
        help="Mutation rate per site per generation. Default: %(default)s",
    )
    parser.add_argument(
        "--ancestral-time-factor",
        type=float,
        default=10.0,
        help=(
            "How far above the fitted split time to start the ancestral deme "
            "for plotting. Default: %(default)s"
        ),
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
    )

    print("\nGenerating visualization...")
    fig = visualize_demography(graph, best_params, args.output_pdf, args.output_png)

    print(f"\nVisualization saved: {args.output_png}")


if __name__ == "__main__":
    main()

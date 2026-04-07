#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import demes
import demesdraw

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

def create_demographic_model(params):
    """Create a demes demographic model from the best-fit parameters"""
    
    # Extract parameters
    N1 = params['N1']  # Relative size of intronerful population
    N2 = params['N2']  # Relative size of intronerless population  
    T = params['T']    # Split time in 2*Ne generations
    M = params['M']    # Migration rate
    theta = params['theta']
    
    # Convert theta to Ne estimate
    # theta = 4*Ne*mu
    mu = 8.2e-10
    Ne_ancestral = theta / (4 * mu)
    
    # Convert relative population sizes to absolute sizes
    N1_abs = int(N1 * Ne_ancestral)
    N2_abs = int(N2 * Ne_ancestral)
    
    # Convert split time from 2*Ne units to generations
    split_time_gens = int(T * 2 * Ne_ancestral)
    
    print(f"\nDemographic model parameters:")
    print(f"Ancestral Ne: {Ne_ancestral:,.0f}")
    print(f"Group 1 Ne: {N1_abs:,.0f}")
    print(f"Group 2 Ne: {N2_abs:,.0f}")
    print(f"Split time: {split_time_gens:,.0f} generations ago")
    print(f"Migration rate: {M:.6f} per generation")
    
    # Create the demes model
    # Note: demes uses time going backwards, so we need to set times appropriately
    b = demes.Builder(
        description="M. pusilla split-migration model from dadi analysis",
        time_units="generations"
    )
    
    # Add ancestral population
    b.add_deme(
        name="ancestral",
        epochs=[
            dict(end_time=split_time_gens, start_size=Ne_ancestral)
        ]
    )
    
    # Add daughter populations after split
    b.add_deme(
        name="Group_1", 
        ancestors=["ancestral"],
        epochs=[
            dict(end_time=0, start_size=N1_abs)
        ]
    )
    
    b.add_deme(
        name="Group_2",
        ancestors=["ancestral"], 
        epochs=[
            dict(end_time=0, start_size=N2_abs)
        ]
    )
    
    # Add symmetric migration between populations
    if M > 0:
        b.add_migration(
            demes=["Group_1", "Group_2"],
            rate=M
        )
    
    return b.resolve()

def visualize_demography(graph, params, output_file):
    """Create visualization of the demographic model"""

    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Define custom colors for each population
    colors = {
        'ancestral': "#69bb87",
        'Group_1': '#1f77b4',
        'Group_2': "#bd2d2d"
    }

    # Demographic history graph with log time scale and custom colors
    demesdraw.tubes(graph, ax=ax, log_time=True, colours=colors, labels='mid')
    ax.set_ylabel('Generations Ago', fontsize=14)

    # Increase tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Increase population label sizes (labels are now mid-positioned on tubes)
    for text in ax.texts:
        text.set_fontsize(14)

    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

    return fig

def main():
    """Main function to create demographic visualization"""
    parser = argparse.ArgumentParser(description="Visualize dadi demographic model")
    parser.add_argument("--params", required=True, help="Model fits TSV file")
    parser.add_argument("--output", required=True, help="Output plot file (pdf or png)")
    args = parser.parse_args()

    print("Parsing best-fit parameters...")
    best_params = parse_best_fit(args.params)

    print("\nCreating demographic model...")
    graph = create_demographic_model(best_params)

    print("\nGenerating visualization...")
    fig = visualize_demography(graph, best_params, args.output)

    print(f"\nVisualization saved: {args.output}")


if __name__ == "__main__":
    main()
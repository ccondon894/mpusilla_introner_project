#!/usr/bin/env python3
"""
Generate bar chart showing counts or proportions of present introners
by family across samples.

Adapted from: /scratch1/chris/introner-genotyping-pipeline/scripts/figure_plotting/introner_family_distributions_panel.py
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Plot introner family distributions across samples')
    parser.add_argument('--genotype-matrix', required=True)
    parser.add_argument('--group1-samples', nargs='+', required=True)
    parser.add_argument('--group2-samples', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--proportions', action='store_true',
                       help='Use log frequency instead of absolute counts')
    return parser.parse_args()


def plot_introner_family_distribution(df, outfile, group1_samples, group2_samples,
                                      use_proportions=False):
    """
    Generate bar charts showing counts or proportions of present introners
    by family across samples. Families are on the x-axis with samples grouped
    for each family.
    """
    if 'family' not in df.columns or df['family'].isnull().all():
        print("No family information available. Skipping.")
        return

    # Filter out missing data and absent family
    filtered_df = df[(df['presence'] != 3) & (df['family'] != -1) &
                     (df['family'] != '') & (~df['family'].isnull())]

    if filtered_df.empty:
        print("No valid data after filtering.")
        return

    # Top N families by frequency
    top_n = 10
    top_families = filtered_df['family'].value_counts().head(top_n).index.tolist()
    family_df = filtered_df[filtered_df['family'].isin(top_families)]

    group1_set = set(group1_samples)
    group2_set = set(group2_samples)
    all_samples = set(family_df['sample'].unique())

    # Custom sort: CCMP1545 first, then other G1, then G2
    samples = []
    if 'CCMP1545' in all_samples:
        samples.append('CCMP1545')
    samples.extend(sorted(s for s in (all_samples & group1_set) if s != 'CCMP1545'))
    samples.extend(sorted(s for s in (all_samples & group2_set)))

    # Calculate total present introners per sample
    sample_totals = {}
    for sample in samples:
        sample_data = filtered_df[(filtered_df['sample'] == sample) & (filtered_df['presence'] == 1)]
        sample_totals[sample] = len(sample_data)

    # Calculate values for each family
    result_data = {}
    for family in top_families:
        family_data = family_df[family_df['family'] == family]
        sample_values = []
        for sample in samples:
            sample_family_data = family_data[(family_data['sample'] == sample) &
                                             (family_data['presence'] == 1)]
            present_count = len(sample_family_data)
            if use_proportions:
                value = present_count if present_count > 0 else 1
            else:
                value = present_count
            sample_values.append(value)
        result_data[family] = sample_values

    plot_data = pd.DataFrame(result_data, index=samples)

    # Print proportion statistics
    print("\nIntroner Family Proportions by Sample:")
    print("=" * 60)

    proportion_data = {}
    for family in top_families:
        family_proportions = []
        family_data = family_df[family_df['family'] == family]
        for sample in samples:
            sample_family_data = family_data[(family_data['sample'] == sample) &
                                             (family_data['presence'] == 1)]
            present_count = len(sample_family_data)
            total_present = sample_totals[sample]
            proportion = present_count / total_present if total_present > 0 else 0
            family_proportions.append(proportion)
        proportion_data[family] = family_proportions

    proportion_df = pd.DataFrame(proportion_data, index=samples).T

    print(f"{'Family':<20}", end="")
    for sample in samples:
        print(f"{sample:>10}", end="")
    print()
    print("-" * (20 + 10 * len(samples)))
    for family in top_families:
        print(f"{family:<20}", end="")
        for sample in samples:
            print(f"{proportion_df.loc[family, sample]:>10.4f}", end="")
        print()

    print(f"\nTotal present introners per sample:")
    for sample in samples:
        print(f"{sample}: {sample_totals[sample]}")
    print()

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4.5))

    x = np.arange(len(top_families))
    width = 0.8 / len(samples)

    g1_list = [s for s in samples if s in group1_set]
    g2_list = [s for s in samples if s in group2_set]

    if len(g1_list) > 0:
        viridis_colors = cm.viridis(np.linspace(0.2, 0.8, len(g1_list)))
    else:
        viridis_colors = []

    red_colors = ['#B91C1C', '#DC2626']

    sample_colors = {}
    for i, sample in enumerate(g1_list):
        sample_colors[sample] = viridis_colors[i]
    for i, sample in enumerate(g2_list):
        sample_colors[sample] = red_colors[i % len(red_colors)]

    for i, sample in enumerate(samples):
        positions = x + (i - len(samples)/2 + 0.5) * width
        color = sample_colors[sample]
        ax.bar(positions, plot_data.loc[sample], width=width, label=sample,
               alpha=0.8, color=color, edgecolor='black', linewidth=1.2)

    ax.set_xlabel('Introner Family', fontsize=14)

    if use_proportions:
        ax.set_ylabel('Log(Frequency)', fontsize=14)
        ax.set_yscale('log')
        ax.set_ylim(bottom=1)
    else:
        ax.set_ylabel('Number of Present Introners', fontsize=14)
        ax.set_ylim(bottom=0)

    ax.set_xticks(x)
    ax.set_xticklabels(top_families, rotation=45, ha='right', fontsize=12)
    ax.tick_params(labelsize=12)

    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(title="Sample", fontsize=7, title_fontsize=7, ncol=3,
              frameon=True, edgecolor='black', loc='upper right')

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Plot saved to: {outfile}")


def main():
    args = parse_arguments()
    df = pd.read_csv(args.genotype_matrix, sep='\t')
    df = df[df["family"] != -1]
    plot_introner_family_distribution(
        df, args.output, args.group1_samples, args.group2_samples,
        use_proportions=args.proportions)


if __name__ == "__main__":
    main()

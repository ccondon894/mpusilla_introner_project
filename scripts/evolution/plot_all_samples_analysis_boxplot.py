#!/usr/bin/env python3
"""
Create a combined violin plot for all-samples flank/body diversity analysis.

The figure shows population-specific flanking Dxy and ancestral flanking/body
Dxy on one axis, then polymorphic and Group 1 fixed introner-body pi on a
second axis so the different value ranges remain legible.
"""

import argparse
import os
os.environ.setdefault('MPLCONFIGDIR', '/scratch1/chris/tmp/matplotlib')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create Dxy boxplot panels for all-samples diversity analysis')
    parser.add_argument('--input', required=True,
                       help='All-samples diversity metrics TSV (flanking Dxy)')
    parser.add_argument('--body-dxy', required=True,
                       help='Fixed shared introner body Dxy TSV')
    parser.add_argument('--body-decay', required=True,
                       help='Introner body decay per-locus TSV with body pi')
    parser.add_argument('--output', required=True,
                       help='Output PNG file')
    parser.add_argument('--flank_length', required=True,
                       help='Flanking sequence length for plot title')
    return parser.parse_args()


def add_significance_bar(ax, x1, x2, y, p_value, height_offset=0.02):
    """Add a significance bar between two positions on the plot."""
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'

    ax.plot([x1, x1, x2, x2], [y, y + height_offset, y + height_offset, y],
            linewidth=1.5, color='black')
    ax.text((x1 + x2) / 2, y + height_offset, sig_text,
            ha='center', va='bottom', fontsize=12, fontweight='bold')


def perform_pairwise_tests(data_dict, comparisons, correction='bonferroni'):
    """
    Perform Mann-Whitney U tests for specified pairwise comparisons.

    Returns dict with comparison keys and p-values (Bonferroni-corrected).
    """
    results = {}
    n_comparisons = len(comparisons)

    for cat1, cat2 in comparisons:
        if cat1 in data_dict and cat2 in data_dict:
            data1 = data_dict[cat1]
            data2 = data_dict[cat2]

            if len(data1) > 0 and len(data2) > 0:
                statistic, p_value = mannwhitneyu(data1, data2, alternative='two-sided')

                if correction == 'bonferroni':
                    p_value_corrected = min(p_value * n_comparisons, 1.0)
                else:
                    p_value_corrected = p_value

                results[(cat1, cat2)] = {
                    'statistic': statistic,
                    'p_value': p_value,
                    'p_value_corrected': p_value_corrected,
                    'n1': len(data1),
                    'n2': len(data2)
                }
            else:
                results[(cat1, cat2)] = None

    return results


def print_group_summary(key, vals):
    """Print a compact summary for one plotted group."""
    if vals:
        print(f"  {key}: n={len(vals)}, median={np.median(vals):.4f}, "
              f"mean={np.mean(vals):.4f}")
    else:
        print(f"  {key}: n=0")


def make_boxplot_panel(ax, data_dict, labels, colors, ylabel, comparisons=None):
    """
    Create a single boxplot panel with seaborn styling.

    Args:
        ax: matplotlib axes
        data_dict: ordered dict of {key: [values]}
        labels: list of display labels matching data_dict keys
        colors: list of colors matching data_dict keys
        ylabel: y-axis label
        comparisons: list of (key1, key2) tuples for significance bars
    """
    keys = list(data_dict.keys())
    data_list = [data_dict[k] for k in keys]

    # Build long-form dataframe for seaborn
    plot_records = []
    for key, label, values in zip(keys, labels, data_list):
        for v in values:
            plot_records.append({'category': label, 'value': v})

    plot_df = pd.DataFrame(plot_records)

    if plot_df.empty:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes)
        return

    sns.violinplot(x='category', y='value', hue='category', data=plot_df,
                   palette=dict(zip(labels, colors)), ax=ax, width=0.8,
                   inner='quartile', cut=0, density_norm='width',
                   legend=False, linewidth=1.2)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=13)

    # Mean diamonds
    means = [np.mean(d) if len(d) > 0 else 0 for d in data_list]
    ax.scatter(range(len(data_list)), means, color='black', marker='D', s=100,
              label='Mean', zorder=4, edgecolor='white', linewidth=1)

    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(labelsize=13)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Significance bars
    if comparisons:
        test_results = perform_pairwise_tests(data_dict, comparisons)

        all_vals = [v for d in data_list for v in d]
        if all_vals:
            max_val = max(all_vals)
        else:
            max_val = 1

        y_margin = max_val * 0.12
        bar_spacing = max_val * 0.09

        for bar_idx, (cat1, cat2) in enumerate(comparisons):
            if (cat1, cat2) in test_results and test_results[(cat1, cat2)] is not None:
                result = test_results[(cat1, cat2)]
                pos1 = keys.index(cat1)
                pos2 = keys.index(cat2)
                bar_y = max_val + y_margin + bar_idx * bar_spacing
                add_significance_bar(ax, pos1, pos2, bar_y, result['p_value_corrected'],
                                    height_offset=max_val * 0.025)

                print(f"  {cat1} vs {cat2}: p={result['p_value']:.2e} "
                      f"(corrected={result['p_value_corrected']:.2e}), "
                      f"n1={result['n1']}, n2={result['n2']}")

        # Adjust y-axis
        n_bars = len(comparisons)
        needed_ylim = max_val + y_margin + n_bars * bar_spacing + max_val * 0.08
        ax.set_ylim(bottom=0, top=needed_ylim)

    ax.legend(loc='upper right', fontsize=10, frameon=True, edgecolor='black')


def make_violin(ax, data_dict, labels, colors, ylabel, title):
    """Create one violin plot with jittered points and mean diamonds."""
    records = [
        {'category': label, 'value': value}
        for label, values in data_dict.items()
        for value in values
    ]
    plot_df = pd.DataFrame(records)
    if plot_df.empty:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes)
        return

    palette = dict(zip(labels, colors))
    sns.violinplot(x='category', y='value', hue='category', data=plot_df,
                   order=labels, palette=palette, ax=ax, width=0.82,
                   inner='quartile', cut=0, density_norm='width',
                   legend=False, linewidth=1.2)
    sns.stripplot(x='category', y='value', data=plot_df, order=labels, ax=ax,
                  color='black', alpha=0.22, size=2.2, jitter=0.18, zorder=3)

    means = [plot_df.loc[plot_df['category'] == label, 'value'].mean()
             for label in labels]
    ax.scatter(range(len(labels)), means, color='black', marker='D', s=72,
               label='Mean', zorder=5, edgecolor='white', linewidth=0.9)

    for idx, label in enumerate(labels):
        n = plot_df.loc[plot_df['category'] == label, 'value'].notna().sum()
        ax.text(idx, 0.985, f'n={n:,}', transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=9, rotation=90)

    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.grid(axis='y', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.legend(loc='upper right', fontsize=10, frameon=True, edgecolor='black')


def main():
    args = parse_arguments()

    # Load data
    print(f"Loading flanking diversity metrics from {args.input}")
    flanking_df = pd.read_csv(args.input, sep='\t')

    print(f"Loading fixed shared introner body Dxy from {args.body_dxy}")
    body_df = pd.read_csv(args.body_dxy, sep='\t')

    print(f"Loading introner body pi from {args.body_decay}")
    body_decay_df = pd.read_csv(args.body_decay, sep='\t')

    shared_flank_df = flanking_df[flanking_df['category'] == 'group1_fixed_group2_fixed']
    ancestral_flank_df = shared_flank_df[
        shared_flank_df['cross_group_status'].isin([
            'ancestral', 'likely_ancestral', 'ancestral_low_identity'
        ])
    ]
    ancestral_body_df = body_df[body_df['ancestry_class'] == 'ancestral']

    g1_present_g2_absent = flanking_df[flanking_df['category'] == 'group1_fixed_group2_absent']['dxy_group1_group2'].dropna().tolist()
    g1_absent_g2_present = flanking_df[flanking_df['category'] == 'group1_absent_group2_fixed']['dxy_group1_group2'].dropna().tolist()
    ancestral_flank = ancestral_flank_df['dxy_group1_group2'].dropna().tolist()
    ancestral_body = ancestral_body_df['dxy_introner'].dropna().tolist()

    group1_body_df = body_decay_df[
        body_decay_df['analysis_scope'].eq('group1_primary')
    ]
    polymorphic_body = group1_body_df[
        group1_body_df['analysis_class'].eq('polymorphic')
    ]['pi_introner_body'].dropna().tolist()
    fixed_body = group1_body_df[
        group1_body_df['analysis_class'].eq('fixed_present')
    ]['pi_introner_body'].dropna().tolist()

    dxy_data = {
        'Population 1\nspecific flanks': g1_present_g2_absent,
        'Population 2\nspecific flanks': g1_absent_g2_present,
        'ancestral\nflanks': ancestral_flank,
        'ancestral\nbody': ancestral_body,
    }
    dxy_labels = list(dxy_data)
    dxy_colors = [
        '#3b528b',
        '#B91C1C',
        '#2d6a4f',
        '#1f9e89',
    ]

    pi_data = {
        'polymorphic\nbody': polymorphic_body,
        'Population 1 fixed\nbody': fixed_body,
    }
    pi_labels = list(pi_data)
    pi_colors = [
        '#F58518',
        '#7B2CBF',
    ]

    print("\n=== Dxy plot groups ===")
    for key, vals in dxy_data.items():
        print_group_summary(key, vals)
    print("\n=== Pi plot groups ===")
    for key, vals in pi_data.items():
        print_group_summary(key, vals)

    # ---- Create figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))
    make_violin(ax1, dxy_data, dxy_labels, dxy_colors, 'Dxy',
                f'Flank and ancestral-body Dxy ({args.flank_length} bp flanks)')
    make_violin(ax2, pi_data, pi_labels, pi_colors, 'pi',
                'Introner-body pi')

    plt.tight_layout()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nFigure saved to: {args.output}")


if __name__ == "__main__":
    main()

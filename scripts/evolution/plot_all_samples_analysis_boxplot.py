#!/usr/bin/env python3
"""
Create a 2-panel boxplot figure for all-samples diversity analysis.

Panel A: Flanking region Dxy across fixation categories
  - Group 1 present & Group 2 absent
  - Group 1 absent & Group 2 present
  - All shared loci
  - Fixed shared loci

Panel B: Introner body Dxy for shared loci, split by family concordance
  - All shared, same family
  - All shared, different family
  - Fixed shared, same family
  - Fixed shared, different family
"""

import argparse
import os
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
    parser.add_argument('--shared-dxy', required=True,
                       help='Shared introner divergence TSV (introner body + flank Dxy)')
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

    sns.boxplot(x='category', y='value', hue='category', data=plot_df,
               palette=dict(zip(labels, colors)), ax=ax, width=0.6,
               showfliers=False, fliersize=0, legend=False, linewidth=1.5)

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


def main():
    args = parse_arguments()

    # Load data
    print(f"Loading flanking diversity metrics from {args.input}")
    flanking_df = pd.read_csv(args.input, sep='\t')

    print(f"Loading shared introner Dxy from {args.shared_dxy}")
    shared_df = pd.read_csv(args.shared_dxy, sep='\t')

    # ---- Panel A: Flanking Dxy ----
    print("\n=== Panel A: Flanking Dxy ===")

    # Boxes 1-2 from flanking diversity metrics
    g1_present_g2_absent = flanking_df[flanking_df['category'] == 'group1_fixed_group2_absent']['dxy_group1_group2'].dropna().tolist()
    g1_absent_g2_present = flanking_df[flanking_df['category'] == 'group1_absent_group2_fixed']['dxy_group1_group2'].dropna().tolist()

    # Boxes 3-4 from shared introner Dxy (flanking component)
    all_shared_flank = shared_df['dxy_flank_mean'].dropna().tolist()
    fixed_shared_flank = shared_df[shared_df['category'] == 'fixed_shared']['dxy_flank_mean'].dropna().tolist()

    panel_a_data = {
        'g1_present_g2_absent': g1_present_g2_absent,
        'g1_absent_g2_present': g1_absent_g2_present,
        'all_shared': all_shared_flank,
        'fixed_shared': fixed_shared_flank,
    }
    panel_a_labels = [
        'Group 1 present\nGroup 2 absent',
        'Group 1 absent\nGroup 2 present',
        'All shared',
        'Fixed shared',
    ]
    panel_a_colors = ['#3b528b', '#B91C1C', '#5ec962', '#2d6a4f']

    panel_a_comparisons = [
        ('g1_present_g2_absent', 'g1_absent_g2_present'),
        ('all_shared', 'fixed_shared'),
        ('g1_present_g2_absent', 'all_shared'),
    ]

    for key, vals in panel_a_data.items():
        print(f"  {key}: n={len(vals)}, median={np.median(vals):.4f}" if vals else f"  {key}: n=0")

    # ---- Panel B: Introner Body Dxy ----
    print("\n=== Panel B: Introner Body Dxy ===")

    valid_shared = shared_df.dropna(subset=['dxy_introner'])

    all_same = valid_shared[valid_shared['family_concordance'] == 'same_family']['dxy_introner'].tolist()
    all_diff = valid_shared[valid_shared['family_concordance'] == 'different_family']['dxy_introner'].tolist()
    fixed_same = valid_shared[(valid_shared['category'] == 'fixed_shared') &
                               (valid_shared['family_concordance'] == 'same_family')]['dxy_introner'].tolist()
    fixed_diff = valid_shared[(valid_shared['category'] == 'fixed_shared') &
                               (valid_shared['family_concordance'] == 'different_family')]['dxy_introner'].tolist()

    panel_b_data = {
        'all_same': all_same,
        'all_diff': all_diff,
        'fixed_same': fixed_same,
        'fixed_diff': fixed_diff,
    }
    panel_b_labels = [
        'All shared\nsame family',
        'All shared\ndiff. family',
        'Fixed shared\nsame family',
        'Fixed shared\ndiff. family',
    ]
    panel_b_colors = ['#2ca02c', '#d62728', '#2d6a4f', '#8b0000']

    panel_b_comparisons = [
        ('all_same', 'all_diff'),
        ('fixed_same', 'fixed_diff'),
        ('all_same', 'fixed_same'),
    ]

    for key, vals in panel_b_data.items():
        print(f"  {key}: n={len(vals)}, median={np.median(vals):.4f}" if vals else f"  {key}: n=0")

    # ---- Create figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    print("\nPanel A significance tests:")
    make_boxplot_panel(ax1, panel_a_data, panel_a_labels, panel_a_colors,
                       f'Flanking Dxy ({args.flank_length}bp)', panel_a_comparisons)
    ax1.set_title('A. Flanking Region Divergence', fontsize=16, fontweight='bold', pad=15)

    print("\nPanel B significance tests:")
    make_boxplot_panel(ax2, panel_b_data, panel_b_labels, panel_b_colors,
                       'Introner Body Dxy', panel_b_comparisons)
    ax2.set_title('B. Introner Body Divergence (Shared Loci)', fontsize=16, fontweight='bold', pad=15)

    plt.tight_layout()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nFigure saved to: {args.output}")


if __name__ == "__main__":
    main()

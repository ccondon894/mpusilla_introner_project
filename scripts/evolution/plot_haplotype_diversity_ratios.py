#!/usr/bin/env python3
"""
Create jitter plots of haplotype diversity deltas (π_introner - π_4d)
with sign tests for each frequency bin.

Panel A: Present Haplotype Diversity (PHDR delta)
Panel B: Absent Haplotype Diversity (AHDR delta)

Points below y=0 indicate flanking diversity lower than neutral expectation.
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binomtest


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create jitter plots for haplotype diversity ratios with sign tests')
    parser.add_argument('--diversity_metrics', required=True,
                       help='Input TSV with introner diversity metrics (frequency info)')
    parser.add_argument('--haplotype_ratios', required=True,
                       help='Input TSV with haplotype diversity ratios (PHDR/AHDR)')
    parser.add_argument('--output', required=True,
                       help='Output PNG file for the plot')
    parser.add_argument('--summary_stats', required=True,
                       help='Output TSV file for sign test summary')
    parser.add_argument('--figsize', default='14,7',
                       help='Figure size as width,height (default: 14,7)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for output figure (default: 300)')
    parser.add_argument('--verbose', action='store_true',
                       help='Print verbose output')
    return parser.parse_args()


def load_and_merge_data(diversity_file, ratios_file, verbose=False):
    """Load and merge diversity metrics with haplotype ratios."""
    if verbose:
        print("Loading and merging data...")

    diversity_df = pd.read_csv(diversity_file, sep='\t')
    ratios_df = pd.read_csv(ratios_file, sep='\t')

    merged_df = pd.merge(
        diversity_df[['ortholog_id', 'frequency', 'present_count', 'absent_count']],
        ratios_df[['ortholog_id', 'PHDR', 'AHDR',
                    'pi_present_introner', 'pi_absent_introner',
                    'avg_4D_pi_present', 'avg_4D_pi_absent']],
        on='ortholog_id', how='inner')

    # Compute deltas
    merged_df['delta_present'] = merged_df['pi_present_introner'] - merged_df['avg_4D_pi_present']
    merged_df['delta_absent'] = merged_df['pi_absent_introner'] - merged_df['avg_4D_pi_absent']

    if verbose:
        print(f"  Merged dataset: {len(merged_df)} rows")

    return merged_df


def run_sign_tests(merged_df, verbose=False):
    """Run binomial sign tests for each frequency bin."""
    if verbose:
        print("Running sign tests...")

    results = []

    # PHDR sign tests (by present frequency, freq >= 2 for pi calculation)
    for freq in sorted(merged_df['frequency'].unique()):
        sub = merged_df[merged_df['frequency'] == freq]
        valid = sub['delta_present'].dropna()

        if len(valid) == 0:
            continue

        n_below = (valid < 0).sum()
        n_total = len(valid)
        pct_below = 100 * n_below / n_total

        test = binomtest(n_below, n_total, 0.5, alternative='two-sided')

        phdr_vals = sub['PHDR'].dropna()

        results.append({
            'ratio_type': 'PHDR',
            'frequency': freq,
            'n_total': n_total,
            'n_below_neutral': n_below,
            'pct_below_neutral': round(pct_below, 1),
            'sign_test_p': test.pvalue,
            'median_delta': round(valid.median(), 6),
            'median_ratio': round(phdr_vals.median(), 4) if len(phdr_vals) > 0 else None,
            'mean_ratio': round(phdr_vals.mean(), 4) if len(phdr_vals) > 0 else None,
        })

    # AHDR sign tests (by absent frequency)
    for absent_freq in sorted(merged_df['absent_count'].unique()):
        sub = merged_df[merged_df['absent_count'] == absent_freq]
        valid = sub['delta_absent'].dropna()

        if len(valid) == 0:
            continue

        n_below = (valid < 0).sum()
        n_total = len(valid)
        pct_below = 100 * n_below / n_total

        test = binomtest(n_below, n_total, 0.5, alternative='two-sided')

        ahdr_vals = sub['AHDR'].dropna()

        results.append({
            'ratio_type': 'AHDR',
            'frequency': absent_freq,
            'n_total': n_total,
            'n_below_neutral': n_below,
            'pct_below_neutral': round(pct_below, 1),
            'sign_test_p': test.pvalue,
            'median_delta': round(valid.median(), 6),
            'median_ratio': round(ahdr_vals.median(), 4) if len(ahdr_vals) > 0 else None,
            'mean_ratio': round(ahdr_vals.mean(), 4) if len(ahdr_vals) > 0 else None,
        })

    results_df = pd.DataFrame(results)

    # Bonferroni correction across all tests
    n_tests = len(results_df)
    if n_tests > 0:
        results_df['sign_test_p_corrected'] = (results_df['sign_test_p'] * n_tests).clip(upper=1.0)

    if verbose:
        print(f"  Applied Bonferroni correction across {n_tests} tests")

    return results_df


def sig_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'ns'


def create_jitter_panel(ax, merged_df, freq_col, delta_col, ratio_label, color_rare,
                        color_common, sign_test_df, ratio_type):
    """Create a single jitter plot panel."""
    # Prepare long-form data
    plot_df = merged_df[[freq_col, delta_col]].dropna().copy()
    plot_df.columns = ['frequency', 'delta']
    plot_df['frequency'] = plot_df['frequency'].astype(int)

    frequencies = sorted(plot_df['frequency'].unique())

    if not frequencies:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return

    # Color by rare vs common
    # Rare present = freq 1-2, rare absent = freq 1-2 (of absent count)
    rare_freqs = {1, 2}
    plot_df['color_group'] = plot_df['frequency'].apply(
        lambda f: 'rare' if f in rare_freqs else 'common')

    palette = {'rare': color_rare, 'common': color_common}

    sns.stripplot(x='frequency', y='delta', hue='color_group', data=plot_df,
                  palette=palette, ax=ax, alpha=0.4, size=4, jitter=0.3,
                  dodge=False, legend=False, zorder=2)

    # Neutral baseline
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7, zorder=3)

    # Add sign test annotations
    sub_tests = sign_test_df[sign_test_df['ratio_type'] == ratio_type]
    y_max = plot_df['delta'].quantile(0.98) if len(plot_df) > 0 else 0.01
    y_min = plot_df['delta'].quantile(0.02) if len(plot_df) > 0 else -0.01

    for i, freq in enumerate(frequencies):
        row = sub_tests[sub_tests['frequency'] == freq]
        if len(row) == 0:
            continue
        row = row.iloc[0]

        stars = sig_stars(row['sign_test_p_corrected'])
        pct = row['pct_below_neutral']
        label = f'{pct:.0f}%\n{stars}'

        ax.text(i, y_max * 1.15, label, ha='center', va='bottom', fontsize=8,
                fontweight='bold' if stars != 'ns' else 'normal')

    # Styling
    ax.set_xlabel(f'{ratio_label} Frequency', fontsize=14)
    ax.set_ylabel(f'$\\pi_{{introner}}$ $-$ $\\pi_{{4D}}$', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.2, axis='y')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set y-limits with room for annotations
    ax.set_ylim(y_min * 1.3, y_max * 1.5)

    # Add sample size labels at bottom
    for i, freq in enumerate(frequencies):
        n = len(plot_df[plot_df['frequency'] == freq])
        ax.text(i, y_min * 1.2, f'n={n}', ha='center', va='top', fontsize=8, color='gray')


def main():
    args = parse_arguments()
    figsize = tuple(map(float, args.figsize.split(',')))

    # Load data
    merged_df = load_and_merge_data(args.diversity_metrics, args.haplotype_ratios, args.verbose)

    # Run sign tests
    sign_test_df = run_sign_tests(merged_df, args.verbose)

    # Save sign test summary
    os.makedirs(os.path.dirname(args.summary_stats), exist_ok=True)
    sign_test_df.to_csv(args.summary_stats, sep='\t', index=False)
    print(f"Sign test summary saved to: {args.summary_stats}")

    if args.verbose:
        print(f"\nSign test results (Bonferroni-corrected, n={len(sign_test_df)} tests):")
        for _, row in sign_test_df.iterrows():
            print(f"  {row['ratio_type']} freq {int(row['frequency'])}: "
                  f"{row['pct_below_neutral']:.0f}% below neutral, "
                  f"p={row['sign_test_p']:.2e}, "
                  f"p_corrected={row['sign_test_p_corrected']:.2e} "
                  f"{sig_stars(row['sign_test_p_corrected'])}")

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Panel A: PHDR delta (by present frequency)
    create_jitter_panel(ax1, merged_df, 'frequency', 'delta_present',
                       'Present', '#2166ac', '#92c5de',
                       sign_test_df, 'PHDR')
    ax1.set_title('A. Present Haplotype Diversity', fontsize=16, fontweight='bold', pad=15)

    # Panel B: AHDR delta (by absent frequency)
    create_jitter_panel(ax2, merged_df, 'absent_count', 'delta_absent',
                       'Absent', '#b2182b', '#f4a582',
                       sign_test_df, 'AHDR')
    ax2.set_title('B. Absent Haplotype Diversity', fontsize=16, fontweight='bold', pad=15)

    plt.tight_layout()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    plt.close()

    print(f"Plot saved to: {args.output}")


if __name__ == "__main__":
    main()

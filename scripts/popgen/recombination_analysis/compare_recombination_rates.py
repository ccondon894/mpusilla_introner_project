#!/usr/bin/env python3
"""
Compare recombination rates between exonic regions and introner-containing regions.
Uses pyrho output, GTF annotations, and introner BED file.
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns


def parse_pyrho_file(filepath):
    """Parse a single pyrho output file."""
    df = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        names=['start', 'end', 'rate']
    )
    return df


def parse_gtf_exons(gtf_file):
    """Extract exon coordinates from GTF file."""
    exons = defaultdict(list)

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])

            if feature == 'exon':
                exons[chrom].append((start, end))

    return exons


def parse_bed_file(bed_file):
    """Parse BED file for introner regions."""
    regions = defaultdict(list)

    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            regions[chrom].append((start, end))

    return regions


def get_overlap_length(start1, end1, start2, end2):
    """Calculate overlap length between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def classify_window(chrom, start, end, exons, introners, method='midpoint', exclude_mating_type=True):
    """
    Classify a pyrho window as exonic, introner, or neither.

    Methods:
    - 'midpoint': Use midpoint of window
    - 'majority': Use majority overlap (>50%)
    - 'any': Any overlap counts
    """
    # Exclude mating-type region on scaffold_2
    if exclude_mating_type and chrom == 'scaffold_2':
        midpoint = (start + end) // 2
        if 49808 <= midpoint <= 1730591:
            return 'mating_type_excluded'

    window_len = end - start

    if method == 'midpoint':
        midpoint = (start + end) // 2

        # Check exons
        for ex_start, ex_end in exons.get(chrom, []):
            if ex_start <= midpoint <= ex_end:
                return 'exon'

        # Check introners
        for in_start, in_end in introners.get(chrom, []):
            if in_start <= midpoint <= in_end:
                return 'introner'

        return 'intergenic'

    elif method == 'majority':
        exon_overlap = 0
        introner_overlap = 0

        for ex_start, ex_end in exons.get(chrom, []):
            exon_overlap += get_overlap_length(start, end, ex_start, ex_end)

        for in_start, in_end in introners.get(chrom, []):
            introner_overlap += get_overlap_length(start, end, in_start, in_end)

        if exon_overlap > window_len / 2:
            return 'exon'
        elif introner_overlap > window_len / 2:
            return 'introner'
        else:
            return 'intergenic'

    elif method == 'any':
        # Check any overlap
        for ex_start, ex_end in exons.get(chrom, []):
            if get_overlap_length(start, end, ex_start, ex_end) > 0:
                return 'exon'

        for in_start, in_end in introners.get(chrom, []):
            if get_overlap_length(start, end, in_start, in_end) > 0:
                return 'introner'

        return 'intergenic'

    return 'other'


def weighted_median(values, weights):
    """Calculate weighted median."""
    values = np.array(values)
    weights = np.array(weights)

    sorted_indices = np.argsort(values)
    sorted_values = values[sorted_indices]
    sorted_weights = weights[sorted_indices]

    cumsum = np.cumsum(sorted_weights)
    cutoff = cumsum[-1] / 2.0

    return sorted_values[cumsum >= cutoff][0]


def main():
    parser = argparse.ArgumentParser(
        description="Compare recombination rates: exons vs introner regions"
    )
    parser.add_argument(
        '--pyrho_dir',
        required=True,
        help='Directory with pyrho output files'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='GTF file with exon annotations'
    )
    parser.add_argument(
        '--introner_bed',
        required=True,
        help='BED file with introner regions'
    )
    parser.add_argument(
        '--method',
        default='midpoint',
        choices=['midpoint', 'majority', 'any'],
        help='Classification method (default: midpoint)'
    )
    parser.add_argument(
        '--output_prefix',
        default='recombination_comparison',
        help='Output file prefix'
    )
    args = parser.parse_args()

    print(" Loading annotations...")
    exons = parse_gtf_exons(args.gtf)
    introners = parse_bed_file(args.introner_bed)

    print(f"   Found {sum(len(v) for v in exons.values())} exons")
    print(f"   Found {sum(len(v) for v in introners.values())} introner regions")

    print("\n📊 Processing pyrho files...")

    # Collect rates AND window lengths by category
    data = defaultdict(lambda: {'rates': [], 'lengths': []})

    import os
    import glob

    pyrho_files = glob.glob(os.path.join(args.pyrho_dir, '*.pyrho.out'))

    for pyrho_file in sorted(pyrho_files):
        # Extract chromosome/scaffold name
        basename = os.path.basename(pyrho_file)
        chrom = basename.replace('.pyrho.out', '')

        # Read pyrho data
        df = parse_pyrho_file(pyrho_file)

        # Classify each window
        for _, row in df.iterrows():
            start, end, rate = int(row['start']), int(row['end']), row['rate']
            length = end - start

            category = classify_window(chrom, start, end, exons, introners, args.method)
            data[category]['rates'].append(rate)
            data[category]['lengths'].append(length)

        print(f"   Processed {chrom}: {len(df)} windows")

    # Convert to arrays
    for category in data:
        data[category]['rates'] = np.array(data[category]['rates'])
        data[category]['lengths'] = np.array(data[category]['lengths'])

    # Summary statistics
    print("\n📈 Summary Statistics:")
    print(f"\n{'Category':<25} {'N Windows':<12} {'Total bp':<15} {'Mean Length (bp)':<20}")
    print("-" * 80)

    for category in ['exon', 'introner', 'intergenic', 'mating_type_excluded']:
        n = len(data[category]['rates'])
        total_bp = np.sum(data[category]['lengths'])
        mean_len = np.mean(data[category]['lengths']) if n > 0 else 0
        print(f"{category:<25} {n:<12} {total_bp:<15,} {mean_len:<20.1f}")

    # Convert to arrays for tests
    exon_rates = data['exon']['rates']
    introner_rates = data['introner']['rates']
    exon_lengths = data['exon']['lengths']
    introner_lengths = data['introner']['lengths']

    print("\n📊 Recombination Rates:")
    print(f"   Exon - Unweighted median: {np.median(exon_rates):.2e}")
    print(f"   Exon - Weighted median:   {weighted_median(exon_rates, exon_lengths):.2e}")
    print(f"   Introner - Unweighted median: {np.median(introner_rates):.2e}")
    print(f"   Introner - Weighted median:   {weighted_median(introner_rates, introner_lengths):.2e}")

    # Check for systematic window length differences
    print(f"\n📏 Window Length Comparison:")
    print(f"   Exon windows: mean = {np.mean(exon_lengths):.1f} bp, median = {np.median(exon_lengths):.1f} bp")
    print(f"   Introner windows: mean = {np.mean(introner_lengths):.1f} bp, median = {np.median(introner_lengths):.1f} bp")

    if len(exon_lengths) > 0 and len(introner_lengths) > 0:
        stat_len, pval_len = mannwhitneyu(exon_lengths, introner_lengths, alternative='two-sided')
        print(f"   Mann-Whitney test on lengths: p = {pval_len:.2e}")
        if pval_len < 0.05:
            print(f"   ⚠️  Window lengths differ significantly between categories!")
        else:
            print(f"   ✅ Window lengths are similar between categories")

    # Mann-Whitney U test (UNWEIGHTED)
    if len(exon_rates) > 0 and len(introner_rates) > 0:
        statistic, pvalue = mannwhitneyu(exon_rates, introner_rates, alternative='two-sided')

        print(f"\n🧪 Mann-Whitney U Test (UNWEIGHTED - each window counts once):")
        print(f"   U-statistic: {statistic:.2e}")
        print(f"   p-value: {pvalue:.2e}")

        if pvalue < 0.05:
            print(f"   ✅ SIGNIFICANT difference (p < 0.05)")
        else:
            print(f"   ❌ NO significant difference (p >= 0.05)")

    # Mann-Whitney U test (WEIGHTED by window length)
    print(f"\n⚖️  Mann-Whitney U Test (WEIGHTED - each bp counts once):")

    # Expand data: repeat each rate by its window length
    exon_rates_expanded = np.repeat(exon_rates, exon_lengths.astype(int))
    introner_rates_expanded = np.repeat(introner_rates, introner_lengths.astype(int))

    if len(exon_rates_expanded) > 0 and len(introner_rates_expanded) > 0:
        statistic_w, pvalue_w = mannwhitneyu(exon_rates_expanded, introner_rates_expanded, alternative='two-sided')

        print(f"   U-statistic: {statistic_w:.2e}")
        print(f"   p-value: {pvalue_w:.2e}")

        if pvalue_w < 0.05:
            print(f"   ✅ SIGNIFICANT difference (p < 0.05)")
        else:
            print(f"   ❌ NO significant difference (p >= 0.05)")


    # Save results
    results_df = pd.DataFrame({
        'category': ['exon', 'introner', 'intergenic', 'mating_type_excluded'],
        'n_windows': [
            len(data['exon']['rates']),
            len(data['introner']['rates']),
            len(data['intergenic']['rates']),
            len(data['mating_type_excluded']['rates'])
        ],
        'total_bp': [
            np.sum(data['exon']['lengths']),
            np.sum(data['introner']['lengths']),
            np.sum(data['intergenic']['lengths']),
            np.sum(data['mating_type_excluded']['lengths'])
        ],
        'mean_window_length': [
            np.mean(data['exon']['lengths']) if len(data['exon']['rates']) > 0 else np.nan,
            np.mean(data['introner']['lengths']) if len(data['introner']['rates']) > 0 else np.nan,
            np.mean(data['intergenic']['lengths']) if len(data['intergenic']['rates']) > 0 else np.nan,
            np.mean(data['mating_type_excluded']['lengths']) if len(data['mating_type_excluded']['rates']) > 0 else np.nan
        ],
        'unweighted_median_rate': [
            np.median(data['exon']['rates']) if len(data['exon']['rates']) > 0 else np.nan,
            np.median(data['introner']['rates']) if len(data['introner']['rates']) > 0 else np.nan,
            np.median(data['intergenic']['rates']) if len(data['intergenic']['rates']) > 0 else np.nan,
            np.median(data['mating_type_excluded']['rates']) if len(data['mating_type_excluded']['rates']) > 0 else np.nan
        ],
        'weighted_median_rate': [
            weighted_median(data['exon']['rates'], data['exon']['lengths']) if len(data['exon']['rates']) > 0 else np.nan,
            weighted_median(data['introner']['rates'], data['introner']['lengths']) if len(data['introner']['rates']) > 0 else np.nan,
            weighted_median(data['intergenic']['rates'], data['intergenic']['lengths']) if len(data['intergenic']['rates']) > 0 else np.nan,
            weighted_median(data['mating_type_excluded']['rates'], data['mating_type_excluded']['lengths']) if len(data['mating_type_excluded']['rates']) > 0 else np.nan
        ]
    })

    results_file = f"{args.output_prefix}_summary.tsv"
    results_df.to_csv(results_file, sep='\t', index=False)
    print(f"\n Saved summary to {results_file}")

    # Create visualization
    print("\n📊 Creating plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Unweighted rates
    plot_data = []
    for category in ['exon', 'introner', 'intergenic']:
        for rate in data[category]['rates']:
            plot_data.append({'category': category, 'rate': rate})

    plot_df = pd.DataFrame(plot_data)

    sns.violinplot(
        data=plot_df[plot_df['category'].isin(['exon', 'introner'])],
        x='category',
        y='rate',
        ax=axes[0, 0],
        inner='box'
    )
    axes[0, 0].set_ylabel('Recombination Rate (ρ)')
    axes[0, 0].set_xlabel('Genomic Category')
    axes[0, 0].set_title(f'Unweighted (each window counts once)\np={pvalue:.2e}')
    axes[0, 0].set_yscale('log')

    # Plot 2: Window length distributions
    length_data = []
    for category in ['exon', 'introner']:
        for length in data[category]['lengths']:
            length_data.append({'category': category, 'length': length})

    length_df = pd.DataFrame(length_data)

    sns.violinplot(
        data=length_df,
        x='category',
        y='length',
        ax=axes[0, 1],
        inner='box'
    )
    axes[0, 1].set_ylabel('Window Length (bp)')
    axes[0, 1].set_xlabel('Genomic Category')
    axes[0, 1].set_title(f'Window Length Distribution\np={pval_len:.2e}')
    axes[0, 1].set_yscale('log')

    # Plot 3: Comparison of weighted vs unweighted medians
    weighted_comparison = pd.DataFrame({
        'category': ['exon', 'exon', 'introner', 'introner'],
        'median_rate': [
            np.median(data['exon']['rates']),
            weighted_median(data['exon']['rates'], data['exon']['lengths']),
            np.median(data['introner']['rates']),
            weighted_median(data['introner']['rates'], data['introner']['lengths'])
        ],
        'type': ['Unweighted', 'Weighted', 'Unweighted', 'Weighted']
    })

    sns.barplot(
        data=weighted_comparison,
        x='category',
        y='median_rate',
        hue='type',
        ax=axes[1, 0]
    )
    axes[1, 0].set_ylabel('Median Recombination Rate (ρ)')
    axes[1, 0].set_xlabel('Genomic Category')
    axes[1, 0].set_title('Unweighted vs Weighted Medians')
    axes[1, 0].set_yscale('log')
    axes[1, 0].legend(title='Method')

    # Plot 4: Weighted analysis summary
    axes[1, 1].text(0.1, 0.9, 'Statistical Tests Summary', fontsize=12, fontweight='bold', transform=axes[1, 1].transAxes)

    summary_text = f"""
Unweighted Test (window-level):
  p-value: {pvalue:.2e}
  {'Significant' if pvalue < 0.05 else 'Not significant'}

Weighted Test (bp-level):
  p-value: {pvalue_w:.2e}
  {'Significant' if pvalue_w < 0.05 else 'Not significant'}

Window Lengths:
  Exon mean: {np.mean(exon_lengths):.1f} bp
  Introner mean: {np.mean(introner_lengths):.1f} bp
  Difference p: {pval_len:.2e}
"""

    axes[1, 1].text(0.1, 0.7, summary_text, fontsize=9, transform=axes[1, 1].transAxes,
                    verticalalignment='top', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    axes[1, 1].axis('off')

    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file)
    print(f"💾 Saved plot to {plot_file}")

    # Save detailed data with lengths
    detail_data = []
    for category in ['exon', 'introner', 'intergenic']:
        for rate, length in zip(data[category]['rates'], data[category]['lengths']):
            detail_data.append({'category': category, 'rate': rate, 'window_length': length})

    detail_df = pd.DataFrame(detail_data)
    detail_file = f"{args.output_prefix}_detailed.tsv"
    detail_df.to_csv(detail_file, sep='\t', index=False)
    print(f"💾 Saved detailed data to {detail_file}")


if __name__ == '__main__':
    main()

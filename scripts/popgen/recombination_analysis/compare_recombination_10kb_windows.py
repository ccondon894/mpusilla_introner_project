#!/usr/bin/env python3
"""
Compare recombination rates between introner-containing and non-introner-containing
genomic regions using fixed 10kb windows.

Each 10kb window is classified as:
- introner_containing: overlaps ≥1 introner
- non_introner_containing: overlaps 0 introners
- excluded: overlaps mating-type region

Weighted mean recombination rate is calculated per window from pyrho output.
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob


def parse_pyrho_file(filepath):
    """Parse a single pyrho output file."""
    df = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        names=['start', 'end', 'rate']
    )
    return df


def parse_bed_file(bed_file):
    """Parse BED file for introner regions."""
    regions = defaultdict(list)

    with open(bed_file, 'r') as f:
        for line in f:
            # Handle both tab and space-delimited BED files
            fields = line.strip().split()
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


def calculate_weighted_mean_rate(window_start, window_end, pyrho_windows):
    """
    Calculate weighted mean recombination rate for a window.

    Args:
        window_start: Start of 10kb window
        window_end: End of 10kb window
        pyrho_windows: DataFrame with pyrho output (start, end, rate)

    Returns:
        Weighted mean rate, or np.nan if no overlap
    """
    # Vectorized approach for efficiency
    # Find pyrho windows that overlap with our window
    mask = (pyrho_windows['end'] > window_start) & (pyrho_windows['start'] < window_end)
    overlapping = pyrho_windows[mask]

    if len(overlapping) == 0:
        return np.nan

    # Calculate overlaps
    overlap_starts = np.maximum(overlapping['start'].values, window_start)
    overlap_ends = np.minimum(overlapping['end'].values, window_end)
    overlap_lengths = overlap_ends - overlap_starts

    # Weighted mean
    total_weighted_rate = np.sum(overlapping['rate'].values * overlap_lengths)
    total_length = np.sum(overlap_lengths)

    if total_length == 0:
        return np.nan

    return total_weighted_rate / total_length


def check_introner_overlap(window_start, window_end, introner_regions):
    """Check if window overlaps any introner region."""
    for introner_start, introner_end in introner_regions:
        if get_overlap_length(window_start, window_end, introner_start, introner_end) > 0:
            return True
    return False


def count_introners_in_window(window_start, window_end, introner_regions):
    """Count number of introners that overlap with window."""
    count = 0
    for introner_start, introner_end in introner_regions:
        if get_overlap_length(window_start, window_end, introner_start, introner_end) > 0:
            count += 1
    return count


def classify_window(chrom, window_start, window_end, introner_regions,
                   exclude_mating_type=True):
    """
    Classify a 10kb window.

    Returns:
        'introner_containing', 'non_introner_containing', or 'excluded'
    """
    # Check if overlaps mating-type region
    if exclude_mating_type and chrom == 'scaffold_2':
        # Mating-type region: scaffold_2:49,808-1,730,591
        if get_overlap_length(window_start, window_end, 49808, 1730591) > 0:
            return 'excluded'

    # Check for introner overlap
    if check_introner_overlap(window_start, window_end, introner_regions.get(chrom, [])):
        return 'introner_containing'

    return 'non_introner_containing'


def main():
    parser = argparse.ArgumentParser(
        description="Compare recombination rates using 10kb windows: introner-containing vs non-introner-containing"
    )
    parser.add_argument(
        '--pyrho_dir',
        required=True,
        help='Directory with pyrho output files'
    )
    parser.add_argument(
        '--introner_bed',
        required=True,
        help='BED file with introner regions'
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=10000,
        help='Window size in bp (default: 10000)'
    )
    parser.add_argument(
        '--output_prefix',
        default='recombination_10kb_windows',
        help='Output file prefix'
    )
    args = parser.parse_args()

    print("📂 Loading introner annotations...")
    introners = parse_bed_file(args.introner_bed)
    print(f"   Found {sum(len(v) for v in introners.values())} introner regions")

    print(f"\n📊 Processing pyrho files with {args.window_size}bp windows...")

    # Collect pyrho data by chromosome
    pyrho_data = {}
    pyrho_files = glob.glob(os.path.join(args.pyrho_dir, '*.pyrho.out'))

    for pyrho_file in sorted(pyrho_files):
        basename = os.path.basename(pyrho_file)
        chrom = basename.replace('.pyrho.out', '')
        pyrho_data[chrom] = parse_pyrho_file(pyrho_file)
        print(f"   Loaded {chrom}: {len(pyrho_data[chrom])} pyrho windows")

    print(f"\n🔲 Creating {args.window_size}bp windows and calculating rates...")

    # Store window data
    windows = []

    for chrom in sorted(pyrho_data.keys()):
        pyrho_windows = pyrho_data[chrom]

        # Get chromosome length from pyrho data
        chrom_length = int(pyrho_windows['end'].max())

        # Create 10kb windows
        for window_start in range(0, chrom_length, args.window_size):
            window_end = min(window_start + args.window_size, chrom_length)

            # Calculate weighted mean rate
            mean_rate = calculate_weighted_mean_rate(window_start, window_end, pyrho_windows)

            if np.isnan(mean_rate):
                continue

            # Classify window
            classification = classify_window(chrom, window_start, window_end, introners)

            # Count introners
            n_introners = count_introners_in_window(
                window_start, window_end, introners.get(chrom, [])
            )

            windows.append({
                'scaffold': chrom,
                'start': window_start,
                'end': window_end,
                'weighted_mean_rate': mean_rate,
                'category': classification,
                'n_introners': n_introners
            })

    # Convert to DataFrame
    windows_df = pd.DataFrame(windows)

    print(f"   Created {len(windows_df)} windows total")

    # Summary statistics
    print("\n📈 Summary Statistics:")
    print(f"\n{'Category':<30} {'N Windows':<12} {'Mean Rate':<20} {'Median Rate':<20}")
    print("-" * 85)

    stats_summary = []
    for category in ['non_introner_containing', 'introner_containing', 'excluded']:
        subset = windows_df[windows_df['category'] == category]
        n = len(subset)
        mean_rate = np.mean(subset['weighted_mean_rate']) if n > 0 else np.nan
        median_rate = np.median(subset['weighted_mean_rate']) if n > 0 else np.nan

        print(f"{category:<30} {n:<12} {mean_rate:<20.2e} {median_rate:<20.2e}")

        stats_summary.append({
            'category': category,
            'n_windows': n,
            'mean_rate': mean_rate,
            'median_rate': median_rate
        })

    # Statistical test
    introner_windows = windows_df[windows_df['category'] == 'introner_containing']['weighted_mean_rate'].values
    non_introner_windows = windows_df[windows_df['category'] == 'non_introner_containing']['weighted_mean_rate'].values

    # Initialize pvalue for plotting
    pvalue = None
    statistic = None

    if len(introner_windows) > 0 and len(non_introner_windows) > 0:
        statistic, pvalue = mannwhitneyu(introner_windows, non_introner_windows, alternative='two-sided')

        print(f"\n🧪 Mann-Whitney U Test:")
        print(f"   Introner-containing windows (n={len(introner_windows)}): median = {np.median(introner_windows):.2e}")
        print(f"   Non-introner-containing windows (n={len(non_introner_windows)}): median = {np.median(non_introner_windows):.2e}")
        print(f"   U-statistic: {statistic:.2e}")
        print(f"   p-value: {pvalue:.2e}")

        if pvalue < 0.05:
            print(f"   ✅ SIGNIFICANT difference (p < 0.05)")
            fold_change = np.median(introner_windows) / np.median(non_introner_windows)
            print(f"   Fold change: {fold_change:.3f}")
        else:
            print(f"   ❌ NO significant difference (p >= 0.05)")
    else:
        print(f"\n⚠️  Cannot perform statistical test:")
        print(f"   Introner-containing windows: {len(introner_windows)}")
        print(f"   Non-introner-containing windows: {len(non_introner_windows)}")

    # Save summary
    summary_df = pd.DataFrame(stats_summary)
    summary_file = f"{args.output_prefix}_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n💾 Saved summary to {summary_file}")

    # Save detailed windows
    detail_file = f"{args.output_prefix}_windows.tsv"
    windows_df.to_csv(detail_file, sep='\t', index=False)
    print(f"💾 Saved window details to {detail_file}")

    # Create visualizations
    print("\n📊 Creating plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Prepare data for plotting (exclude 'excluded' category)
    plot_df = windows_df[windows_df['category'] != 'excluded'].copy()

    # Plot 1: Violin plot
    sns.violinplot(
        data=plot_df,
        x='category',
        y='weighted_mean_rate',
        ax=axes[0, 0],
        inner='box'
    )
    axes[0, 0].set_ylabel('Recombination Rate (r)')
    axes[0, 0].set_xlabel('Window Category')
    axes[0, 0].set_title(f'Recombination Rate Distribution\n({args.window_size}bp windows)')
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_xticklabels(['Non-introner', 'Introner'], rotation=45)

    # Plot 2: Box plot with individual points
    sns.boxplot(
        data=plot_df,
        x='category',
        y='weighted_mean_rate',
        ax=axes[0, 1],
        showfliers=False
    )
    axes[0, 1].set_ylabel('Recombination Rate (r)')
    axes[0, 1].set_xlabel('Window Category')
    title = 'Median Comparison'
    if pvalue is not None:
        title += f'\np={pvalue:.2e}'
    axes[0, 1].set_title(title)
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xticklabels(['Non-introner', 'Introner'], rotation=45)

    # Plot 3: Histogram overlay
    for category, color in [('non_introner_containing', 'blue'),
                            ('introner_containing', 'red')]:
        subset = plot_df[plot_df['category'] == category]['weighted_mean_rate']
        label = 'Introner' if category == 'introner_containing' else 'Non-introner'
        axes[1, 0].hist(np.log10(subset), bins=50, alpha=0.5, label=label, color=color)

    axes[1, 0].set_xlabel('log10(Recombination Rate)')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Rate Distribution (log scale)')
    axes[1, 0].legend()

    # Plot 4: Summary statistics panel
    axes[1, 1].axis('off')

    # Build summary text based on available data
    summary_parts = [
        "STATISTICAL SUMMARY",
        "="*40,
        "",
        "Sample Sizes:",
        f"  Non-introner windows: {len(non_introner_windows):,}",
        f"  Introner windows: {len(introner_windows):,}",
        f"  Excluded windows: {len(windows_df[windows_df['category'] == 'excluded']):,}",
        ""
    ]

    if len(introner_windows) > 0 and len(non_introner_windows) > 0:
        summary_parts.extend([
            "Median Rates:",
            f"  Non-introner: {np.median(non_introner_windows):.2e}",
            f"  Introner: {np.median(introner_windows):.2e}",
            f"  Ratio: {np.median(introner_windows)/np.median(non_introner_windows):.3f}",
            "",
            "Mann-Whitney U Test:",
            f"  U-statistic: {statistic:.2e}",
            f"  p-value: {pvalue:.2e}",
            f"  Result: {'Significant' if pvalue < 0.05 else 'Not significant'}",
            ""
        ])
    else:
        summary_parts.append("Statistical test not performed (insufficient data)")
        summary_parts.append("")

    summary_parts.extend([
        f"Window Size: {args.window_size:,} bp",
        "Mating-type region excluded: Yes"
    ])

    summary_text = "\n".join(summary_parts)

    axes[1, 1].text(0.1, 0.9, summary_text, transform=axes[1, 1].transAxes,
                    fontsize=9, verticalalignment='top', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file)
    print(f"💾 Saved plot to {plot_file}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

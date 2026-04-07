#!/usr/bin/env python3
"""
Compare recombination rates for polymorphic-loss vs polymorphic-gain introners.

Uses genotype matrix to classify introners as:
- Polymorphic-Loss: present in reference (CCMP1545=1), absent in ≥1 other sample (=2)
- Polymorphic-Gain: absent in reference (CCMP1545=2), present in ≥1 other sample (=1)

Performs two independent analyses:
1. Polymorphic-Loss introner windows vs non-introner windows
2. Polymorphic-Gain introner windows vs non-introner windows

Each 10kb window is classified as polymorphic-loss-containing, polymorphic-gain-containing,
or non-introner-containing.
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


def load_genotype_matrix(matrix_file, sample_name='CCMP1545', exclude_samples=['RCC1749', 'RCC3052']):
    """
    Load genotype matrix (long format) and classify introners as polymorphic-loss or polymorphic-gain.

    Matrix format: one row per sample-introner pair
    - ortholog_id: Introner group identifier
    - sample: Sample name
    - contig: Chromosome/contig name
    - start, end: Genomic coordinates
    - presence: 1=present, 2=absent, 3=missing

    Polymorphic-Loss: present in reference (=1), absent in ≥1 other sample (=2)
    Polymorphic-Gain: absent in reference (=2), present in ≥1 other sample (=1)

    Args:
        matrix_file: Path to genotype matrix TSV
        sample_name: Reference sample for classification
        exclude_samples: List of samples to exclude from analysis

    Returns:
        Dict mapping (contig, start, end) → {polymorphic_loss, polymorphic_gain}
    """
    introner_classes = {}

    print(f"   Reading genotype matrix from {matrix_file}...")
    df = pd.read_csv(matrix_file, sep='\t')

    print(f"   Filtering for presence ∈ {{1, 2}} and excluding {exclude_samples}...")

    # Filter for valid presence values (1 or 2) and exclude specified samples
    df = df[df['presence'].isin([1, 2])]
    df = df[~df['sample'].isin(exclude_samples)]

    # Get unique samples in filtered data
    unique_samples = df['sample'].unique()
    print(f"   Found {len(unique_samples)} samples for analysis")

    # Group by ortholog_id to get unique introner regions and their presence across samples
    for orth_id, group in df.groupby('ortholog_id'):
        if len(group) == 0:
            continue

        # Get genomic coordinates (should be same for all rows of same ortholog)
        row = group.iloc[0]
        contig = row['contig']
        start = int(row['start'])
        end = int(row['end'])
        coord_key = (contig, start, end)

        # Get presence status for reference sample in this ortholog group
        ref_rows = group[group['sample'] == sample_name]
        if len(ref_rows) == 0:
            # Reference sample not present for this ortholog
            continue

        ref_presence = int(ref_rows.iloc[0]['presence'])

        # Get presence values for other samples in this ortholog
        other_rows = group[group['sample'] != sample_name]
        other_presences = [int(p) for p in other_rows['presence'].values]

        if len(other_presences) == 0:
            # No other samples with data
            continue

        # Classify
        # Polymorphic-Loss: reference present (1), ≥1 other sample absent (2)
        polymorphic_loss = ref_presence == 1 and 2 in other_presences

        # Polymorphic-Gain: reference absent (2), ≥1 other sample present (1)
        polymorphic_gain = ref_presence == 2 and 1 in other_presences

        introner_classes[coord_key] = {
            'polymorphic_loss': polymorphic_loss,
            'polymorphic_gain': polymorphic_gain
        }

    print(f"   Classified {len(introner_classes)} introner regions")
    return introner_classes


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


def normalize_contig_name(contig_name):
    """
    Normalize contig name to match pyrho output format.

    Converts:
      - "CCMP1545#0#scaffold_10" → "scaffold_10"
      - "RCC114-intronerized#0#10#0" → needs contig extraction
    """
    # Try to extract scaffold number
    if 'scaffold_' in contig_name:
        # Format: "CCMP1545#0#scaffold_10" or similar
        parts = contig_name.split('#')
        for part in parts:
            if part.startswith('scaffold_'):
                return part
    # If no scaffold_ found, return as-is (might already be normalized)
    return contig_name


def classify_window_introner_types(window_start, window_end, chrom, introner_classes):
    """
    Classify window based on polymorphic-loss and polymorphic-gain introners.

    Returns dict with: has_polymorphic_loss, has_polymorphic_gain, n_polymorphic_loss, n_polymorphic_gain
    """
    has_polymorphic_loss = False
    has_polymorphic_gain = False
    n_polymorphic_loss = 0
    n_polymorphic_gain = 0

    for (contig, start, end), classes in introner_classes.items():
        # Normalize contig names for comparison
        norm_contig = normalize_contig_name(contig)
        if norm_contig != chrom:
            continue

        # Check if this introner overlaps the window
        if get_overlap_length(window_start, window_end, start, end) > 0:
            if classes['polymorphic_loss']:
                has_polymorphic_loss = True
                n_polymorphic_loss += 1
            if classes['polymorphic_gain']:
                has_polymorphic_gain = True
                n_polymorphic_gain += 1

    return {
        'has_polymorphic_loss': has_polymorphic_loss,
        'has_polymorphic_gain': has_polymorphic_gain,
        'n_polymorphic_loss': n_polymorphic_loss,
        'n_polymorphic_gain': n_polymorphic_gain
    }


def classify_window(chrom, window_start, window_end, introner_regions,
                   exclude_mating_type=True):
    """
    Classify a 10kb window for original analysis (from BED file).

    Returns:
        'introner_containing', 'non_introner_containing', or 'excluded'
    """
    # Check if overlaps mating-type region (only for CCMP1545)
    if exclude_mating_type and chrom == 'scaffold_2':
        # Mating-type region: scaffold_2:49,808-1,730,591
        if get_overlap_length(window_start, window_end, 49808, 1730591) > 0:
            return 'excluded'

    # Check for introner overlap
    if check_introner_overlap(window_start, window_end, introner_regions.get(chrom, [])):
        return 'introner_containing'

    return 'non_introner_containing'


def perform_mwu_test(test_group, other_group, test_name):
    """
    Perform Mann-Whitney U test and print results.

    Returns: (statistic, pvalue, fold_change) or (None, None, None) if insufficient data
    """
    if len(test_group) == 0 or len(other_group) == 0:
        print(f"\n⚠️  Cannot perform {test_name} test (insufficient data)")
        print(f"   Test group: {len(test_group)}, Other group: {len(other_group)}")
        return None, None, None

    statistic, pvalue = mannwhitneyu(test_group, other_group, alternative='two-sided')
    median_test = np.median(test_group)
    median_other = np.median(other_group)
    fold_change = median_test / median_other if median_other != 0 else np.nan

    print(f"\n🧪 {test_name}:")
    print(f"   Test group (n={len(test_group)}): median = {median_test:.2e}")
    print(f"   Other group (n={len(other_group)}): median = {median_other:.2e}")
    print(f"   U-statistic: {statistic:.2e}")
    print(f"   p-value: {pvalue:.2e}")

    if pvalue < 0.05:
        print(f"   ✅ SIGNIFICANT difference (p < 0.05)")
        print(f"   Fold change: {fold_change:.3f}")
    else:
        print(f"   ❌ NO significant difference (p >= 0.05)")

    return statistic, pvalue, fold_change


def main():
    parser = argparse.ArgumentParser(
        description="Compare recombination rates for polymorphic-loss and polymorphic-gain introners"
    )
    parser.add_argument(
        '--pyrho_dir',
        required=True,
        help='Directory with pyrho output files'
    )
    parser.add_argument(
        '--introner_bed',
        required=True,
        help='BED file with introner regions (for original analysis)'
    )
    parser.add_argument(
        '--introner_matrix',
        required=True,
        help='Genotype matrix TSV file'
    )
    parser.add_argument(
        '--sample_name',
        default='CCMP1545',
        help='Reference sample for polymorphic-loss/gain classification (default: CCMP1545)'
    )
    parser.add_argument(
        '--exclude_samples',
        default='RCC1749,RCC3052',
        help='Comma-separated samples to exclude (default: RCC1749,RCC3052)'
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=5000,
        help='Window size in bp (default: 5000)'
    )
    parser.add_argument(
        '--output_prefix',
        default='polymorphic_loss_gain_introners',
        help='Output file prefix'
    )
    args = parser.parse_args()

    # Parse exclude samples
    exclude_samples = [s.strip() for s in args.exclude_samples.split(',')]

    print("📂 Loading introner annotations from BED...")
    introners = parse_bed_file(args.introner_bed)
    print(f"   Found {sum(len(v) for v in introners.values())} introner regions")

    print("\n📊 Loading and classifying introners from genotype matrix...")
    introner_classes = load_genotype_matrix(
        args.introner_matrix,
        sample_name=args.sample_name,
        exclude_samples=exclude_samples
    )

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

            # Original classification (from BED file)
            classification = classify_window(chrom, window_start, window_end, introners)

            # Classify for polymorphic-loss/gain (from matrix)
            introner_types = classify_window_introner_types(
                window_start, window_end, chrom, introner_classes
            )

            windows.append({
                'scaffold': chrom,
                'start': window_start,
                'end': window_end,
                'weighted_mean_rate': mean_rate,
                'category': classification,
                'has_polymorphic_loss': introner_types['has_polymorphic_loss'],
                'has_polymorphic_gain': introner_types['has_polymorphic_gain'],
                'n_polymorphic_loss': introner_types['n_polymorphic_loss'],
                'n_polymorphic_gain': introner_types['n_polymorphic_gain']
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
            'analysis': category,
            'n_windows': n,
            'mean_rate': mean_rate,
            'median_rate': median_rate
        })

    # Test A: Polymorphic-Loss introners vs non-introner windows
    print("\n" + "="*85)
    print("POLYMORPHIC-LOSS INTRONER ANALYSIS")
    print("="*85)
    print("Introners: present in CCMP1545, absent in ≥1 other Group1 sample")

    non_introner_windows = windows_df[windows_df['category'] == 'non_introner_containing']['weighted_mean_rate'].values
    loss_windows = windows_df[windows_df['has_polymorphic_loss']]['weighted_mean_rate'].values

    stat_loss, pval_loss, fc_loss = perform_mwu_test(
        loss_windows, non_introner_windows,
        "Polymorphic-Loss vs Non-introner"
    )

    if stat_loss is not None:
        stats_summary.append({
            'analysis': 'polymorphic_loss_vs_non',
            'n_windows': f"{len(loss_windows)},{len(non_introner_windows)}",
            'mean_rate': f"{np.mean(loss_windows):.2e}",
            'median_rate': f"{np.median(loss_windows):.2e}"
        })

    # Test B: Polymorphic-Gain introners vs non-introner windows
    print("\n" + "="*85)
    print("POLYMORPHIC-GAIN INTRONER ANALYSIS")
    print("="*85)
    print("Introners: absent in CCMP1545, present in ≥1 other Group1 sample")

    gain_windows = windows_df[windows_df['has_polymorphic_gain']]['weighted_mean_rate'].values

    stat_gain, pval_gain, fc_gain = perform_mwu_test(
        gain_windows, non_introner_windows,
        "Polymorphic-Gain vs Non-introner"
    )

    if stat_gain is not None:
        stats_summary.append({
            'analysis': 'polymorphic_gain_vs_non',
            'n_windows': f"{len(gain_windows)},{len(non_introner_windows)}",
            'mean_rate': f"{np.mean(gain_windows):.2e}",
            'median_rate': f"{np.median(gain_windows):.2e}"
        })

    # Save summary
    summary_df = pd.DataFrame(stats_summary)
    summary_file = f"{args.output_prefix}_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n💾 Saved summary to {summary_file}")

    # Save detailed windows
    detail_file = f"{args.output_prefix}_windows.tsv"
    windows_df.to_csv(detail_file, sep='\t', index=False)
    print(f"💾 Saved window details to {detail_file}")

    # Create 2-panel visualization
    print("\n📊 Creating 2-panel visualization...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Prepare data for plotting (exclude 'excluded' category)
    plot_df = windows_df[windows_df['category'] != 'excluded'].copy()

    # Panel 1: Polymorphic-Loss vs Non-introner
    plot_data_loss = plot_df[['has_polymorphic_loss', 'weighted_mean_rate']].copy()
    plot_data_loss['category'] = plot_data_loss['has_polymorphic_loss'].map({True: 'Loss', False: 'Non-introner'})
    sns.boxplot(
        data=plot_data_loss,
        x='category',
        y='weighted_mean_rate',
        ax=axes[0],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[0].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[0].set_xlabel('')
    title = 'Polymorphic-Loss Introners'
    if pval_loss is not None:
        title += f'\np={pval_loss:.2e}, FC={fc_loss:.3f}'
    axes[0].set_title(title, fontsize=12, fontweight='bold')
    axes[0].set_yscale('log')
    axes[0].set_xticklabels(['Non-introner', 'Loss'])

    # Panel 2: Polymorphic-Gain vs Non-introner
    plot_data_gain = plot_df[['has_polymorphic_gain', 'weighted_mean_rate']].copy()
    plot_data_gain['category'] = plot_data_gain['has_polymorphic_gain'].map({True: 'Gain', False: 'Non-introner'})
    sns.boxplot(
        data=plot_data_gain,
        x='category',
        y='weighted_mean_rate',
        ax=axes[1],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[1].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[1].set_xlabel('')
    title = 'Polymorphic-Gain Introners'
    if pval_gain is not None:
        title += f'\np={pval_gain:.2e}, FC={fc_gain:.3f}'
    axes[1].set_title(title, fontsize=12, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].set_xticklabels(['Non-introner', 'Gain'])

    plt.suptitle(f'Recombination Rates: Polymorphic Loss/Gain Analysis\n{args.window_size}bp windows',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"💾 Saved plot to {plot_file}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

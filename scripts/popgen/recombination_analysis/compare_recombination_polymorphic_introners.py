#!/usr/bin/env python3
"""
Compare recombination rates for polymorphic, gained, and lost introners.

Uses genotype matrix to classify introners as:
- Polymorphic: presence varies across samples
- Gained: absent in reference (CCMP1545=2), present in ≥1 other sample (=1)
- Lost: present in reference (CCMP1545=1), absent in ≥1 other sample (=2)

Performs four independent analyses:
1. Original: Introner-containing vs non-introner-containing (from BED file)
2. Polymorphic introner windows vs others
3. Gained introner windows vs others
4. Lost introner windows vs others

Each 10kb window is classified with all applicable types.
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
    Load genotype matrix (long format) and classify introners as polymorphic/gained/lost.

    Matrix format: one row per sample-introner pair
    - ortholog_id: Introner group identifier
    - sample: Sample name
    - contig: Chromosome/contig name
    - start, end: Genomic coordinates
    - presence: 1=present, 2=absent, 3=missing

    Args:
        matrix_file: Path to genotype matrix TSV
        sample_name: Reference sample for gained/lost classification
        exclude_samples: List of samples to exclude from analysis

    Returns:
        Dict mapping (contig, start, end) → {polymorphic, gained, lost, sample_count}
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
        # Polymorphic: presence values vary (multiple different values)
        all_presences = [ref_presence] + other_presences
        polymorphic = len(set(all_presences)) > 1

        # Gained: reference absent (2), ≥1 other sample present (1)
        gained = ref_presence == 2 and 1 in other_presences

        # Lost: reference present (1), ≥1 other sample absent (2)
        lost = ref_presence == 1 and 2 in other_presences

        introner_classes[coord_key] = {
            'polymorphic': polymorphic,
            'gained': gained,
            'lost': lost,
            'sample_count': len(other_presences)
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


def count_introners_in_window(window_start, window_end, introner_regions):
    """Count number of introners that overlap with window."""
    count = 0
    for introner_start, introner_end in introner_regions:
        if get_overlap_length(window_start, window_end, introner_start, introner_end) > 0:
            count += 1
    return count


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
    Classify window based on introner types (polymorphic, gained, lost).

    Returns dict with: has_polymorphic, has_gained, has_lost, n_polymorphic, n_gained, n_lost
    """
    has_polymorphic = False
    has_gained = False
    has_lost = False
    n_polymorphic = 0
    n_gained = 0
    n_lost = 0

    for (contig, start, end), classes in introner_classes.items():
        # Normalize contig names for comparison
        norm_contig = normalize_contig_name(contig)
        if norm_contig != chrom:
            continue

        # Check if this introner overlaps the window
        if get_overlap_length(window_start, window_end, start, end) > 0:
            if classes['polymorphic']:
                has_polymorphic = True
                n_polymorphic += 1
            if classes['gained']:
                has_gained = True
                n_gained += 1
            if classes['lost']:
                has_lost = True
                n_lost += 1

    return {
        'has_polymorphic': has_polymorphic,
        'has_gained': has_gained,
        'has_lost': has_lost,
        'n_polymorphic': n_polymorphic,
        'n_gained': n_gained,
        'n_lost': n_lost
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
        description="Compare recombination rates for polymorphic, gained, and lost introners"
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
        help='Reference sample for gained/lost classification (default: CCMP1545)'
    )
    parser.add_argument(
        '--exclude_samples',
        default='RCC1749,RCC3052',
        help='Comma-separated samples to exclude (default: RCC1749,RCC3052)'
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=10000,
        help='Window size in bp (default: 10000)'
    )
    parser.add_argument(
        '--output_prefix',
        default='polymorphic_introners',
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

            # Count introners from BED
            n_introners = count_introners_in_window(
                window_start, window_end, introners.get(chrom, [])
            )

            # Classify for polymorphic/gained/lost (from matrix)
            introner_types = classify_window_introner_types(
                window_start, window_end, chrom, introner_classes
            )

            windows.append({
                'scaffold': chrom,
                'start': window_start,
                'end': window_end,
                'weighted_mean_rate': mean_rate,
                'category': classification,
                'n_introners': n_introners,
                'has_polymorphic': introner_types['has_polymorphic'],
                'has_gained': introner_types['has_gained'],
                'has_lost': introner_types['has_lost'],
                'n_polymorphic': introner_types['n_polymorphic'],
                'n_gained': introner_types['n_gained'],
                'n_lost': introner_types['n_lost']
            })

    # Convert to DataFrame
    windows_df = pd.DataFrame(windows)

    print(f"   Created {len(windows_df)} windows total")

    # Summary statistics for original analysis
    print("\n📈 Summary Statistics (Original - from BED file):")
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

    # Original statistical test
    introner_windows_orig = windows_df[windows_df['category'] == 'introner_containing']['weighted_mean_rate'].values
    non_introner_windows_orig = windows_df[windows_df['category'] == 'non_introner_containing']['weighted_mean_rate'].values

    stat_orig, pval_orig, fc_orig = perform_mwu_test(
        introner_windows_orig, non_introner_windows_orig,
        "Original Analysis (Introner-containing vs Non-introner-containing)"
    )

    if stat_orig is not None:
        stats_summary.append({
            'analysis': 'original_introner_vs_non',
            'n_windows': f"{len(introner_windows_orig)},{len(non_introner_windows_orig)}",
            'mean_rate': f"{np.mean(introner_windows_orig):.2e}",
            'median_rate': f"{np.median(introner_windows_orig):.2e}"
        })

    # Test A: Polymorphic introners
    print("\n" + "="*85)
    print("POLYMORPHIC INTRONER ANALYSIS")
    print("="*85)
    polymorphic_windows = windows_df[windows_df['has_polymorphic']]['weighted_mean_rate'].values
    non_polymorphic_windows = windows_df[~windows_df['has_polymorphic']]['weighted_mean_rate'].values
    stat_poly, pval_poly, fc_poly = perform_mwu_test(
        polymorphic_windows, non_polymorphic_windows,
        "Polymorphic Introner Windows vs Non-polymorphic"
    )

    if stat_poly is not None:
        stats_summary.append({
            'analysis': 'polymorphic_vs_non',
            'n_windows': f"{len(polymorphic_windows)},{len(non_polymorphic_windows)}",
            'mean_rate': f"{np.mean(polymorphic_windows):.2e}",
            'median_rate': f"{np.median(polymorphic_windows):.2e}"
        })

    # Test B: Gained introners
    print("\n" + "="*85)
    print("GAINED INTRONER ANALYSIS")
    print("="*85)
    gained_windows = windows_df[windows_df['has_gained']]['weighted_mean_rate'].values
    non_gained_windows = windows_df[~windows_df['has_gained']]['weighted_mean_rate'].values
    stat_gain, pval_gain, fc_gain = perform_mwu_test(
        gained_windows, non_gained_windows,
        "Gained Introner Windows vs Non-gained"
    )

    if stat_gain is not None:
        stats_summary.append({
            'analysis': 'gained_vs_non',
            'n_windows': f"{len(gained_windows)},{len(non_gained_windows)}",
            'mean_rate': f"{np.mean(gained_windows):.2e}",
            'median_rate': f"{np.median(gained_windows):.2e}"
        })

    # Test C: Lost introners
    print("\n" + "="*85)
    print("LOST INTRONER ANALYSIS")
    print("="*85)
    lost_windows = windows_df[windows_df['has_lost']]['weighted_mean_rate'].values
    non_lost_windows = windows_df[~windows_df['has_lost']]['weighted_mean_rate'].values
    stat_lost, pval_lost, fc_lost = perform_mwu_test(
        lost_windows, non_lost_windows,
        "Lost Introner Windows vs Non-lost"
    )

    if stat_lost is not None:
        stats_summary.append({
            'analysis': 'lost_vs_non',
            'n_windows': f"{len(lost_windows)},{len(non_lost_windows)}",
            'mean_rate': f"{np.mean(lost_windows):.2e}",
            'median_rate': f"{np.median(lost_windows):.2e}"
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

    # Create 4-panel visualization
    print("\n📊 Creating 4-panel visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Prepare data for plotting (exclude 'excluded' category)
    plot_df = windows_df[windows_df['category'] != 'excluded'].copy()

    # Panel 1: Original analysis (Introner-containing vs Non)
    plot_data_orig = plot_df[['category', 'weighted_mean_rate']].copy()
    plot_data_orig = plot_data_orig[plot_data_orig['category'] != 'excluded']
    sns.boxplot(
        data=plot_data_orig,
        x='category',
        y='weighted_mean_rate',
        ax=axes[0, 0],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[0, 0].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[0, 0].set_xlabel('')
    title = 'Original: Introner-containing vs Non'
    if pval_orig is not None:
        title += f'\np={pval_orig:.2e}, FC={fc_orig:.3f}'
    axes[0, 0].set_title(title, fontsize=12, fontweight='bold')
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_xticklabels(['Non-introner', 'Introner'])

    # Panel 2: Polymorphic introners
    plot_data_poly = plot_df[['has_polymorphic', 'weighted_mean_rate']].copy()
    plot_data_poly['category'] = plot_data_poly['has_polymorphic'].map({True: 'Polymorphic', False: 'Non-polymorphic'})
    sns.boxplot(
        data=plot_data_poly,
        x='category',
        y='weighted_mean_rate',
        ax=axes[0, 1],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[0, 1].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[0, 1].set_xlabel('')
    title = 'Polymorphic Introners'
    if pval_poly is not None:
        title += f'\np={pval_poly:.2e}, FC={fc_poly:.3f}'
    axes[0, 1].set_title(title, fontsize=12, fontweight='bold')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xticklabels(['Non-polymorphic', 'Polymorphic'])

    # Panel 3: Gained introners
    plot_data_gain = plot_df[['has_gained', 'weighted_mean_rate']].copy()
    plot_data_gain['category'] = plot_data_gain['has_gained'].map({True: 'Gained', False: 'Non-gained'})
    sns.boxplot(
        data=plot_data_gain,
        x='category',
        y='weighted_mean_rate',
        ax=axes[1, 0],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[1, 0].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[1, 0].set_xlabel('')
    title = 'Gained Introners'
    if pval_gain is not None:
        title += f'\np={pval_gain:.2e}, FC={fc_gain:.3f}'
    axes[1, 0].set_title(title, fontsize=12, fontweight='bold')
    axes[1, 0].set_yscale('log')
    axes[1, 0].set_xticklabels(['Non-gained', 'Gained'])

    # Panel 4: Lost introners
    plot_data_lost = plot_df[['has_lost', 'weighted_mean_rate']].copy()
    plot_data_lost['category'] = plot_data_lost['has_lost'].map({True: 'Lost', False: 'Non-lost'})
    sns.boxplot(
        data=plot_data_lost,
        x='category',
        y='weighted_mean_rate',
        ax=axes[1, 1],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[1, 1].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[1, 1].set_xlabel('')
    title = 'Lost Introners'
    if pval_lost is not None:
        title += f'\np={pval_lost:.2e}, FC={fc_lost:.3f}'
    axes[1, 1].set_title(title, fontsize=12, fontweight='bold')
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_xticklabels(['Non-lost', 'Lost'])

    plt.suptitle(f'Recombination Rates: Polymorphic Introner Analysis\n{args.window_size}bp windows',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"💾 Saved plot to {plot_file}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

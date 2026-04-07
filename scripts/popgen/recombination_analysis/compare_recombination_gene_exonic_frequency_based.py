#!/usr/bin/env python3
"""
Gene-centric frequency-based introner recombination analysis.

Restricts analysis to exonic regions only (gene background).
Classifies polymorphic introners by allele frequency (count of Group1 samples with presence=1):
- Recent Gains: 1-3 samples (9-27% frequency) - brand new acquisitions
- Recent Losses: 8-10 samples (73-91% frequency) - well-established introners being purged

Compares recent-gain and recent-loss exonic windows vs non-introner-containing windows.

Algorithm:
1. Parse GTF to get exon positions per gene
2. Extract and classify introners from genotype matrix by frequency
3. Cluster introners per gene (merge <5kb apart)
4. Generate windows per gene:
   - For recent gain clusters: center 5kb window on cluster midpoint
   - For recent loss clusters: center 5kb window on cluster midpoint
   - For non-introner exonic regions: tile with 5kb windows
5. Calculate weighted mean recombination rates from pyrho
6. Perform two Mann-Whitney U tests: recent_gain vs non, recent_loss vs non
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


def load_gtf(gtf_file):
    """
    Parse GTF file and extract exon positions per gene.

    Returns:
        dict: {gene_id: {'exons': [(start, end), ...], 'contig': contig}}
    """
    genes = defaultdict(lambda: {'exons': [], 'contig': None})

    print(f"   Parsing GTF file: {gtf_file}...")
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'exon':
                continue

            contig = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])  # End is already exclusive in 0-based

            # Extract gene_id from attributes
            attributes = fields[8]
            gene_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break

            if gene_id:
                genes[gene_id]['exons'].append((start, end))
                genes[gene_id]['contig'] = contig

    # Merge overlapping exons for each gene
    for gene_id in genes:
        exons = sorted(genes[gene_id]['exons'])
        merged = []
        for start, end in exons:
            if merged and start <= merged[-1][1]:
                # Overlapping or adjacent exons - merge
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        genes[gene_id]['exons'] = merged

    print(f"   Found {len(genes)} genes with exon annotations")
    return genes


def extract_introners_by_frequency(matrix_file, exclude_samples=['RCC1749', 'RCC3052']):
    """
    Extract introner positions from genotype matrix and classify by allele frequency in Group1.

    Classifies polymorphic introners by count of Group1 samples with presence=1:
    - Recent Gains: 1-3 samples (9-27% frequency)
    - Recent Losses: 8-10 samples (73-91% frequency)

    Args:
        matrix_file: Path to genotype matrix TSV
        exclude_samples: List of samples to exclude (Group2: RCC1749, RCC3052)

    Returns:
        dict: {(contig, start, end): {'type': 'recent_gain' or 'recent_loss', 'frequency': count}}
    """
    introners = {}

    # Define Group1 samples
    GROUP1_SAMPLES = {
        'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482',
        'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
    }

    print(f"   Reading genotype matrix: {matrix_file}...")
    df = pd.read_csv(matrix_file, sep='\t')

    print(f"   Filtering for presence ∈ {{1, 2}} and excluding {exclude_samples}...")

    # Filter for valid presence values and exclude Group2 samples
    df = df[df['presence'].isin([1, 2])]
    df = df[~df['sample'].isin(exclude_samples)]

    print(f"   Analyzing {len(GROUP1_SAMPLES)} Group1 samples")

    # Classify each ortholog group by frequency
    for orth_id, group in df.groupby('ortholog_id'):
        if len(group) == 0:
            continue

        # Filter to Group1 samples only
        group_g1 = group[group['sample'].isin(GROUP1_SAMPLES)]

        if len(group_g1) != len(GROUP1_SAMPLES):
            # Not all Group1 samples have data for this ortholog
            continue

        presences = group_g1['presence'].values

        # Check for polymorphism (need both 1 and 2)
        has_present = 1 in presences
        has_absent = 2 in presences

        if not (has_present and has_absent):
            continue

        # Count samples with presence=1
        count_present = int(np.sum(presences == 1))

        # Classify by frequency
        if count_present in [1, 2, 3]:
            introner_type = 'recent_gain'
        elif count_present in [8, 9, 10]:
            introner_type = 'recent_loss'
        else:
            # Skip intermediate frequencies
            continue

        # Get genomic coordinates from first row (coordinates are same for all samples)
        row = group_g1.iloc[0]
        contig = row['contig']
        start = int(row['start'])
        end = int(row['end'])
        coord_key = (contig, start, end)

        introners[coord_key] = {
            'type': introner_type,
            'frequency': count_present
        }

    print(f"   Classified {len(introners)} introner regions (recent_gain and recent_loss)")
    return introners


def normalize_contig_name(contig_name):
    """
    Normalize contig name to match pyrho output format.
    Converts: "CCMP1545#0#scaffold_10" → "scaffold_10"
    """
    if 'scaffold_' in contig_name:
        parts = contig_name.split('#')
        for part in parts:
            if part.startswith('scaffold_'):
                return part
    return contig_name


def get_overlap_length(start1, end1, start2, end2):
    """Calculate overlap length between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def cluster_introners(introners, gene_id, gene_data, merge_distance=5000):
    """
    Cluster introners within a gene that are <merge_distance apart.
    Tracks whether each cluster is recent_gain or recent_loss.

    Returns:
        dict: {recent_gain: [(cluster_start, cluster_end, midpoint), ...],
               recent_loss: [(cluster_start, cluster_end, midpoint), ...]}
    """
    gene_contig = gene_data['contig']
    gene_exons = gene_data['exons']

    # Separate recent_gain and recent_loss introners for this gene
    recent_gain_introners = []
    recent_loss_introners = []

    for (contig, start, end), info in introners.items():
        if normalize_contig_name(contig) != normalize_contig_name(gene_contig):
            continue

        # Check if introner overlaps with gene exons
        overlaps_exon = False
        for exon_start, exon_end in gene_exons:
            if get_overlap_length(start, end, exon_start, exon_end) > 0:
                overlaps_exon = True
                break

        if overlaps_exon:
            if info['type'] == 'recent_gain':
                recent_gain_introners.append((start, end))
            elif info['type'] == 'recent_loss':
                recent_loss_introners.append((start, end))

    # Cluster separately for recent_gain and recent_loss
    recent_gain_clusters = _cluster_list(recent_gain_introners, merge_distance)
    recent_loss_clusters = _cluster_list(recent_loss_introners, merge_distance)

    return {'recent_gain': recent_gain_clusters, 'recent_loss': recent_loss_clusters}


def _cluster_list(introner_list, merge_distance):
    """Helper function to cluster a list of introners."""
    if len(introner_list) == 0:
        return []

    introner_list = sorted(introner_list)
    clusters = []
    current_cluster = [introner_list[0]]

    for i in range(1, len(introner_list)):
        start, end = introner_list[i]
        prev_start, prev_end = current_cluster[-1]

        # Check distance between end of previous and start of current
        if start - prev_end <= merge_distance:
            # Within merge distance - add to cluster
            current_cluster.append((start, end))
        else:
            # Too far apart - finalize current cluster
            cluster_start = current_cluster[0][0]
            cluster_end = current_cluster[-1][1]
            midpoint = (cluster_start + cluster_end) // 2
            clusters.append((cluster_start, cluster_end, midpoint))

            # Start new cluster
            current_cluster = [(start, end)]

    # Finalize last cluster
    if current_cluster:
        cluster_start = current_cluster[0][0]
        cluster_end = current_cluster[-1][1]
        midpoint = (cluster_start + cluster_end) // 2
        clusters.append((cluster_start, cluster_end, midpoint))

    return clusters


def is_window_in_exons(window_start, window_end, exons):
    """Check if window has any overlap with exonic regions."""
    for exon_start, exon_end in exons:
        if get_overlap_length(window_start, window_end, exon_start, exon_end) > 0:
            return True
    return False


def generate_gene_windows(gene_id, gene_data, introner_clusters, window_size=5000):
    """
    Generate 5kb windows for a gene, classifying by recent_gain/recent_loss/non-introner.

    For recent_gain/recent_loss clusters: center window on cluster midpoint
    For non-introner exonic regions: tile with 5kb windows

    Returns:
        list: [{'start': int, 'end': int, 'type': str, ...}, ...]
    """
    windows = []
    exons = gene_data['exons']

    # Windows centered on recent_gain clusters
    for cluster_start, cluster_end, midpoint in introner_clusters.get('recent_gain', []):
        window_start = max(0, midpoint - window_size // 2)
        window_end = window_start + window_size

        if is_window_in_exons(window_start, window_end, exons):
            windows.append({
                'gene_id': gene_id,
                'start': window_start,
                'end': window_end,
                'type': 'recent_gain_containing',
                'cluster_midpoint': midpoint
            })

    # Windows centered on recent_loss clusters
    for cluster_start, cluster_end, midpoint in introner_clusters.get('recent_loss', []):
        window_start = max(0, midpoint - window_size // 2)
        window_end = window_start + window_size

        if is_window_in_exons(window_start, window_end, exons):
            windows.append({
                'gene_id': gene_id,
                'start': window_start,
                'end': window_end,
                'type': 'recent_loss_containing',
                'cluster_midpoint': midpoint
            })

    # Get non-introner exonic regions
    remaining_regions = list(exons)

    # Remove regions covered by recent_gain clusters
    for cluster_start, cluster_end, _ in introner_clusters.get('recent_gain', []):
        new_remaining = []
        for region_start, region_end in remaining_regions:
            if cluster_end <= region_start or cluster_start >= region_end:
                new_remaining.append((region_start, region_end))
            else:
                if region_start < cluster_start:
                    new_remaining.append((region_start, cluster_start))
                if region_end > cluster_end:
                    new_remaining.append((cluster_end, region_end))
        remaining_regions = new_remaining

    # Remove regions covered by recent_loss clusters
    for cluster_start, cluster_end, _ in introner_clusters.get('recent_loss', []):
        new_remaining = []
        for region_start, region_end in remaining_regions:
            if cluster_end <= region_start or cluster_start >= region_end:
                new_remaining.append((region_start, region_end))
            else:
                if region_start < cluster_start:
                    new_remaining.append((region_start, cluster_start))
                if region_end > cluster_end:
                    new_remaining.append((cluster_end, region_end))
        remaining_regions = new_remaining

    # Tile non-introner regions with 5kb windows
    for region_start, region_end in remaining_regions:
        for window_start in range(region_start, region_end, window_size):
            window_end = min(window_start + window_size, region_end)

            # Skip windows that are too small (< 50% of window size)
            if window_end - window_start < window_size * 0.5:
                continue

            if is_window_in_exons(window_start, window_end, exons):
                windows.append({
                    'gene_id': gene_id,
                    'start': window_start,
                    'end': window_end,
                    'type': 'non_introner_containing',
                    'cluster_midpoint': None
                })

    return windows


def parse_pyrho_file(filepath):
    """Parse a single pyrho output file."""
    df = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        names=['start', 'end', 'rate']
    )
    return df


def calculate_weighted_mean_rate(window_start, window_end, pyrho_windows):
    """Calculate weighted mean recombination rate for a window."""
    mask = (pyrho_windows['end'] > window_start) & (pyrho_windows['start'] < window_end)
    overlapping = pyrho_windows[mask]

    if len(overlapping) == 0:
        return np.nan

    overlap_starts = np.maximum(overlapping['start'].values, window_start)
    overlap_ends = np.minimum(overlapping['end'].values, window_end)
    overlap_lengths = overlap_ends - overlap_starts

    total_weighted_rate = np.sum(overlapping['rate'].values * overlap_lengths)
    total_length = np.sum(overlap_lengths)

    if total_length == 0:
        return np.nan

    return total_weighted_rate / total_length


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
        description="Gene-centric loss/gain introner recombination analysis"
    )
    parser.add_argument(
        '--pyrho_dir',
        required=True,
        help='Directory with pyrho output files'
    )
    parser.add_argument(
        '--gtf_file',
        required=True,
        help='GTF file with gene annotations'
    )
    parser.add_argument(
        '--introner_matrix',
        required=True,
        help='Genotype matrix TSV file'
    )
    parser.add_argument(
        '--sample_name',
        default='CCMP1545',
        help='Reference sample for loss/gain classification (default: CCMP1545)'
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
        '--merge_distance',
        type=int,
        default=5000,
        help='Distance to merge nearby introners (default: 5000)'
    )
    parser.add_argument(
        '--output_prefix',
        default='gene_exonic_loss_gain_introners',
        help='Output file prefix'
    )
    args = parser.parse_args()

    # Parse exclude samples
    exclude_samples = [s.strip() for s in args.exclude_samples.split(',')]

    print("📂 Loading gene annotations from GTF...")
    genes = load_gtf(args.gtf_file)

    print(f"\n📂 Extracting frequency-based introners from genotype matrix...")
    introners = extract_introners_by_frequency(
        args.introner_matrix,
        exclude_samples=exclude_samples
    )

    print("\n📊 Clustering introners per gene...")
    gene_introner_clusters = {}
    n_genes_with_gain = 0
    n_genes_with_loss = 0
    for gene_id, gene_data in genes.items():
        clusters = cluster_introners(introners, gene_id, gene_data, args.merge_distance)
        if len(clusters['recent_gain']) > 0:
            n_genes_with_gain += 1
        if len(clusters['recent_loss']) > 0:
            n_genes_with_loss += 1
        gene_introner_clusters[gene_id] = clusters

    print(f"   {n_genes_with_gain} genes contain recent_gain introners")
    print(f"   {n_genes_with_loss} genes contain recent_loss introners")

    print("\n🔲 Generating gene-based windows...")
    all_windows = []
    for gene_id, gene_data in genes.items():
        clusters = gene_introner_clusters[gene_id]
        windows = generate_gene_windows(gene_id, gene_data, clusters, args.window_size)
        all_windows.extend(windows)

    print(f"   Generated {len(all_windows)} exonic windows total")

    print("\n📊 Loading pyrho data...")
    pyrho_data = {}
    pyrho_files = glob.glob(os.path.join(args.pyrho_dir, '*.pyrho.out'))

    for pyrho_file in sorted(pyrho_files):
        basename = os.path.basename(pyrho_file)
        chrom = basename.replace('.pyrho.out', '')
        pyrho_data[chrom] = parse_pyrho_file(pyrho_file)
        print(f"   Loaded {chrom}: {len(pyrho_data[chrom])} pyrho windows")

    print("\n⚡ Calculating recombination rates...")
    excluded_count = 0
    for window in all_windows:
        gene_id = window['gene_id']
        contig = genes[gene_id]['contig']
        norm_contig = normalize_contig_name(contig)

        # Check if in mating-type region
        if norm_contig == 'scaffold_2':
            if 49808 <= window['start'] < 1730591:
                window['rate'] = np.nan
                excluded_count += 1
                continue

        # Get recombination rate from pyrho
        if norm_contig not in pyrho_data:
            window['rate'] = np.nan
            continue

        window['rate'] = calculate_weighted_mean_rate(
            window['start'], window['end'], pyrho_data[norm_contig]
        )

    if excluded_count > 0:
        print(f"   Excluded {excluded_count} windows in mating-type region (scaffold_2:49,808-1,730,591)")

    # Convert to DataFrame and remove NaN rates
    windows_df = pd.DataFrame(all_windows)
    windows_df = windows_df[~windows_df['rate'].isna()].copy()
    print(f"   Calculated rates for {len(windows_df)} windows")

    # Summary statistics
    print("\n📈 Summary Statistics:")
    print(f"\n{'Window Type':<30} {'N Windows':<12} {'Mean Rate':<20} {'Median Rate':<20}")
    print("-" * 85)

    stats_summary = []
    for wtype in ['non_introner_containing', 'recent_loss_containing', 'recent_gain_containing']:
        subset = windows_df[windows_df['type'] == wtype]
        n = len(subset)
        mean_rate = np.mean(subset['rate']) if n > 0 else np.nan
        median_rate = np.median(subset['rate']) if n > 0 else np.nan

        print(f"{wtype:<30} {n:<12} {mean_rate:<20.2e} {median_rate:<20.2e}")

        stats_summary.append({
            'window_type': wtype,
            'n_windows': n,
            'mean_rate': mean_rate,
            'median_rate': median_rate
        })

    # Test A: Recent Loss vs Non-introner
    print("\n" + "="*85)
    print("RECENT-LOSS INTRONER ANALYSIS")
    print("="*85)
    print("Introners: present in 8-10 Group1 samples (73-91% frequency)")

    non_introner_windows = windows_df[windows_df['type'] == 'non_introner_containing']['rate'].values
    recent_loss_windows = windows_df[windows_df['type'] == 'recent_loss_containing']['rate'].values

    stat_loss, pval_loss, fc_loss = perform_mwu_test(
        recent_loss_windows, non_introner_windows,
        "Recent-Loss vs Non-introner"
    )

    if stat_loss is not None:
        stats_summary.append({
            'window_type': 'recent_loss_vs_non',
            'n_windows': f"{len(recent_loss_windows)},{len(non_introner_windows)}",
            'mean_rate': f"{np.mean(recent_loss_windows):.2e}",
            'median_rate': f"{np.median(recent_loss_windows):.2e}"
        })

    # Test B: Recent Gain vs Non-introner
    print("\n" + "="*85)
    print("RECENT-GAIN INTRONER ANALYSIS")
    print("="*85)
    print("Introners: present in 1-3 Group1 samples (9-27% frequency)")

    recent_gain_windows = windows_df[windows_df['type'] == 'recent_gain_containing']['rate'].values

    stat_gain, pval_gain, fc_gain = perform_mwu_test(
        recent_gain_windows, non_introner_windows,
        "Recent-Gain vs Non-introner"
    )

    if stat_gain is not None:
        stats_summary.append({
            'window_type': 'recent_gain_vs_non',
            'n_windows': f"{len(recent_gain_windows)},{len(non_introner_windows)}",
            'mean_rate': f"{np.mean(recent_gain_windows):.2e}",
            'median_rate': f"{np.median(recent_gain_windows):.2e}"
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

    # Panel 1: Recent Loss vs Non-introner
    plot_data_loss = windows_df[windows_df['type'].isin(['recent_loss_containing', 'non_introner_containing'])].copy()
    plot_data_loss['category'] = plot_data_loss['type'].map({
        'recent_loss_containing': 'Recent Loss',
        'non_introner_containing': 'Non-introner'
    })
    sns.boxplot(
        data=plot_data_loss,
        x='category',
        y='rate',
        ax=axes[0],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[0].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[0].set_xlabel('')
    title = 'Recent-Loss Introners (8-10 samples)'
    if pval_loss is not None:
        title += f'\np={pval_loss:.2e}, FC={fc_loss:.3f}'
    axes[0].set_title(title, fontsize=12, fontweight='bold')
    axes[0].set_yscale('log')
    axes[0].set_xticklabels(['Non-introner', 'Recent Loss'])

    # Panel 2: Recent Gain vs Non-introner
    plot_data_gain = windows_df[windows_df['type'].isin(['recent_gain_containing', 'non_introner_containing'])].copy()
    plot_data_gain['category'] = plot_data_gain['type'].map({
        'recent_gain_containing': 'Recent Gain',
        'non_introner_containing': 'Non-introner'
    })
    sns.boxplot(
        data=plot_data_gain,
        x='category',
        y='rate',
        ax=axes[1],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[1].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[1].set_xlabel('')
    title = 'Recent-Gain Introners (1-3 samples)'
    if pval_gain is not None:
        title += f'\np={pval_gain:.2e}, FC={fc_gain:.3f}'
    axes[1].set_title(title, fontsize=12, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].set_xticklabels(['Non-introner', 'Recent Gain'])

    plt.suptitle(f'Gene-Exonic Recombination: Frequency-Based Allele Classification\n{args.window_size}bp windows',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"💾 Saved plot to {plot_file}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Gene-centric loss/gain introner recombination analysis.

Restricts analysis to exonic regions only (gene background).
Classifies polymorphic introners as:
- Polymorphic-Loss: present in reference (CCMP1545=1), absent in ≥1 other Group1 sample (=2)
- Polymorphic-Gain: absent in reference (CCMP1545=2), present in ≥1 other Group1 sample (=1)

Compares loss-containing and gain-containing exonic windows vs non-introner-containing windows.

Algorithm:
1. Parse GTF to get exon positions per gene
2. Extract and classify introners from genotype matrix (loss or gain)
3. Cluster introners per gene (merge <10kb apart)
4. Generate windows per gene:
   - For loss clusters: center 10kb window on cluster midpoint
   - For gain clusters: center 10kb window on cluster midpoint
   - For non-introner exonic regions: tile with 10kb windows
5. Calculate weighted mean recombination rates from pyrho
6. Perform two Mann-Whitney U tests: loss vs non, gain vs non
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


def extract_introners_loss_gain(matrix_file, sample_name='CCMP1545', exclude_samples=['RCC1749', 'RCC3052']):
    """
    Extract introner positions from genotype matrix and classify as loss or gain.

    Polymorphic-Loss: present in reference (=1), absent in ≥1 other Group1 sample (=2)
    Polymorphic-Gain: absent in reference (=2), present in ≥1 other Group1 sample (=1)

    Args:
        matrix_file: Path to genotype matrix TSV
        sample_name: Reference sample (CCMP1545)
        exclude_samples: List of samples to exclude (Group2)

    Returns:
        dict: {(contig, start, end): {'type': 'loss' or 'gain'}}
    """
    introners = {}

    print(f"   Reading genotype matrix: {matrix_file}...")
    df = pd.read_csv(matrix_file, sep='\t')

    print(f"   Filtering for presence ∈ {{1, 2}} and excluding {exclude_samples}...")

    # Filter for valid presence values and exclude Group2 samples
    df = df[df['presence'].isin([1, 2])]
    df = df[~df['sample'].isin(exclude_samples)]

    # Get Group1 samples
    group1_samples = set(df['sample'].unique())
    print(f"   Analyzing Group1 samples: {len(group1_samples)} samples")

    # Classify each ortholog group as loss or gain
    for orth_id, group in df.groupby('ortholog_id'):
        if len(group) == 0:
            continue

        # Get reference sample presence
        ref_rows = group[group['sample'] == sample_name]
        if len(ref_rows) == 0:
            continue

        ref_presence = int(ref_rows.iloc[0]['presence'])

        # Get other samples' presence
        other_rows = group[group['sample'] != sample_name]
        other_presences = [int(p) for p in other_rows['presence'].values]

        if len(other_presences) == 0:
            continue

        # Classify as loss or gain
        polymorphic_loss = ref_presence == 1 and 2 in other_presences
        polymorphic_gain = ref_presence == 2 and 1 in other_presences

        if polymorphic_loss or polymorphic_gain:
            # Get genomic coordinates from reference sample row
            row = ref_rows.iloc[0]
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            coord_key = (contig, start, end)

            introner_type = 'loss' if polymorphic_loss else 'gain'
            introners[coord_key] = {'type': introner_type}

    print(f"   Classified {len(introners)} introner regions (loss and gain)")
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


def cluster_introners(introners, gene_id, gene_data, merge_distance=10000):
    """
    Cluster introners within a gene that are <merge_distance apart.
    Tracks whether each cluster is loss or gain.

    Returns:
        dict: {loss: [(cluster_start, cluster_end, midpoint), ...],
               gain: [(cluster_start, cluster_end, midpoint), ...]}
    """
    gene_contig = gene_data['contig']
    gene_exons = gene_data['exons']

    # Separate loss and gain introners for this gene
    loss_introners = []
    gain_introners = []

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
            if info['type'] == 'loss':
                loss_introners.append((start, end))
            elif info['type'] == 'gain':
                gain_introners.append((start, end))

    # Cluster separately for loss and gain
    loss_clusters = _cluster_list(loss_introners, merge_distance)
    gain_clusters = _cluster_list(gain_introners, merge_distance)

    return {'loss': loss_clusters, 'gain': gain_clusters}


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
    Generate 5kb windows for a gene, classifying by loss/gain/non-introner.

    For loss/gain clusters: center window on cluster midpoint
    For non-introner exonic regions: tile with 5kb windows

    Returns:
        list: [{'start': int, 'end': int, 'type': str, ...}, ...]
    """
    windows = []
    exons = gene_data['exons']

    # Windows centered on loss clusters
    for cluster_start, cluster_end, midpoint in introner_clusters.get('loss', []):
        window_start = max(0, midpoint - window_size // 2)
        window_end = window_start + window_size

        if is_window_in_exons(window_start, window_end, exons):
            windows.append({
                'gene_id': gene_id,
                'start': window_start,
                'end': window_end,
                'type': 'loss_containing',
                'cluster_midpoint': midpoint
            })

    # Windows centered on gain clusters
    for cluster_start, cluster_end, midpoint in introner_clusters.get('gain', []):
        window_start = max(0, midpoint - window_size // 2)
        window_end = window_start + window_size

        if is_window_in_exons(window_start, window_end, exons):
            windows.append({
                'gene_id': gene_id,
                'start': window_start,
                'end': window_end,
                'type': 'gain_containing',
                'cluster_midpoint': midpoint
            })

    # Get non-introner exonic regions
    remaining_regions = list(exons)

    # Remove regions covered by loss clusters
    for cluster_start, cluster_end, _ in introner_clusters.get('loss', []):
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

    # Remove regions covered by gain clusters
    for cluster_start, cluster_end, _ in introner_clusters.get('gain', []):
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

    print(f"\n📂 Extracting loss/gain introners from genotype matrix...")
    introners = extract_introners_loss_gain(
        args.introner_matrix,
        sample_name=args.sample_name,
        exclude_samples=exclude_samples
    )

    print("\n📊 Clustering introners per gene...")
    gene_introner_clusters = {}
    n_genes_with_loss = 0
    n_genes_with_gain = 0
    for gene_id, gene_data in genes.items():
        clusters = cluster_introners(introners, gene_id, gene_data, args.merge_distance)
        if len(clusters['loss']) > 0:
            n_genes_with_loss += 1
        if len(clusters['gain']) > 0:
            n_genes_with_gain += 1
        gene_introner_clusters[gene_id] = clusters

    print(f"   {n_genes_with_loss} genes contain loss introners")
    print(f"   {n_genes_with_gain} genes contain gain introners")

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
    for wtype in ['non_introner_containing', 'loss_containing', 'gain_containing']:
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

    # Test A: Loss vs Non-introner
    print("\n" + "="*85)
    print("POLYMORPHIC-LOSS INTRONER ANALYSIS")
    print("="*85)
    print("Introners: present in CCMP1545, absent in ≥1 other Group1 sample")

    non_introner_windows = windows_df[windows_df['type'] == 'non_introner_containing']['rate'].values
    loss_windows = windows_df[windows_df['type'] == 'loss_containing']['rate'].values

    stat_loss, pval_loss, fc_loss = perform_mwu_test(
        loss_windows, non_introner_windows,
        "Polymorphic-Loss vs Non-introner"
    )

    if stat_loss is not None:
        stats_summary.append({
            'window_type': 'loss_vs_non',
            'n_windows': f"{len(loss_windows)},{len(non_introner_windows)}",
            'mean_rate': f"{np.mean(loss_windows):.2e}",
            'median_rate': f"{np.median(loss_windows):.2e}"
        })

    # Test B: Gain vs Non-introner
    print("\n" + "="*85)
    print("POLYMORPHIC-GAIN INTRONER ANALYSIS")
    print("="*85)
    print("Introners: absent in CCMP1545, present in ≥1 other Group1 sample")

    gain_windows = windows_df[windows_df['type'] == 'gain_containing']['rate'].values

    stat_gain, pval_gain, fc_gain = perform_mwu_test(
        gain_windows, non_introner_windows,
        "Polymorphic-Gain vs Non-introner"
    )

    if stat_gain is not None:
        stats_summary.append({
            'window_type': 'gain_vs_non',
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

    # Panel 1: Loss vs Non-introner
    plot_data_loss = windows_df[windows_df['type'].isin(['loss_containing', 'non_introner_containing'])].copy()
    plot_data_loss['category'] = plot_data_loss['type'].map({
        'loss_containing': 'Loss',
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
    title = 'Polymorphic-Loss Introners'
    if pval_loss is not None:
        title += f'\np={pval_loss:.2e}, FC={fc_loss:.3f}'
    axes[0].set_title(title, fontsize=12, fontweight='bold')
    axes[0].set_yscale('log')
    axes[0].set_xticklabels(['Non-introner', 'Loss'])

    # Panel 2: Gain vs Non-introner
    plot_data_gain = windows_df[windows_df['type'].isin(['gain_containing', 'non_introner_containing'])].copy()
    plot_data_gain['category'] = plot_data_gain['type'].map({
        'gain_containing': 'Gain',
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
    title = 'Polymorphic-Gain Introners'
    if pval_gain is not None:
        title += f'\np={pval_gain:.2e}, FC={fc_gain:.3f}'
    axes[1].set_title(title, fontsize=12, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].set_xticklabels(['Non-introner', 'Gain'])

    plt.suptitle(f'Gene-Exonic Recombination: Polymorphic Loss/Gain Analysis\n{args.window_size}bp windows',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plot_file = f"{args.output_prefix}.pdf"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"💾 Saved plot to {plot_file}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Group 2 Gene-centric, introner-centered recombination analysis.

Analyzes recombination rates in introner-containing regions for Group 2 samples (RCC1749, RCC3052).
Restricts analysis to exonic regions only (gene background).
Uses introner-centered windows via clustering of nearby introners (<10kb apart).
Compares introner-containing exonic windows vs non-introner-containing exonic windows.

Algorithm:
1. Parse GTF to get exon positions per gene
2. Extract introners from genotype matrix where RCC1749 or RCC3052 has presence == 1
3. Exclude configured contigs (including the RCC1749 mating-type contig)
4. Cluster introners per gene (merge <10kb apart)
5. Generate windows per gene:
   - For introner clusters: center 10kb window on cluster midpoint
   - For non-introner exonic regions: tile with 10kb windows
6. Calculate weighted mean recombination rates from pyrho
7. Perform Mann-Whitney U test
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
from pathlib import Path


def save_plot_bundle(output_prefix):
    """Save the current figure as PDF and PNG using the shared prefix."""
    prefix = Path(output_prefix)
    plt.savefig(prefix.with_suffix(".pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches='tight')


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


def extract_group2_introners(matrix_file):
    """
    Extract introner positions from genotype matrix for Group 2 samples.

    Extracts introners where RCC1749 or RCC3052 has presence == 1.

    Args:
        matrix_file: Path to genotype matrix TSV

    Returns:
        dict: {(contig, start, end): {introner_info}} indexed by genomic coordinates
    """
    introners = {}
    group2_samples = ['RCC1749', 'RCC3052']

    print(f"   Reading genotype matrix: {matrix_file}...")
    df = pd.read_csv(matrix_file, sep='\t')

    # Filter for presence == 1 (introner present) AND Group2 samples
    df_filtered = df[(df['presence'] == 1) & (df['sample'].isin(group2_samples))]
    print(f"   Filtered to {len(df_filtered)} present introner records for Group2 samples {group2_samples}")

    # Get unique introner locations (by ortholog_id to avoid duplicates from both samples)
    unique_orthologs = {}
    for _, row in df_filtered.iterrows():
        orth_id = row['ortholog_id']
        if orth_id not in unique_orthologs:
            unique_orthologs[orth_id] = row

    # Get unique introner locations
    for orth_id, row in unique_orthologs.items():
        contig = row['contig']
        start = int(row['start'])
        end = int(row['end'])
        coord_key = (contig, start, end)

        if coord_key not in introners:
            introners[coord_key] = {
                'contig': contig,
                'start': start,
                'end': end,
                'ortholog_id': orth_id
            }

    print(f"   Found {len(introners)} unique introner locations in Group2 (RCC1749/RCC3052)")
    return introners


def normalize_contig_name(contig_name):
    """
    Normalize contig name to match pyrho output format.
    Converts: "RCC1749#0#contig_10" → "contig_10"
    """
    if '#' in contig_name:
        parts = contig_name.split('#')
        # Return the last part which should be the contig name
        return parts[-1]
    return contig_name


def get_overlap_length(start1, end1, start2, end2):
    """Calculate overlap length between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def cluster_introners(introners, gene_id, gene_data, merge_distance=10000):
    """
    Cluster introners within a gene that are <merge_distance apart.

    Returns:
        list: [(cluster_start, cluster_end, midpoint), ...]
    """
    gene_contig = gene_data['contig']
    gene_exons = gene_data['exons']

    # Get introners that fall within this gene's exons
    gene_introners = []
    for (contig, start, end), introner_info in introners.items():
        if normalize_contig_name(contig) != normalize_contig_name(gene_contig):
            continue

        # Check if introner overlaps with gene exons
        for exon_start, exon_end in gene_exons:
            if get_overlap_length(start, end, exon_start, exon_end) > 0:
                gene_introners.append((start, end))
                break

    if len(gene_introners) == 0:
        return []

    # Sort introners
    gene_introners = sorted(gene_introners)

    # Cluster introners that are within merge_distance
    clusters = []
    current_cluster = [gene_introners[0]]

    for i in range(1, len(gene_introners)):
        start, end = gene_introners[i]
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


def generate_gene_windows(gene_id, gene_data, introner_clusters, window_size=10000):
    """
    Generate 10kb windows for a gene.

    For introner clusters: center window on cluster midpoint
    For non-introner exonic regions: tile with 10kb windows

    Returns:
        list: [{'start': int, 'end': int, 'type': 'introner' or 'non_introner'}, ...]
    """
    windows = []
    exons = gene_data['exons']

    # Windows centered on introner clusters
    for cluster_start, cluster_end, midpoint in introner_clusters:
        window_start = max(0, midpoint - window_size // 2)
        window_end = window_start + window_size

        # Make sure window is within exonic regions (some overlap okay)
        if is_window_in_exons(window_start, window_end, exons):
            windows.append({
                'gene_id': gene_id,
                'start': window_start,
                'end': window_end,
                'type': 'introner_containing',
                'cluster_midpoint': midpoint
            })

    # Get non-introner exonic regions using interval subtraction
    # Start with all exonic regions, remove parts covered by introner clusters
    remaining_regions = list(exons)  # [(start, end), ...]

    # For each introner cluster, subtract it from remaining regions
    for cluster_start, cluster_end, _ in introner_clusters:
        new_remaining = []
        for region_start, region_end in remaining_regions:
            # No overlap
            if cluster_end <= region_start or cluster_start >= region_end:
                new_remaining.append((region_start, region_end))
            else:
                # Overlap - split region around cluster
                if region_start < cluster_start:
                    new_remaining.append((region_start, cluster_start))
                if region_end > cluster_end:
                    new_remaining.append((cluster_end, region_end))
        remaining_regions = new_remaining

    # Tile non-introner regions with 10kb windows
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


def main():
    parser = argparse.ArgumentParser(
        description="Group 2 gene-centric, introner-centered recombination analysis"
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
        '--exclude-contig',
        action='append',
        default=[],
        help=(
            'Contig to exclude from the analysis. May be supplied more than once; '
            'full and normalized RCC1749 contig names are both accepted.'
        )
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=10000,
        help='Window size in bp (default: 10000)'
    )
    parser.add_argument(
        '--merge_distance',
        type=int,
        default=10000,
        help='Distance to merge nearby introners (default: 10000)'
    )
    parser.add_argument(
        '--output_prefix',
        default='gene_exonic_group2_introners',
        help='Output file prefix'
    )
    parser.add_argument('--exclude-start', type=int, default=None, help='1-based inclusive start of the excluded interval on --exclude-contig.')
    parser.add_argument('--exclude-end', type=int, default=None, help='1-based inclusive end of the excluded interval on --exclude-contig.')
    args = parser.parse_args()
    if (args.exclude_start is None) != (args.exclude_end is None):
        parser.error('--exclude-start and --exclude-end must be supplied together')
    if args.exclude_start is not None and not 1 <= args.exclude_start <= args.exclude_end:
        parser.error('Invalid excluded interval')

    excluded_contigs = {
        normalize_contig_name(contig) for contig in args.exclude_contig
    }
    print("📂 Loading gene annotations from GTF...")
    genes = load_gtf(args.gtf_file)
    if excluded_contigs and args.exclude_start is None:
        genes = {gene_id: value for gene_id, value in genes.items()
                 if normalize_contig_name(value['contig']) not in excluded_contigs}
    print(f"\n📂 Extracting Group 2 introners from genotype matrix...")
    introners = extract_group2_introners(args.introner_matrix)
    if excluded_contigs:
        def excluded(contig, start, end):
            return normalize_contig_name(contig) in excluded_contigs and (
                args.exclude_start is None or (start < args.exclude_end and end > args.exclude_start - 1))
        introners = {key: value for key, value in introners.items()
                     if not excluded(value['contig'], value['start'], value['end'])}
    print("\n📊 Clustering introners per gene...")
    gene_introner_clusters = {}
    n_genes_with_introners = 0
    for gene_id, gene_data in genes.items():
        clusters = cluster_introners(introners, gene_id, gene_data, args.merge_distance)
        if len(clusters) > 0:
            n_genes_with_introners += 1
        gene_introner_clusters[gene_id] = clusters

    print(f"   {n_genes_with_introners} genes contain Group 2 introners")

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
        # Strip the pyrho.out extension
        chrom_full = basename.replace('.pyrho.out', '')
        # For RCC1749 pyrho files named like "RCC1749_0_intronerless_contig_10"
        # Extract just the contig part (everything after the second underscore)
        parts = chrom_full.split('_')
        if len(parts) > 2 and parts[0] == 'RCC1749':
            # Join parts from index 2 onwards (skip 'RCC1749' and '0')
            chrom = '_'.join(parts[2:])
        else:
            chrom = chrom_full

        pyrho_data[chrom] = parse_pyrho_file(pyrho_file)
        print(f"   Loaded {chrom}: {len(pyrho_data[chrom])} pyrho windows")

    print("\n⚡ Calculating recombination rates...")

    # We need to add contig info to windows
    for window in all_windows:
        gene_id = window['gene_id']
        gene_data = genes[gene_id]
        contig = normalize_contig_name(gene_data['contig'])

        if contig not in pyrho_data:
            window['rate'] = np.nan
            continue

        pyrho_windows = pyrho_data[contig]
        rate = calculate_weighted_mean_rate(window['start'], window['end'], pyrho_windows)
        window['rate'] = rate
        window['contig'] = contig

    # Filter out windows with NaN rates
    windows_df = pd.DataFrame(all_windows)
    windows_df = windows_df[~windows_df['rate'].isna()].copy()

    if excluded_contigs and 'contig' in windows_df.columns:
        excluded_output = windows_df['contig'].isin(excluded_contigs)
        if args.exclude_start is not None:
            excluded_output &= windows_df.start.lt(args.exclude_end) & windows_df.end.gt(args.exclude_start - 1)
        windows_df = windows_df.loc[~excluded_output].copy()
    print(f"   Calculated rates for {len(windows_df)} windows")

    # Summary statistics
    print("\n📈 Summary Statistics:")
    print(f"\n{'Window Type':<30} {'N Windows':<12} {'Mean Rate':<20} {'Median Rate':<20}")
    print("-" * 85)

    stats_summary = []
    for wtype in ['non_introner_containing', 'introner_containing']:
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

    # Statistical test
    introner_windows = windows_df[windows_df['type'] == 'introner_containing']['rate'].values
    non_introner_windows = windows_df[windows_df['type'] == 'non_introner_containing']['rate'].values

    pvalue = None
    statistic = None

    print(f"\n🧪 Mann-Whitney U Test:")
    if len(introner_windows) > 0 and len(non_introner_windows) > 0:
        statistic, pvalue = mannwhitneyu(introner_windows, non_introner_windows, alternative='two-sided')

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
        print(f"   ⚠️  Cannot perform test (insufficient data)")
        print(f"   Introner windows: {len(introner_windows)}, Non-introner: {len(non_introner_windows)}")

    # Save outputs
    summary_df = pd.DataFrame(stats_summary)
    summary_file = f"{args.output_prefix}_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n💾 Saved summary to {summary_file}")

    detail_file = f"{args.output_prefix}_windows.tsv"
    windows_df.to_csv(detail_file, sep='\t', index=False)
    print(f"💾 Saved window details to {detail_file}")

    # Create visualization
    print("\n📊 Creating visualization...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Box plot
    sns.boxplot(
        data=windows_df,
        x='type',
        y='rate',
        ax=axes[0],
        showfliers=False,
        palette=['lightblue', 'lightcoral']
    )
    axes[0].set_ylabel('Recombination Rate (r)', fontsize=11)
    axes[0].set_xlabel('')
    title = 'Group 2 Exonic Windows: Introner-containing vs Non-introner'
    if pvalue is not None:
        title += f'\np={pvalue:.2e}'
    axes[0].set_title(title, fontsize=12, fontweight='bold')
    axes[0].set_yscale('log')
    axes[0].set_xticklabels(['Non-introner', 'Introner'])

    # Histogram
    for wtype, color, label in [('non_introner_containing', 'blue', 'Non-introner'),
                                  ('introner_containing', 'red', 'Introner')]:
        subset = windows_df[windows_df['type'] == wtype]['rate']
        axes[1].hist(np.log10(subset), bins=40, alpha=0.5, label=label, color=color)

    axes[1].set_xlabel('log10(Recombination Rate)', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('Rate Distribution (log scale)', fontsize=12)
    axes[1].legend()

    plt.suptitle(f'Group 2 Gene-Exonic Recombination Analysis\n{args.window_size}bp windows, {args.merge_distance}bp merge distance',
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    save_plot_bundle(args.output_prefix)
    print(f"Saved plots to {Path(args.output_prefix).with_suffix('.pdf')} and {Path(args.output_prefix).with_suffix('.png')}")

    print("\n✅ Analysis complete!")


if __name__ == '__main__':
    main()

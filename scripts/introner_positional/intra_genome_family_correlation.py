#!/usr/bin/env python3
"""
Intra-genome introner family correlation analysis.

Analyzes the correlation between the spatial distributions of two introner families
within a single genome. This helps determine whether introner families show
coordinated or independent spatial patterns.

Presence encoding in genotype matrix:
  1 = introner present (counted)
  2 = introner absent (not counted)
  3 = not callable / insufficient coverage (not counted, treated as absent)

Coordinate note:
  Genotype matrix coordinates are 0-based with 100bp flanking on each side.
  This script removes the flanking and converts to 1-based.

Usage:
    python intra_genome_family_correlation.py \\
        --genotype-matrix genotype_matrix.final.tsv \\
        --genome-fai assemblies/CCMP1545.vg_paths.fa.fai \\
        --genome-name CCMP1545 \\
        --family1 2 --family2 4 \\
        --window-size 10000 \\
        --output results.json --plot results_hexbin.png
"""

import argparse
import json
import sys
from bisect import bisect_right
from collections import defaultdict
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate intra-genome Spearman correlation between introner families',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--genotype-matrix', required=True,
                       help='Genotype matrix TSV file with all introner positions')
    parser.add_argument('--genome-fai', required=True,
                       help='FASTA index file for the genome')
    parser.add_argument('--genome-name', required=True,
                       help='Sample/genome name (e.g., CCMP1545, RCC1749)')
    parser.add_argument('--family1', type=int, default=2,
                       help='First introner family to compare (default: 2)')
    parser.add_argument('--family2', type=int, default=4,
                       help='Second introner family to compare (default: 4)')
    parser.add_argument('--window-size', type=int, default=10000,
                       help='Window size in base pairs (default: 10000)')
    parser.add_argument('--exclude-fixed', action='store_true',
                       help='Exclude fixed introners (present in both CCMP1545 and RCC1749)')
    parser.add_argument('--output', required=True,
                       help='Output JSON file for detailed results')
    parser.add_argument('--plot', required=True,
                       help='Output PNG file for hexbin plot')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip generating the hexbin plot')

    return parser.parse_args()


def parse_fai_file(fai_path):
    """Parse FASTA index file to get chromosome/scaffold lengths.

    Normalizes chromosome names (# -> _) to match genotype matrix convention.
    """
    genome_info = {}
    with open(fai_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom_name = normalize_chrom_name(fields[0])
            genome_info[chrom_name] = int(fields[1])
    return genome_info


def normalize_chrom_name(chrom_name):
    """Normalize chromosome names (# -> _ separator)."""
    return chrom_name.replace('#', '_')


def parse_introner_positions_from_genotype_matrix(genotype_matrix_file, sample_name,
                                                   family, exclude_fixed=False):
    """
    Parse genotype matrix to extract introner positions for a specific sample and family.

    Only rows with presence == 1 are counted. Presence 2 (absent) and 3 (not callable)
    are skipped.
    """
    print(f"Parsing introners for {sample_name}, family {family}...")

    df = pd.read_csv(genotype_matrix_file, sep='\t')

    # Identify fixed introners if needed
    fixed_ortholog_ids = set()
    if exclude_fixed:
        print(f"  Identifying fixed introners to exclude...")
        for ortholog_id, group in df.groupby('ortholog_id'):
            ccmp1545_rows = group[group['sample'] == 'CCMP1545']
            rcc1749_rows = group[group['sample'] == 'RCC1749']

            if (len(ccmp1545_rows) > 0 and len(rcc1749_rows) > 0 and
                (ccmp1545_rows['presence'] == 1).any() and
                (rcc1749_rows['presence'] == 1).any()):
                fixed_ortholog_ids.add(ortholog_id)

        print(f"  Found {len(fixed_ortholog_ids)} fixed introners")

    # Filter for this sample and presence == 1
    sample_df = df[(df['sample'] == sample_name) & (df['presence'] == 1)]

    introners = []
    for _, row in sample_df.iterrows():
        if row['start'] == -1 or row['end'] == -1:
            continue

        if exclude_fixed and row['ortholog_id'] in fixed_ortholog_ids:
            continue

        # Check family match
        try:
            introner_family = int(row['family'])
        except (ValueError, TypeError):
            continue

        if introner_family != family:
            continue

        chrom = normalize_chrom_name(row['contig'])

        # Remove 100bp flanking and convert 0-based to 1-based
        start = row['start'] + 100 + 1
        end = row['end'] - 100

        gene = row['gene'] if pd.notna(row['gene']) and row['gene'] != '' else 'intergenic'

        introners.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'id': row['ortholog_id'],
            'gene': gene
        })

    print(f"  Found {len(introners)} introners present in {sample_name}, family {family}")
    return introners


def build_introner_index(introners):
    """Build sorted index of introner positions per chromosome for fast overlap queries."""
    index = defaultdict(list)
    for intr in introners:
        index[intr['chrom']].append((intr['start'], intr['end']))
    for chrom in index:
        index[chrom].sort()
    return index


def count_introners_in_region(index, chrom, start, end):
    """Count introners overlapping a region using bisect for efficiency."""
    positions = index.get(chrom, [])
    if not positions:
        return 0

    right_idx = bisect_right(positions, (end, float('inf')))

    count = 0
    for i in range(right_idx - 1, -1, -1):
        intr_start, intr_end = positions[i]
        if intr_end < start:
            if intr_start < start - 10000:
                break
            continue
        count += 1

    return count


def create_genome_windows(genome_info, window_size):
    """Create fixed-size non-overlapping windows across the genome (1-based)."""
    windows = []
    for chrom, chrom_length in sorted(genome_info.items()):
        for window_start in range(1, chrom_length + 1, window_size):
            window_end = min(window_start + window_size - 1, chrom_length)
            windows.append({
                'chrom': chrom,
                'start': window_start,
                'end': window_end,
                'window_id': f"{chrom}:{window_start}-{window_end}",
                'length': window_end - window_start + 1
            })
    return windows


def calculate_window_family_densities(windows, family1_index, family2_index):
    """Calculate family densities for windows in a single genome."""
    window_data = []

    print(f"\nCalculating family densities for {len(windows)} windows...")

    for i, window in enumerate(windows):
        if (i + 1) % 500 == 0:
            print(f"  Processed {i + 1}/{len(windows)} windows...")

        family1_count = count_introners_in_region(
            family1_index, window['chrom'], window['start'], window['end']
        )
        family2_count = count_introners_in_region(
            family2_index, window['chrom'], window['start'], window['end']
        )

        family1_density = (family1_count / window['length']) * 1000
        family2_density = (family2_count / window['length']) * 1000

        window_data.append({
            'window_id': window['window_id'],
            'chrom': window['chrom'],
            'start': window['start'],
            'end': window['end'],
            'family1_count': family1_count,
            'family1_density': family1_density,
            'family2_count': family2_count,
            'family2_density': family2_density
        })

    # Log summary stats
    n_both_zero = sum(1 for w in window_data
                      if w['family1_count'] == 0 and w['family2_count'] == 0)
    n_either = sum(1 for w in window_data
                   if w['family1_count'] > 0 or w['family2_count'] > 0)
    print(f"  Processed {len(windows)} windows")
    print(f"  Windows with at least one introner: {n_either}")
    print(f"  Windows with zero introners from both families: {n_both_zero}")

    return window_data


def calculate_spearman_correlation(window_data):
    """Calculate Spearman correlation between family densities."""
    family1_densities = [w['family1_density'] for w in window_data]
    family2_densities = [w['family2_density'] for w in window_data]

    correlation, p_value = spearmanr(family1_densities, family2_densities)

    return {
        'correlation': correlation,
        'p_value': p_value,
        'n_windows': len(window_data),
        'family1_densities': family1_densities,
        'family2_densities': family2_densities
    }


def create_hexbin_plot(results, window_data, output_file, window_size, family1, family2):
    """Create hexbin plot showing correlation between family densities."""
    family1_densities = results['family1_densities']
    family2_densities = results['family2_densities']

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    hb = ax.hexbin(family1_densities, family2_densities, gridsize=30, cmap='YlOrRd', mincnt=1)

    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Number of Windows', fontsize=12)

    max_family1 = max(family1_densities) if family1_densities else 1
    max_family2 = max(family2_densities) if family2_densities else 1
    max_val = max(max_family1, max_family2)

    ax.set_xlim(-0.05, max_val * 1.05)
    ax.set_ylim(-0.05, max_val * 1.05)

    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1,
            label='Perfect correlation')

    ax.set_xlabel(f'Family {family1} Introner Density (introners/kb)', fontsize=12)
    ax.set_ylabel(f'Family {family2} Introner Density (introners/kb)', fontsize=12)
    ax.set_title(f'Intra-genome Family {family1} vs {family2} Correlation ({window_size/1000:.0f}kb windows)\n'
                f'Spearman \u03c1 = {results["correlation"]:.3f}, p = {results["p_value"]:.2e}',
                fontsize=14, pad=20)

    ax.grid(True, alpha=0.3)
    ax.legend()

    n_zero_zero = sum(1 for i in range(len(family1_densities))
                     if family1_densities[i] == 0 and family2_densities[i] == 0)

    stats_text = (f'N = {results["n_windows"]:,} windows\n'
                 f'(0,0) density windows: {n_zero_zero:,} ({100*n_zero_zero/results["n_windows"]:.1f}%)\n'
                 f'Family {family1} mean density: {np.mean(family1_densities):.4f}/kb\n'
                 f'Family {family2} mean density: {np.mean(family2_densities):.4f}/kb')

    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nHexbin plot saved to: {output_file}")


def save_detailed_results(results, window_data, output_file, window_size, family1, family2, genome_name):
    """Save detailed results to JSON file."""
    output_data = {
        "metadata": {
            "analysis": "Intra-genome Introner Family Correlation Analysis",
            "genome": genome_name,
            "family1": family1,
            "family2": family2,
            "window_size_bp": window_size,
            "n_windows": results['n_windows']
        },
        "statistics": {
            "spearman_correlation": round(results['correlation'], 4),
            "p_value": results['p_value'],
            "family1_mean_density": round(np.mean(results['family1_densities']), 4),
            "family2_mean_density": round(np.mean(results['family2_densities']), 4),
            "family1_median_density": round(np.median(results['family1_densities']), 4),
            "family2_median_density": round(np.median(results['family2_densities']), 4)
        },
        "window_data": []
    }

    for window in window_data:
        output_data["window_data"].append({
            "window_id": window['window_id'],
            "chrom": window['chrom'],
            "start": window['start'],
            "end": window['end'],
            "family1_count": window['family1_count'],
            "family1_density": round(window['family1_density'], 4),
            "family2_count": window['family2_count'],
            "family2_density": round(window['family2_density'], 4)
        })

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"Detailed results saved to: {output_file}")


def print_summary_statistics(results, family1, family2, genome_name):
    """Print summary statistics."""
    print("\n" + "="*80)
    print(f"INTRA-GENOME INTRONER FAMILY CORRELATION ANALYSIS")
    print(f"Genome: {genome_name} | Comparing Family {family1} vs Family {family2}")
    print("="*80)

    print(f"\nDataset Summary:")
    print(f"  Number of windows analyzed: {results['n_windows']:,}")

    family1_densities = results['family1_densities']
    family2_densities = results['family2_densities']

    print(f"\nFamily {family1} introner density distribution:")
    print(f"  Mean: {np.mean(family1_densities):.4f} introners/kb")
    print(f"  Median: {np.median(family1_densities):.4f} introners/kb")
    print(f"  Min: {min(family1_densities):.4f}, Max: {max(family1_densities):.4f}")
    print(f"  Windows with 0 introners: {family1_densities.count(0)} "
          f"({100*family1_densities.count(0)/len(family1_densities):.1f}%)")

    print(f"\nFamily {family2} introner density distribution:")
    print(f"  Mean: {np.mean(family2_densities):.4f} introners/kb")
    print(f"  Median: {np.median(family2_densities):.4f} introners/kb")
    print(f"  Min: {min(family2_densities):.4f}, Max: {max(family2_densities):.4f}")
    print(f"  Windows with 0 introners: {family2_densities.count(0)} "
          f"({100*family2_densities.count(0)/len(family2_densities):.1f}%)")

    print(f"\nSpearman's Rank Correlation:")
    print(f"  Correlation coefficient (\u03c1): {results['correlation']:.4f}")
    print(f"  P-value: {results['p_value']:.2e}")

    if results['p_value'] < 0.001:
        significance = "highly significant (p < 0.001)"
    elif results['p_value'] < 0.01:
        significance = "very significant (p < 0.01)"
    elif results['p_value'] < 0.05:
        significance = "significant (p < 0.05)"
    else:
        significance = "not significant (p >= 0.05)"

    print(f"  Statistical significance: {significance}")

    if results['correlation'] > 0.5:
        pattern = "strong positive correlation (families show similar spatial patterns)"
    elif results['correlation'] > 0.2:
        pattern = "moderate positive correlation"
    elif results['correlation'] > -0.2:
        pattern = "weak/no correlation (families distribute independently)"
    elif results['correlation'] > -0.5:
        pattern = "moderate negative correlation"
    else:
        pattern = "strong negative correlation (families exclude each other)"

    print(f"  Interpretation: {pattern}")


def main():
    """Main function."""
    args = parse_arguments()

    print("="*80)
    print("INTRA-GENOME INTRONER FAMILY CORRELATION ANALYSIS")
    print("="*80)

    # Parse genome information
    print("\nLoading genome information...")
    genome_info = parse_fai_file(args.genome_fai)
    total_bp = sum(genome_info.values())
    print(f"  {args.genome_name}: {total_bp:,} bp across {len(genome_info)} sequences")

    # Parse introner positions
    family1_introners = parse_introner_positions_from_genotype_matrix(
        args.genotype_matrix, args.genome_name, args.family1, args.exclude_fixed
    )
    family2_introners = parse_introner_positions_from_genotype_matrix(
        args.genotype_matrix, args.genome_name, args.family2, args.exclude_fixed
    )

    # Build spatial indices
    print("\nBuilding spatial indices...")
    family1_index = build_introner_index(family1_introners)
    family2_index = build_introner_index(family2_introners)

    # Create windows
    print(f"\nCreating {args.window_size}bp windows...")
    windows = create_genome_windows(genome_info, args.window_size)
    print(f"  Created {len(windows)} windows")

    # Calculate family densities
    window_data = calculate_window_family_densities(windows, family1_index, family2_index)

    if len(window_data) < 10:
        print(f"\nError: Only {len(window_data)} windows found. "
              "Try increasing --window-size", file=sys.stderr)
        sys.exit(1)

    # Calculate Spearman correlation
    print("\nCalculating Spearman correlation...")
    results = calculate_spearman_correlation(window_data)

    # Print results
    print_summary_statistics(results, args.family1, args.family2, args.genome_name)

    # Save detailed results
    save_detailed_results(results, window_data, args.output,
                         args.window_size, args.family1, args.family2, args.genome_name)

    # Generate hexbin plot
    if not args.no_plot:
        create_hexbin_plot(results, window_data, args.plot,
                          args.window_size, args.family1, args.family2)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == "__main__":
    main()

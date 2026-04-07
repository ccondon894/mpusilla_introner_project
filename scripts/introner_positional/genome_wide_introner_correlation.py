#!/usr/bin/env python3
"""
Genome-wide introner density correlation analysis between two strains.

Performs a position-based Spearman rank correlation of introner density
across the entire genome using MUMmer alignment blocks to map corresponding
regions between two assemblies.

Approach:
1. Parse MUMmer alignments to identify syntenic regions
2. Divide reference genome into fixed-size windows
3. Map corresponding windows to query genome using alignment blocks
4. Calculate introner density in aligned window pairs
5. Generate hexbin scatter plot and calculate Spearman correlation

Presence encoding in genotype matrix:
  1 = introner present (counted)
  2 = introner absent (not counted)
  3 = not callable / insufficient coverage (not counted, treated as absent)

Coordinate note:
  Genotype matrix coordinates are 0-based with 100bp flanking on each side.
  This script removes the flanking and converts to 1-based for alignment with
  FAI-based genomic coordinates.

Usage:
    python genome_wide_introner_correlation.py \\
        --delta mummer-output/CCMP1545_vs_RCC1749.1delta \\
        --genotype-matrix genotype_matrix.final.tsv \\
        --g1-fai assemblies/CCMP1545.vg_paths.fa.fai \\
        --g2-fai assemblies/RCC1749.vg_paths.fa.fai \\
        --g1-name CCMP1545 --g2-name RCC1749 \\
        --window-size 10000 \\
        --output results.json --plot results_hexbin.png
"""

import argparse
import json
import sys
from bisect import bisect_left, bisect_right
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
        description='Calculate genome-wide Spearman correlation of introner densities',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--delta', required=True,
                       help='MUMmer delta file for genome alignment')
    parser.add_argument('--genotype-matrix', required=True,
                       help='Genotype matrix TSV file with all introner positions')
    parser.add_argument('--g1-fai', required=True,
                       help='FASTA index file for reference (group 1) genome')
    parser.add_argument('--g2-fai', required=True,
                       help='FASTA index file for query (group 2) genome')
    parser.add_argument('--g1-name', required=True,
                       help='Reference strain name (e.g., CCMP1545)')
    parser.add_argument('--g2-name', required=True,
                       help='Query strain name (e.g., RCC1749)')
    parser.add_argument('--window-size', type=int, default=10000,
                       help='Window size in base pairs (default: 10000)')
    parser.add_argument('--introner-family', type=str, default='all', dest='family',
                       help='Filter introners by family (default: all)')
    parser.add_argument('--min-alignment-coverage', type=float, default=0.5,
                       help='Minimum fraction of window that must be aligned (default: 0.5)')
    parser.add_argument('--exclude-fixed', action='store_true',
                       help='Exclude fixed introners (present in both strains)')
    parser.add_argument('--output', required=True,
                       help='Output JSON file for detailed results')
    parser.add_argument('--plot', required=True,
                       help='Output PNG file for hexbin plot')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip generating the hexbin plot')

    return parser.parse_args()


def parse_fai_file(fai_path):
    """Parse FASTA index file to get chromosome/scaffold lengths."""
    genome_info = {}
    with open(fai_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom_name = normalize_chrom_name(fields[0])
            chrom_length = int(fields[1])
            genome_info[chrom_name] = chrom_length
    return genome_info


def parse_delta_file(delta_path):
    """
    Parse MUMmer delta file to extract alignment blocks.

    Returns list of alignment blocks with ref/query coordinates.
    Note: Proportional mapping is used to translate reference positions
    to query positions within each alignment block. This is an approximation
    that assumes roughly collinear alignment within each block.
    """
    alignments = []

    print(f"Parsing MUMmer delta file: {delta_path}")

    with open(delta_path, 'r') as f:
        # Skip header lines
        f.readline()
        f.readline()

        current_ref_seq = None
        current_query_seq = None
        reading_alignment_block = False

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                fields = line[1:].split()
                current_ref_seq = normalize_chrom_name(fields[0])
                current_query_seq = normalize_chrom_name(fields[1])
                reading_alignment_block = False

            elif line == "0":
                reading_alignment_block = False

            elif not reading_alignment_block:
                fields = line.split()
                if len(fields) >= 7:
                    try:
                        block = {
                            'ref_seq': current_ref_seq,
                            'query_seq': current_query_seq,
                            'ref_start': int(fields[0]),
                            'ref_end': int(fields[1]),
                            'query_start': int(fields[2]),
                            'query_end': int(fields[3]),
                            'errors': int(fields[4]),
                            'matches': int(fields[5]),
                            'indels': int(fields[6])
                        }

                        # Normalize coordinates (always start < end)
                        if block['ref_start'] > block['ref_end']:
                            block['ref_start'], block['ref_end'] = block['ref_end'], block['ref_start']
                            block['ref_reversed'] = True
                        else:
                            block['ref_reversed'] = False

                        if block['query_start'] > block['query_end']:
                            block['query_start'], block['query_end'] = block['query_end'], block['query_start']
                            block['query_reversed'] = True
                        else:
                            block['query_reversed'] = False

                        alignments.append(block)
                        reading_alignment_block = True

                    except (ValueError, IndexError):
                        continue

    print(f"  Found {len(alignments)} alignment blocks")
    return alignments


def normalize_chrom_name(chrom_name):
    """
    Normalize chromosome names to match between genotype matrix and FASTA files.
    Genotype matrix uses '#' as separator (e.g., CCMP1545#0#scaffold_1)
    FASTA uses '_' as separator (e.g., CCMP1545_0_scaffold_1)
    """
    return chrom_name.replace('#', '_')


def parse_introner_positions_from_genotype_matrix(genotype_matrix_file, sample_name,
                                                    g1_name, g2_name,
                                                    exclude_fixed=False, family='all'):
    """
    Parse genotype matrix to extract introner positions for a specific sample.

    Coordinates in the matrix are 0-based and include 100bp flanking sequence
    on each side. We remove flanking and convert to 1-based.

    Only rows with presence == 1 are counted as present introners.
    Rows with presence == 2 (absent) or 3 (not callable) are skipped.
    """
    print(f"Parsing introners for {sample_name} from genotype matrix...")

    df = pd.read_csv(genotype_matrix_file, sep='\t')

    # Filter by family if specified
    if family != 'all':
        try:
            family_int = int(family)
            df = df[df['family'] == family_int]
        except (ValueError, TypeError):
            df = df[df['family'] == family]

    # Identify fixed introners (present in both reference strains) if needed
    fixed_ortholog_ids = set()
    if exclude_fixed:
        print(f"  Identifying fixed introners to exclude...")
        for ortholog_id, group in df.groupby('ortholog_id'):
            g1_rows = group[group['sample'] == g1_name]
            g2_rows = group[group['sample'] == g2_name]

            if (len(g1_rows) > 0 and len(g2_rows) > 0 and
                (g1_rows['presence'] == 1).any() and
                (g2_rows['presence'] == 1).any()):
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

    print(f"  Found {len(introners)} introners present in {sample_name}")
    return introners


def build_introner_index(introners):
    """
    Build a sorted index of introner positions per chromosome for fast overlap queries.

    Uses bisect-based lookups instead of O(n) linear scans.
    For each chromosome, stores sorted list of (start, end) tuples.
    """
    index = defaultdict(list)
    for intr in introners:
        index[intr['chrom']].append((intr['start'], intr['end']))

    # Sort each chromosome's introners by start position
    for chrom in index:
        index[chrom].sort()

    return index


def count_introners_in_region(index, chrom, start, end):
    """
    Count introners overlapping a region using bisect for efficiency.

    An introner overlaps [start, end] if introner.start <= end AND introner.end >= start.
    We use bisect to find candidates whose start <= end, then check the end >= start condition.
    """
    positions = index.get(chrom, [])
    if not positions:
        return 0

    # Find rightmost introner whose start <= end
    # All introners with index < right_idx have start <= end
    right_idx = bisect_right(positions, (end, float('inf')))

    count = 0
    # Check candidates from the end backwards; stop when start is too far left
    for i in range(right_idx - 1, -1, -1):
        intr_start, intr_end = positions[i]
        if intr_end < start:
            # Since list is sorted by start, all earlier entries also have start < intr_start
            # But their end could still overlap... however if this introner's end < region start
            # and we're going backwards, earlier introners have smaller starts.
            # We can't break early on end alone since ends aren't sorted.
            # But we CAN break if intr_start + max_introner_length < start.
            # For safety, just continue checking.
            # Actually, for typical introner sizes (~200bp), if intr_start is far below
            # start, we can break. Use a generous cutoff.
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


def calculate_overlap(range1_start, range1_end, range2_start, range2_end):
    """Calculate overlap between two ranges."""
    overlap_start = max(range1_start, range2_start)
    overlap_end = min(range1_end, range2_end)
    return max(0, overlap_end - overlap_start + 1)


def map_window_to_aligned_region(window, alignments, min_coverage):
    """
    Map a reference window to aligned query regions using proportional mapping.

    For each alignment block overlapping the window, the overlapping portion is
    proportionally mapped to query coordinates. This is an approximation that
    assumes roughly linear correspondence within each alignment block.
    """
    window_chrom = window['chrom']
    window_start = window['start']
    window_end = window['end']
    window_length = window['length']

    overlapping_alignments = []
    for aln in alignments:
        if aln['ref_seq'] != window_chrom:
            continue

        overlap = calculate_overlap(window_start, window_end,
                                    aln['ref_start'], aln['ref_end'])
        if overlap > 0:
            overlapping_alignments.append({
                'alignment': aln,
                'overlap': overlap
            })

    if not overlapping_alignments:
        return []

    total_coverage = sum(oa['overlap'] for oa in overlapping_alignments)
    coverage_fraction = total_coverage / window_length

    if coverage_fraction < min_coverage:
        return []

    query_regions = []
    for oa in overlapping_alignments:
        aln = oa['alignment']

        overlap_ref_start = max(window_start, aln['ref_start'])
        overlap_ref_end = min(window_end, aln['ref_end'])

        ref_rel_start = overlap_ref_start - aln['ref_start']
        ref_rel_end = overlap_ref_end - aln['ref_start']

        aln_ref_length = aln['ref_end'] - aln['ref_start']
        aln_query_length = aln['query_end'] - aln['query_start']

        if aln_ref_length > 0:
            query_rel_start = int(ref_rel_start * aln_query_length / aln_ref_length)
            query_rel_end = int(ref_rel_end * aln_query_length / aln_ref_length)

            query_start = aln['query_start'] + query_rel_start
            query_end = aln['query_start'] + query_rel_end

            query_regions.append({
                'query_seq': aln['query_seq'],
                'query_start': query_start,
                'query_end': query_end,
                'coverage': oa['overlap'] / window_length
            })

    return query_regions


def calculate_window_introner_densities(ref_windows, alignments, ref_introner_index,
                                        query_introner_index, min_coverage):
    """
    Calculate introner densities for aligned window pairs.

    Returns list of window pair data with densities.
    Also logs how many windows were dropped due to insufficient alignment coverage.
    """
    window_data = []
    n_no_alignment = 0
    n_low_coverage = 0

    print(f"\nMapping {len(ref_windows)} reference windows to query regions...")

    for i, ref_window in enumerate(ref_windows):
        if (i + 1) % 500 == 0:
            print(f"  Processed {i + 1}/{len(ref_windows)} windows...")

        query_regions = map_window_to_aligned_region(ref_window, alignments, min_coverage)

        if not query_regions:
            # Check if there was any alignment at all
            has_any = any(
                aln['ref_seq'] == ref_window['chrom'] and
                calculate_overlap(ref_window['start'], ref_window['end'],
                                  aln['ref_start'], aln['ref_end']) > 0
                for aln in alignments
            )
            if has_any:
                n_low_coverage += 1
            else:
                n_no_alignment += 1
            continue

        ref_introner_count = count_introners_in_region(
            ref_introner_index, ref_window['chrom'],
            ref_window['start'], ref_window['end']
        )

        query_introner_count = 0
        total_query_length = 0
        for query_region in query_regions:
            query_count = count_introners_in_region(
                query_introner_index, query_region['query_seq'],
                query_region['query_start'], query_region['query_end']
            )
            query_introner_count += query_count
            total_query_length += query_region['query_end'] - query_region['query_start'] + 1

        ref_density = (ref_introner_count / ref_window['length']) * 1000
        query_density = (query_introner_count / total_query_length) * 1000 if total_query_length > 0 else 0

        window_data.append({
            'ref_window_id': ref_window['window_id'],
            'ref_chrom': ref_window['chrom'],
            'ref_start': ref_window['start'],
            'ref_end': ref_window['end'],
            'ref_introner_count': ref_introner_count,
            'ref_density': ref_density,
            'query_introner_count': query_introner_count,
            'query_density': query_density,
            'num_query_regions': len(query_regions)
        })

    n_total = len(ref_windows)
    n_kept = len(window_data)
    n_dropped = n_total - n_kept
    print(f"  Kept {n_kept} aligned window pairs out of {n_total} total windows")
    print(f"  Filtered out {n_dropped} windows ({100*n_dropped/n_total:.1f}%):")
    print(f"    No alignment overlap: {n_no_alignment}")
    print(f"    Below {min_coverage*100:.0f}% alignment coverage: {n_low_coverage}")

    return window_data


def calculate_spearman_correlation(window_data):
    """Calculate Spearman correlation for window densities."""
    ref_densities = [w['ref_density'] for w in window_data]
    query_densities = [w['query_density'] for w in window_data]

    correlation, p_value = spearmanr(ref_densities, query_densities)

    return {
        'correlation': correlation,
        'p_value': p_value,
        'n_windows': len(window_data),
        'ref_densities': ref_densities,
        'query_densities': query_densities
    }


def create_hexbin_plot(results, window_data, output_file, window_size, g1_name, g2_name):
    """Create hexbin plot showing correlation between introner densities."""
    ref_densities = results['ref_densities']
    query_densities = results['query_densities']

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    hb = ax.hexbin(ref_densities, query_densities, gridsize=30, cmap='YlOrRd', mincnt=1)

    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Number of Windows', fontsize=12)

    max_ref = max(ref_densities) if ref_densities else 1
    max_query = max(query_densities) if query_densities else 1
    max_val = max(max_ref, max_query)

    ax.set_xlim(-0.05, max_val * 1.05)
    ax.set_ylim(-0.05, max_val * 1.05)

    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1,
            label='Perfect correlation')

    ax.set_xlabel(f'{g1_name} Introner Density (introners/kb)', fontsize=12)
    ax.set_ylabel(f'{g2_name} Introner Density (introners/kb)', fontsize=12)
    ax.set_title(f'Genome-wide Introner Density Correlation ({window_size/1000:.0f}kb windows)\n'
                f'Spearman \u03c1 = {results["correlation"]:.3f}, p = {results["p_value"]:.2e}',
                fontsize=14, pad=20)

    ax.grid(True, alpha=0.3)
    ax.legend()

    n_zero_zero = sum(1 for i in range(len(ref_densities))
                     if ref_densities[i] == 0 and query_densities[i] == 0)

    stats_text = (f'N = {results["n_windows"]:,} window pairs\n'
                 f'(0,0) density windows: {n_zero_zero:,} ({100*n_zero_zero/results["n_windows"]:.1f}%)\n'
                 f'{g1_name} mean density: {np.mean(ref_densities):.4f}/kb\n'
                 f'{g2_name} mean density: {np.mean(query_densities):.4f}/kb')

    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nHexbin plot saved to: {output_file}")


def save_detailed_results(results, window_data, output_file, window_size, min_coverage):
    """Save detailed results to JSON file."""
    output_data = {
        "metadata": {
            "analysis": "Genome-wide Spearman Rank Correlation Analysis: Introner Density",
            "window_size_bp": window_size,
            "min_alignment_coverage": min_coverage,
            "n_windows": results['n_windows']
        },
        "statistics": {
            "spearman_correlation": round(results['correlation'], 4),
            "p_value": results['p_value'],
            "ref_mean_density": round(np.mean(results['ref_densities']), 4),
            "query_mean_density": round(np.mean(results['query_densities']), 4),
            "ref_median_density": round(np.median(results['ref_densities']), 4),
            "query_median_density": round(np.median(results['query_densities']), 4)
        },
        "window_data": []
    }

    for window in window_data:
        output_data["window_data"].append({
            "ref_window_id": window['ref_window_id'],
            "ref_chrom": window['ref_chrom'],
            "ref_start": window['ref_start'],
            "ref_end": window['ref_end'],
            "ref_introner_count": window['ref_introner_count'],
            "ref_density": round(window['ref_density'], 4),
            "query_introner_count": window['query_introner_count'],
            "query_density": round(window['query_density'], 4),
            "num_query_regions": window['num_query_regions']
        })

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"Detailed results saved to: {output_file}")


def print_summary_statistics(results, g1_name, g2_name):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("GENOME-WIDE INTRONER DENSITY CORRELATION ANALYSIS")
    print("="*80)

    print(f"\nDataset Summary:")
    print(f"  Number of aligned window pairs: {results['n_windows']:,}")

    ref_densities = results['ref_densities']
    query_densities = results['query_densities']

    print(f"\n{g1_name} introner density distribution:")
    print(f"  Mean: {np.mean(ref_densities):.4f} introners/kb")
    print(f"  Median: {np.median(ref_densities):.4f} introners/kb")
    print(f"  Min: {min(ref_densities):.4f}, Max: {max(ref_densities):.4f}")
    print(f"  Windows with 0 introners: {ref_densities.count(0)} "
          f"({100*ref_densities.count(0)/len(ref_densities):.1f}%)")

    print(f"\n{g2_name} introner density distribution:")
    print(f"  Mean: {np.mean(query_densities):.4f} introners/kb")
    print(f"  Median: {np.median(query_densities):.4f} introners/kb")
    print(f"  Min: {min(query_densities):.4f}, Max: {max(query_densities):.4f}")
    print(f"  Windows with 0 introners: {query_densities.count(0)} "
          f"({100*query_densities.count(0)/len(query_densities):.1f}%)")

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


def main():
    """Main function."""
    args = parse_arguments()

    print("="*80)
    print("GENOME-WIDE INTRONER DENSITY CORRELATION ANALYSIS")
    print("="*80)

    # Parse genome information
    print("\nLoading genome information...")
    g1_genome = parse_fai_file(args.g1_fai)
    g2_genome = parse_fai_file(args.g2_fai)

    g1_total = sum(g1_genome.values())
    g2_total = sum(g2_genome.values())

    print(f"  {args.g1_name}: {g1_total:,} bp across {len(g1_genome)} sequences")
    print(f"  {args.g2_name}: {g2_total:,} bp across {len(g2_genome)} sequences")

    # Parse alignment data
    alignments = parse_delta_file(args.delta)

    # Parse introner positions
    g1_introners = parse_introner_positions_from_genotype_matrix(
        args.genotype_matrix, args.g1_name, args.g1_name, args.g2_name,
        args.exclude_fixed, family=args.family)
    g2_introners = parse_introner_positions_from_genotype_matrix(
        args.genotype_matrix, args.g2_name, args.g1_name, args.g2_name,
        args.exclude_fixed, family=args.family)

    # Build spatial indices for fast overlap queries
    print("\nBuilding spatial indices...")
    g1_index = build_introner_index(g1_introners)
    g2_index = build_introner_index(g2_introners)

    # Create windows on reference genome
    print(f"\nCreating {args.window_size}bp windows...")
    g1_windows = create_genome_windows(g1_genome, args.window_size)
    print(f"  Created {len(g1_windows)} windows for {args.g1_name}")

    # Calculate introner densities
    window_data = calculate_window_introner_densities(
        g1_windows, alignments, g1_index, g2_index, args.min_alignment_coverage
    )

    if len(window_data) < 10:
        print(f"\nError: Only {len(window_data)} aligned windows found. "
              "Try reducing --min-alignment-coverage or increasing --window-size",
              file=sys.stderr)
        sys.exit(1)

    # Calculate Spearman correlation
    print("\nCalculating Spearman correlation...")
    results = calculate_spearman_correlation(window_data)

    # Print results
    print_summary_statistics(results, args.g1_name, args.g2_name)

    # Save detailed results
    save_detailed_results(results, window_data, args.output,
                         args.window_size, args.min_alignment_coverage)

    # Generate hexbin plot
    if not args.no_plot:
        create_hexbin_plot(results, window_data, args.plot,
                          args.window_size, args.g1_name, args.g2_name)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == "__main__":
    main()

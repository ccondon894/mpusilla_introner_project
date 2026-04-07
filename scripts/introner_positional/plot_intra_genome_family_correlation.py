#!/usr/bin/env python3
"""
Generate a publication-quality figure comparing intra-genome family correlations
across different genomic window sizes for different genomes.
"""

import argparse
import json
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Plot intra-genome family correlations across window sizes',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--result-dir', required=True,
                       help='Directory containing JSON result files')
    parser.add_argument('--genomes', nargs='+', default=['CCMP1545', 'RCC1749'],
                       help='Genomes to include in plot (default: CCMP1545 RCC1749)')
    parser.add_argument('--family-pairs', nargs='+', default=['2_4'],
                       help='Family pairs as F1_F2 strings (default: 2_4)')
    parser.add_argument('--window-sizes', nargs='+', type=int, default=[1, 5, 10, 20, 50],
                       help='Window sizes in kb to include (default: 1 5 10 20 50)')
    parser.add_argument('--output', required=True,
                       help='Output PNG file')

    return parser.parse_args()


def load_result_file(result_file):
    """Load correlation data from JSON result file."""
    try:
        with open(result_file, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError:
        print(f"Warning: Could not parse JSON from {result_file}", file=sys.stderr)
        return None


def get_result_file_path(result_dir, genome, family1, family2, window_size):
    """Generate expected result file path."""
    filename = f"spearman_{genome}_{family1}_vs_{family2}_{window_size}kb.json"
    return os.path.join(result_dir, filename)


def collect_data(result_dir, genomes, family_pairs, window_sizes):
    """Collect correlation data from result files."""
    data = {}

    for genome in genomes:
        for pair_str in family_pairs:
            f1, f2 = pair_str.split('_')

            series_key = f"{genome} (Family {f1} vs {f2})"
            series_correlations = []
            series_sample_sizes = []
            series_windows = []

            for ws in window_sizes:
                result_file = get_result_file_path(result_dir, genome, f1, f2, ws)
                result = load_result_file(result_file)

                if result is None:
                    print(f"Warning: Missing result file: {result_file}", file=sys.stderr)
                    continue

                corr = result['statistics']['spearman_correlation']
                n_windows = result['metadata']['n_windows']

                if n_windows == 0:
                    continue

                series_correlations.append(corr)
                series_sample_sizes.append(n_windows)
                series_windows.append(ws)

            if series_correlations:
                data[series_key] = {
                    'genome': genome,
                    'window_sizes': series_windows,
                    'correlations': series_correlations,
                    'sample_sizes': series_sample_sizes
                }

    return data


def plot_correlations(data, output_file, window_sizes, family_pairs):
    """Create a publication-quality plot of intra-genome correlations."""
    if not data:
        print("Error: No data to plot!", file=sys.stderr)
        sys.exit(1)

    # Generate colors from a colormap for flexibility
    base_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    color_idx = 0

    _, ax = plt.subplots(figsize=(11, 7), dpi=300)

    for series_key in sorted(data.keys()):
        series_data = data[series_key]
        ws = series_data['window_sizes']
        corr = series_data['correlations']

        color = base_colors[color_idx % len(base_colors)]
        color_idx += 1

        ax.plot(ws, corr, marker='o', markersize=8, linewidth=2.5,
                label=series_key, color=color, alpha=0.85)

    # Build title from family pairs
    pair_labels = [f"Family {p.replace('_', ' vs ')}" for p in family_pairs]
    title_pairs = ', '.join(pair_labels)

    ax.set_xlabel('Window Size (kb)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Spearman Correlation Coefficient (\u03c1)', fontsize=13, fontweight='bold')
    ax.set_title(f'Intra-genome {title_pairs} Correlation\nAcross Genomic Scales',
                fontsize=14, fontweight='bold', pad=15)

    ax.set_xticks(window_sizes)
    ax.set_xticklabels([str(ws) for ws in window_sizes])

    # Dynamic y-axis limits (filter NaN from constant-value correlations)
    import math
    all_corrs = [c for sd in data.values() for c in sd['correlations'] if not math.isnan(c)]
    if all_corrs:
        y_max = max(max(all_corrs) * 1.15, 0.1)
        y_min = min(min(all_corrs) * 1.15 if min(all_corrs) < 0 else -0.05, -0.05)
    else:
        y_max, y_min = 0.8, -0.05
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, max(window_sizes) + 5)

    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)

    ax.legend(loc='upper right', fontsize=12, framealpha=0.95, edgecolor='black')

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Figure saved to: {output_file}")


def main():
    """Main function."""
    args = parse_arguments()

    if not os.path.isdir(args.result_dir):
        print(f"Error: Result directory not found: {args.result_dir}", file=sys.stderr)
        sys.exit(1)

    print("Collecting intra-genome correlation data from result files...")
    data = collect_data(args.result_dir, args.genomes, args.family_pairs, args.window_sizes)

    if not data:
        print("Error: No valid data files found!", file=sys.stderr)
        sys.exit(1)

    print(f"Found data for {len(data)} series")
    for key in sorted(data.keys()):
        n_windows = len(data[key]['window_sizes'])
        print(f"  {key}: {n_windows} window sizes")

    print(f"\nGenerating plot...")
    plot_correlations(data, args.output, args.window_sizes, args.family_pairs)


if __name__ == "__main__":
    main()

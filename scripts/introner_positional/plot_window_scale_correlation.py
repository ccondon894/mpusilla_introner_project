#!/usr/bin/env python3
"""
Generate a publication-quality figure comparing Spearman correlation coefficients
across different genomic window sizes for different introner families.
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
        description='Plot window-scale correlations for introner families',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--result-dir', required=True,
                       help='Directory containing JSON result files')
    parser.add_argument('--families', nargs='+', default=['all', '1', '2', '3', '4', '14', '15'],
                       help='Families to include in plot (default: all 1 2 3 4 14 15)')
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


def get_result_file_path(result_dir, family, window_size):
    """Generate expected result file path for a family and window size."""
    family_label = family if family == 'all' else f'family_{family}'
    filename = f"spearman_{family_label}_{window_size}kb.json"
    return os.path.join(result_dir, filename)


def collect_data(result_dir, families, window_sizes):
    """Collect correlation data from result files for all families and window sizes."""
    data = {}

    for family in families:
        family_correlations = []
        family_sample_sizes = []
        family_windows = []

        for ws in window_sizes:
            result_file = get_result_file_path(result_dir, family, ws)
            result = load_result_file(result_file)

            if result is None:
                print(f"Warning: Missing result file: {result_file}", file=sys.stderr)
                continue

            corr = result['statistics']['spearman_correlation']
            n_windows = result['metadata']['n_windows']

            if n_windows == 0:
                continue

            family_correlations.append(corr)
            family_sample_sizes.append(n_windows)
            family_windows.append(ws)

        if family_correlations:
            data[family] = {
                'window_sizes': family_windows,
                'correlations': family_correlations,
                'sample_sizes': family_sample_sizes
            }

    return data


def plot_correlations(data, output_file, window_sizes):
    """Create a publication-quality plot of correlation coefficients across families."""
    import math

    if not data:
        print("Error: No data to plot!", file=sys.stderr)
        sys.exit(1)

    color_map = {
        'all': '#1f77b4',
        '1': '#ff7f0e',
        '2': '#2ca02c',
        '3': '#d62728',
        '4': '#9467bd',
        '14': '#8c564b',
        '15': '#e377c2'
    }

    _, ax = plt.subplots(figsize=(12, 7), dpi=300)

    for family, family_data in sorted(data.items()):
        if family in {'14', '15', '4'}:
            continue # skipping because not enough data points 
            
        ws = family_data['window_sizes']
        corr = family_data['correlations']

        color = color_map.get(family, '#000000')

        ax.plot(ws, corr, marker='o', markersize=7, linewidth=2.5,
                label=f'Family {family}' if family != 'all' else 'All families',
                color=color, alpha=0.8)

    ax.set_xlabel('Window Size (kb)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Spearman Correlation Coefficient (\u03c1)', fontsize=13, fontweight='bold')
    ax.set_title('Introner Density Correlation Across Genomic Scales\nby Introner Family',
                fontsize=14, fontweight='bold', pad=15)

    ax.set_xticks(window_sizes)
    ax.set_xticklabels([str(ws) for ws in window_sizes])

    # Dynamic y-axis limits based on data (filter NaN from constant-value correlations)
    all_corrs = [c for fd in data.values() for c in fd['correlations'] if not math.isnan(c)]
    if all_corrs:
        y_max = max(max(all_corrs) * 1.15, 0.1)
        y_min = min(min(all_corrs) * 1.15 if min(all_corrs) < 0 else -0.05, -0.05)
    else:
        y_max, y_min = 0.8, -0.05
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, max(window_sizes) + 5)

    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)

    ax.legend(loc='upper left', fontsize=11, framealpha=0.95, edgecolor='black', )

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

    print("Collecting correlation data from result files...")
    data = collect_data(args.result_dir, args.families, args.window_sizes)

    if not data:
        print("Error: No valid data files found!", file=sys.stderr)
        sys.exit(1)

    print(f"Found data for {len(data)} families")
    for family in sorted(data.keys()):
        n_windows = len(data[family]['window_sizes'])
        print(f"  Family {'(all)' if family == 'all' else family}: {n_windows} window sizes")

    print(f"\nGenerating plot with {len(data)} families...")
    plot_correlations(data, args.output, args.window_sizes)


if __name__ == "__main__":
    main()

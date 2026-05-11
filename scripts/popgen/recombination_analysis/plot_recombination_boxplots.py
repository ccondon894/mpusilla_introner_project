#!/usr/bin/env python3
"""
Publication-quality box plots comparing recombination rates across introner categories.

Creates a single-panel figure with 5 introner categories:
- Non-introner (baseline)
- All introners
- Polymorphic introners
- Recent loss
- Recent gain

Performs all pairwise Mann-Whitney U tests with Bonferroni correction
and displays significance bars with stars.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from itertools import combinations
import warnings
import argparse

warnings.filterwarnings('ignore')

# Set style for publication-quality plots
sns.set_style("white")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']


def load_recombination_data(base_dir='./', min_rate=1e-14, max_rate=5e-11):
    """
    Load recombination rate data from multiple TSV files.

    Filters to reasonable recombination rate range to exclude artifacts
    (extreme outliers likely from windows with no pyrho data).

    Parameters:
    -----------
    base_dir : str
        Directory containing input files
    min_rate : float
        Minimum recombination rate to include (filters out extreme outliers)
    max_rate : float
        Maximum recombination rate to include

    Returns:
    --------
    dict
        Dictionary of {category_name: rate_array}
    """
    data = {}

    # Load all-introner analysis
    all_file = f"{base_dir}/gene_exonic_introners_all_5kb_updated_windows.tsv"
    all_df = pd.read_csv(all_file, sep='\t')

    # Non-introner baseline
    non_introner = all_df[all_df['type'] == 'non_introner_containing']['rate'].dropna().values
    non_introner = non_introner[(non_introner >= min_rate) & (non_introner <= max_rate)]
    data['non_introner'] = non_introner

    # All-introner
    all_introner = all_df[all_df['type'] == 'introner_containing']['rate'].dropna().values
    all_introner = all_introner[(all_introner >= min_rate) & (all_introner <= max_rate)]
    data['all_introner'] = all_introner

    # Load polymorphic introner analysis
    poly_file = f"{base_dir}/gene_exonic_polymorphic_introners_5kb_updated_windows.tsv"
    poly_df = pd.read_csv(poly_file, sep='\t')
    polymorphic_introner = poly_df[poly_df['type'] == 'introner_containing']['rate'].dropna().values
    polymorphic_introner = polymorphic_introner[(polymorphic_introner >= min_rate) & (polymorphic_introner <= max_rate)]
    data['polymorphic_introner'] = polymorphic_introner

    # Load frequency-based allele classification
    freq_file = f"{base_dir}/gene_exonic_frequency_based_5kb_updated_windows.tsv"
    freq_df = pd.read_csv(freq_file, sep='\t')

    recent_gain = freq_df[freq_df['type'] == 'recent_gain_containing']['rate'].dropna().values
    recent_gain = recent_gain[(recent_gain >= min_rate) & (recent_gain <= max_rate)]
    data['recent_gain'] = recent_gain

    recent_loss = freq_df[freq_df['type'] == 'recent_loss_containing']['rate'].dropna().values
    recent_loss = recent_loss[(recent_loss >= min_rate) & (recent_loss <= max_rate)]
    data['recent_loss'] = recent_loss

    return data


def load_group2_data(base_dir='./', min_rate=1e-14, max_rate=5e-11):
    """
    Load Group 2 (RCC1749/RCC3052) recombination rate data from TSV file.

    Filters to reasonable recombination rate range to exclude artifacts.

    Parameters:
    -----------
    base_dir : str
        Directory containing input files
    min_rate : float
        Minimum recombination rate to include (filters out extreme outliers)
    max_rate : float
        Maximum recombination rate to include

    Returns:
    --------
    dict
        Dictionary of {category_name: rate_array}
    """
    data = {}

    # Load Group 2 introner analysis
    group2_file = f"{base_dir}/gene_exonic_group2_introners_5kb_updated_windows.tsv"
    group2_df = pd.read_csv(group2_file, sep='\t')

    # Non-introner from Group 2
    non_introner = group2_df[group2_df['type'] == 'non_introner_containing']['rate'].dropna().values
    non_introner = non_introner[(non_introner >= min_rate) & (non_introner <= max_rate)]
    data['non_introner'] = non_introner

    # All-introner from Group 2
    all_introner = group2_df[group2_df['type'] == 'introner_containing']['rate'].dropna().values
    all_introner = all_introner[(all_introner >= min_rate) & (all_introner <= max_rate)]
    data['all_introner'] = all_introner

    return data


def perform_pairwise_tests(data_dict, category_order, alpha=0.05):
    """
    Perform all pairwise Mann-Whitney U tests with Bonferroni correction.

    Parameters:
    -----------
    data_dict : dict
        Dictionary of {category: rate_array}
    category_order : list
        Order of categories
    alpha : float
        Significance level before correction

    Returns:
    --------
    list
        List of tuples (cat1_idx, cat2_idx, U_statistic, p_value, is_significant)
    """
    results = []
    n_comparisons = len(category_order) * (len(category_order) - 1) // 2
    bonferroni_alpha = alpha / n_comparisons

    print(f"\n{'='*80}")
    print(f"Mann-Whitney U Tests (Bonferroni-corrected α = {bonferroni_alpha:.4f})")
    print(f"{'='*80}\n")

    # Get all pairwise combinations
    for i, cat1 in enumerate(category_order):
        for j, cat2 in enumerate(category_order):
            if i < j:  # Only test each pair once
                data1 = data_dict[cat1]
                data2 = data_dict[cat2]

                # Perform Mann-Whitney U test (two-sided)
                u_stat, p_value = mannwhitneyu(data1, data2, alternative='two-sided')

                # Check if significant after Bonferroni correction
                is_significant = p_value < bonferroni_alpha

                # Format significance stars
                if p_value < 0.001:
                    stars = '***'
                elif p_value < 0.01:
                    stars = '**'
                elif p_value < 0.05:
                    stars = '*'
                else:
                    stars = 'ns'

                results.append((i, j, u_stat, p_value, is_significant))

                # Print results
                label1 = cat1.replace('_', ' ').title()
                label2 = cat2.replace('_', ' ').title()
                sig_marker = '✓' if is_significant else ' '
                print(f"{sig_marker} {label1:25s} vs {label2:25s}: "
                      f"U={u_stat:10.1f}, p={p_value:.4e} {stars:>3s}")

    return results


def get_significance_stars(p_value):
    """Convert p-value to star annotation."""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'


def add_significance_bar(ax, x1, x2, y, stars, height_factor=1.0):
    """
    Draw a bracket with significance annotation.

    Parameters:
    -----------
    ax : matplotlib axis
        Axis to draw on
    x1, x2 : int
        X positions of the two categories being compared
    y : float
        Y position for the bar (in data coordinates)
    stars : str
        Significance annotation ('*', '**', '***', 'ns')
    height_factor : float
        Multiplier for bar height positioning
    """
    # Draw horizontal line
    ax.plot([x1, x2], [y, y], 'k-', linewidth=1.5, clip_on=False)

    # Draw vertical ticks at ends
    tick_height = y * 0.05  # 5% of y position
    ax.plot([x1, x1], [y, y - tick_height], 'k-', linewidth=1.5, clip_on=False)
    ax.plot([x2, x2], [y, y - tick_height], 'k-', linewidth=1.5, clip_on=False)

    # Add text annotation
    text_y = y * 1.1  # Position text slightly above bar
    ax.text((x1 + x2) / 2, text_y, stars, ha='center', va='bottom',
            fontsize=12, fontweight='bold', clip_on=False)


def create_boxplot_with_significance_on_axis(ax, data_dict, category_order, colors,
                                              test_results, labels=None,
                                              allowed_comparisons=None,
                                              title=None, alpha=0.05,
                                              show_strip=False):
    """
    Create box plot with significance bars on a given axis.

    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on
    data_dict : dict
        Dictionary of {category: rate_array}
    category_order : list
        Order to display categories
    colors : dict
        Mapping of category to color
    test_results : list
        Results from pairwise tests
    labels : dict, optional
        Mapping of category names to display labels. If None, use category names.
    allowed_comparisons : set, optional
        Set of (i, j) tuples for which comparisons to display. If None, show all.
    title : str, optional
        Panel title
    alpha : float
        Significance threshold for displaying bars
    show_strip : bool, optional
        If True, overlay strip plot under box plot
    """
    # Use custom labels if provided
    if labels is None:
        labels = {cat: cat.replace('_', ' ').title() for cat in category_order}

    # Prepare data for box plot
    plot_data = [data_dict[cat] for cat in category_order]

    # Add strip plot if requested (draw first so it appears behind box plot)
    if show_strip:
        for i, cat in enumerate(category_order):
            rates = data_dict[cat]
            # Subsample if too many points (for performance and visibility)
            if len(rates) > 500:
                rates_sample = np.random.choice(rates, 500, replace=False)
            else:
                rates_sample = rates

            # Add jitter to x positions for visibility
            x_jitter = np.random.normal(i, 0.08, size=len(rates_sample))
            ax.scatter(x_jitter, rates_sample,
                      color=colors[cat], alpha=0.3, s=20,
                      edgecolors='black', zorder=1)

    # Create box plot (on top of strip plot)
    bp = ax.boxplot(plot_data, positions=range(len(category_order)),
                    widths=0.6, patch_artist=True,
                    showfliers=not show_strip,  # Hide outliers if showing strip plot
                    flierprops={'marker': 'o', 'markersize': 3, 'alpha': 0.3,
                               'markerfacecolor': 'gray', 'markeredgecolor': 'none'},
                    medianprops={'color': 'black', 'linewidth': 2},
                    boxprops={'linewidth': 1.5},
                    whiskerprops={'linewidth': 1.5},
                    capprops={'linewidth': 1.5},
                    zorder=2)

    # Color the boxes
    for i, (patch, cat) in enumerate(zip(bp['boxes'], category_order)):
        patch.set_facecolor(colors[cat])
        patch.set_alpha(0.7 if not show_strip else 0.8)

    # Set labels
    ax.set_xticks(range(len(category_order)))
    ax.set_xticklabels([labels[cat] for cat in category_order],
                       fontsize=14, rotation=45, ha='right')
    ax.set_ylabel('Recombination Rate (r)', fontsize=14)
    ax.set_yscale('log')
    ax.tick_params(axis='y', labelsize=13)
    # ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)

    # Set y-axis limits
    all_rates = []
    for cat in category_order:
        all_rates.extend(data_dict[cat])
    y_min = np.min(all_rates) * 0.5
    y_max = np.max(all_rates) * 3  # Initial space for significance bars

    # Filter results to only show whitelisted comparisons
    if allowed_comparisons is not None:
        sig_results = [(i, j, p) for i, j, _, p, _ in test_results
                       if (i, j) in allowed_comparisons]
    else:
        # Show all significant results
        sig_results = [(i, j, p) for i, j, _, p, sig in test_results if sig]

    if sig_results:
        # Sort by span (distance), then by starting position for better visual organization
        sig_results_sorted = sorted(sig_results, key=lambda x: (abs(x[1] - x[0]), x[0]))

        # Assign each bar a unique layer with increased spacing
        base_y = y_max
        max_bar_y = base_y
        for idx, (i, j, p_val) in enumerate(sig_results_sorted):
            # Each bar gets its own layer with exponentially increasing height
            # Use larger multiplier (1.5) for better separation
            bar_y = base_y * (1.5 ** idx)

            stars = get_significance_stars(p_val)
            add_significance_bar(ax, i, j, bar_y, stars)

            # Track the maximum bar height
            max_bar_y = max(max_bar_y, bar_y)

        # Set y_max to accommodate all bars plus space for text annotations
        y_max = max_bar_y * 2.5

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(-0.5, len(category_order) - 0.5)


def main():
    """Create publication-quality 2-panel box plot comparison figure."""

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generate recombination rate box plot comparisons',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Box plots only (default)
  python plot_recombination_boxplots.py

  # Box plots with strip plots overlaid
  python plot_recombination_boxplots.py --strip

  # Specify custom output filename
  python plot_recombination_boxplots.py -o my_figure

  # Combine options
  python plot_recombination_boxplots.py --strip -o figure_with_strips
        """)
    parser.add_argument('--strip', action='store_true',
                       help='Overlay strip plots under box plots to show individual data points')
    parser.add_argument('-o', '--output', type=str, default='recombination_boxplots_comparison',
                       help='Output filename prefix (without extension). Default: recombination_boxplots_comparison')
    parser.add_argument('--base_dir', type=str, default='./',
                       help='Directory containing the gene_exonic_*_windows.tsv input files (default: ./)')
    args = parser.parse_args()

    base_dir = args.base_dir

    # Load data
    print("📊 Loading Group 1 data...")
    group1_data = load_recombination_data(base_dir)

    print("📊 Loading Group 2 data...")
    group2_data = load_group2_data(base_dir)

    # Define Group 1 categories and colors
    group1_categories = ['non_introner', 'all_introner', 'polymorphic_introner',
                         'recent_gain', 'recent_loss']
    group1_colors = {
        'non_introner': '#888888',           # grey
        'all_introner': '#2ecc71',           # green
        'polymorphic_introner': '#1abc9c',   # teal
        'recent_gain': '#3498db',            # blue
        'recent_loss': '#8B5AA3'             # purple (viridis-inspired)
    }

    # Customizable labels for Group 1
    group1_labels = {
        'non_introner': 'Non-Introner\nGroup 1',
        'all_introner': 'All Introner\nGroup 1',
        'polymorphic_introner': 'Polymorphic\nIntroner',
        'recent_gain': 'Recent Gain',
        'recent_loss': 'Recent Loss'
    }

    # Define Group 2 categories and colors
    group2_categories = ['non_introner_g2', 'all_introner_g2']
    group2_colors = {
        'non_introner_g2': '#888888',  # grey
        'all_introner_g2': '#e74c3c'   # red
    }

    # Customizable labels for Group 2
    group2_labels = {
        'non_introner_g2': 'Non-Introner\nGroup 2',
        'all_introner_g2': 'All Introner\nGroup 2'
    }

    # Combine all data into single dictionary for unified plot
    combined_data = {}
    for cat in group1_categories:
        combined_data[cat] = group1_data[cat]
    # Add Group 2 data with distinct keys
    combined_data['non_introner_g2'] = group2_data['non_introner']
    combined_data['all_introner_g2'] = group2_data['all_introner']

    # Combine categories in display order
    all_categories = group1_categories + group2_categories

    # Combine colors
    all_colors = {**group1_colors, **group2_colors}

    # Combine labels
    all_labels = {**group1_labels, **group2_labels}

    # Statistical tests for Group 1
    print("\n" + "="*80)
    print("GROUP 1 STATISTICAL TESTS (CCMP1545)")
    print("="*80)
    group1_tests = perform_pairwise_tests(group1_data, group1_categories)

    print("\n" + "="*80)
    print("GROUP 2 STATISTICAL TESTS (RCC1749/RCC3052)")
    print("="*80)
    group2_tests = perform_pairwise_tests(group2_data, ['non_introner', 'all_introner'])

    # Print summary statistics
    print(f"\n{'='*80}")
    print("GROUP 1 Summary Statistics")
    print(f"{'='*80}\n")
    for cat in group1_categories:
        rates = group1_data[cat]
        print(f"{cat.replace('_', ' ').title():25s}: "
              f"n={len(rates):4d}, median={np.median(rates):.2e}, "
              f"mean={np.mean(rates):.2e}")

    print(f"\n{'='*80}")
    print("GROUP 2 Summary Statistics")
    print(f"{'='*80}\n")
    # Use original keys from group2_data, but display with new labels
    for original_key, display_key in zip(['non_introner', 'all_introner'], group2_categories):
        rates = group2_data[original_key]
        print(f"{display_key.replace('_', ' ').title():25s}: "
              f"n={len(rates):4d}, median={np.median(rates):.2e}, "
              f"mean={np.mean(rates):.2e}")

    # Define which comparisons to show in the combined plot
    # all_categories = ['non_introner', 'all_introner', 'polymorphic_introner', 'recent_gain', 'recent_loss', 'non_introner_g2', 'all_introner_g2']
    # Indices:         [0,              1,              2,                      3,             4,              5,                 6]
    allowed_comparisons = {
        (0, 1),  # Group 1: Non-introner vs All introner
        (0, 2),  # Group 1: Non-introner vs Polymorphic introner
        (1, 2),  # Group 1: All introner vs Polymorphic introner
        (3, 4),  # Group 1: Recent gain vs Recent loss
        (5, 6),  # Group 2: Non-introner vs All introner
    }

    # Perform statistical tests on combined data for display
    print("\n" + "="*80)
    print("COMBINED STATISTICAL TESTS (for visualization)")
    print("="*80)
    combined_tests = perform_pairwise_tests(combined_data, all_categories)

    # Create single panel figure
    plot_type = "box plots + strip plots" if args.strip else "box plots"
    print(f"\n🎨 Creating single-panel figure ({plot_type})...")
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    # Create combined plot
    print("   Plotting all 7 categories with 5 significance bars")
    create_boxplot_with_significance_on_axis(
        ax, combined_data, all_categories, all_colors,
        combined_tests, labels=all_labels,
        allowed_comparisons=allowed_comparisons,
        show_strip=args.strip
    )

    plt.tight_layout()

    # Save outputs
    print("\n💾 Saving outputs...")
    pdf_file = f'{args.output}.pdf'
    png_file = f'{args.output}.png'

    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ Saved: {pdf_file}")

    plt.savefig(png_file, dpi=300, bbox_inches='tight', format='png')
    print(f"   ✓ Saved: {png_file}")

    print("\n✅ Figure generation complete!")
    print(f"\nFigure: Single-panel {plot_type} (14 × 8 inches, 300 DPI)")
    print(f"Categories: 5 Group 1 + 2 Group 2 = 7 total categories")
    print(f"Statistical tests: All pairwise Mann-Whitney U tests reported in console")
    print(f"Significance bars: 5 selected comparisons displayed with stars")
    if args.strip:
        print(f"Strip plots: Individual data points shown (max 500 per category)")
    print(f"\nUsage:")
    print(f"  Box plots only:         python plot_recombination_boxplots.py")
    print(f"  Box + strip plots:      python plot_recombination_boxplots.py --strip")
    print(f"  Custom filename:        python plot_recombination_boxplots.py -o my_figure")
    print(f"  Combined options:       python plot_recombination_boxplots.py --strip -o my_figure")
    print(f"\nTo customize labels: Edit the group1_labels and group2_labels dictionaries in main()")


if __name__ == '__main__':
    main()

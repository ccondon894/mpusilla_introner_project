#!/usr/bin/env python3
"""
Publication-quality strip plots comparing recombination rates across introner categories.

Creates a 2-panel figure showing individual data points:
- Panel A: Group 1 (CCMP1545) with 5 introner categories
- Panel B: Group 2 (RCC1749/RCC3052) with 2 introner categories

Each point represents a single genomic window with its measured recombination rate.
Median markers (black diamonds) show the central tendency for each category.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# Set style for publication-quality plots
sns.set_style("whitegrid", {'grid.linewidth': 0.5, 'grid.color': 'lightgrey'})
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']


def load_group1_data(base_dir='./', min_rate=1e-14, max_rate=5e-11):
    """
    Load Group 1 recombination rate data from multiple TSV files.

    Filters to reasonable recombination rate range to exclude artifacts
    (extreme outliers likely from windows with no pyrho data).

    Parameters:
    -----------
    min_rate : float
        Minimum recombination rate to include (filters out extreme outliers)
    max_rate : float
        Maximum recombination rate to include
    """
    data = {}

    # Load all-introner analysis
    all_file = f"{base_dir}/gene_exonic_introners_all_5kb_updated_windows.tsv"
    all_df = pd.read_csv(all_file, sep='\t')

    # Non-introner from all-introner analysis
    non_introner_group1 = all_df[all_df['type'] == 'non_introner_containing']['rate'].dropna().values
    non_introner_group1 = non_introner_group1[(non_introner_group1 >= min_rate) & (non_introner_group1 <= max_rate)]
    data['non_introner'] = non_introner_group1

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
    Load Group 2 recombination rate data from TSV file.

    Filters to reasonable recombination rate range to exclude artifacts.

    Parameters:
    -----------
    min_rate : float
        Minimum recombination rate to include (filters out extreme outliers)
    max_rate : float
        Maximum recombination rate to include
    """
    data = {}

    # Load Group 2 introner analysis
    group2_file = f"{base_dir}/gene_exonic_group2_introners_5kb_updated_windows.tsv"
    group2_df = pd.read_csv(group2_file, sep='\t')

    # Non-introner from Group 2
    non_introner_group2 = group2_df[group2_df['type'] == 'non_introner_containing']['rate'].dropna().values
    non_introner_group2 = non_introner_group2[(non_introner_group2 >= min_rate) & (non_introner_group2 <= max_rate)]
    data['non_introner'] = non_introner_group2

    # All-introner from Group 2
    all_introner_group2 = group2_df[group2_df['type'] == 'introner_containing']['rate'].dropna().values
    all_introner_group2 = all_introner_group2[(all_introner_group2 >= min_rate) & (all_introner_group2 <= max_rate)]
    data['all_introner'] = all_introner_group2

    return data


def create_strip_plot_panel(ax, data_dict, category_order, colors, title, fontsize_title=14, fontsize_label=13, fontsize_tick=12):
    """
    Create a strip plot panel showing individual data points with median markers.

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
    title : str
        Panel title
    fontsize_* : int
        Font sizes for different elements
    """

    # Prepare data for plotting
    plot_data = []
    for category in category_order:
        if category in data_dict:
            rates = data_dict[category]
            for rate in rates:
                plot_data.append({'category': category, 'rate': rate})

    plot_df = pd.DataFrame(plot_data)

    # Create strip plot with jitter
    for i, category in enumerate(category_order):
        if category in plot_df['category'].values:
            cat_data = plot_df[plot_df['category'] == category]['rate'].values
            # Add jitter to x position (increased from 0.04 to 0.09 for more spread)
            x_pos = np.random.normal(i, 0.09, size=len(cat_data))
            ax.scatter(x_pos, cat_data, alpha=0.5, s=25, color=colors[category],
                      edgecolor='black', linewidth=0.5, rasterized=True)

    # Add median markers (larger diamonds)
    for i, category in enumerate(category_order):
        if category in data_dict:
            rates = data_dict[category]
            median_rate = np.median(rates)
            ax.scatter(i, median_rate, marker='D', s=200, color='#111111',
                      edgecolor='white', linewidth=1.5, zorder=10)

    # Set labels and formatting
    ax.set_xticks(range(len(category_order)))
    ax.set_xticklabels([cat.replace('_', ' ').title() for cat in category_order],
                       fontsize=fontsize_tick, rotation=45, ha='right')
    ax.set_ylabel('Recombination Rate (r)', fontsize=fontsize_label, fontweight='bold')
    ax.set_yscale('log')
    ax.set_title(title, fontsize=fontsize_title, fontweight='bold', pad=15)

    # Improve tick label formatting
    ax.tick_params(axis='y', labelsize=fontsize_tick)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)

    # Set y-axis limits and x-axis limits
    all_rates = []
    for cat in category_order:
        if cat in data_dict:
            all_rates.extend(data_dict[cat])
    if all_rates:
        y_min = np.min(all_rates) * 0.5
        y_max = np.max(all_rates) * 2
        ax.set_ylim(y_min, y_max)

    # Set x-axis limits with padding
    ax.set_xlim(-0.5, len(category_order) - 0.5)


def main():
    """Create publication-quality 2-panel comparison figure."""

    base_dir = './'

    print("📊 Loading Group 1 data...")
    group1_data = load_group1_data(base_dir)

    print("📊 Loading Group 2 data...")
    group2_data = load_group2_data(base_dir)

    # Define category orders and colors
    group1_categories = ['non_introner', 'all_introner', 'polymorphic_introner', 'recent_gain', 'recent_loss']
    group1_colors = {
        'non_introner': '#888888',           # grey
        'all_introner': '#2ecc71',            # green
        'polymorphic_introner': '#1abc9c',   # teal
        'recent_gain': '#3498db',             # blue
        'recent_loss': '#e74c3c'              # orange
    }

    group2_categories = ['non_introner', 'all_introner']
    group2_colors = {
        'non_introner': '#888888',  # grey
        'all_introner': '#e74c3c'   # red
    }

    # Create figure with 2 panels
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    print("🎨 Creating Panel A (Group 1)...")
    create_strip_plot_panel(
        axes[0],
        group1_data,
        group1_categories,
        group1_colors,
        title='Panel A: Group 1 (CCMP1545)',
        fontsize_title=14,
        fontsize_label=13,
        fontsize_tick=12
    )

    print("🎨 Creating Panel B (Group 2)...")
    create_strip_plot_panel(
        axes[1],
        group2_data,
        group2_categories,
        group2_colors,
        title='Panel B: Group 2 (RCC1749/RCC3052)',
        fontsize_title=14,
        fontsize_label=13,
        fontsize_tick=12
    )

    # Remove y-axis label from Panel B (already shared/labeled in A)
    axes[1].set_ylabel('')

    # Adjust layout
    plt.tight_layout()

    # Save as high-quality PDF and PNG
    print("💾 Saving outputs...")
    pdf_file = 'recombination_rate_comparison.pdf'
    png_file = 'recombination_rate_comparison.png'

    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', format='pdf')
    print(f"   ✓ Saved: {pdf_file}")

    plt.savefig(png_file, dpi=300, bbox_inches='tight', format='png')
    print(f"   ✓ Saved: {png_file}")

    print("\n✅ Figure generation complete!")
    print(f"\nFigure dimensions: 12 × 5 inches (300 DPI)")
    print(f"Panel A: Strip plot with 5 Group 1 introner categories")
    print(f"Panel B: Strip plot with 2 Group 2 introner categories")
    print(f"\nVisualization: Individual data points (colored by category)")
    print(f"               Median markers shown as black diamonds")

    # Print summary statistics
    print("\n📈 Summary Statistics:")
    print("\nGROUP 1:")
    for cat in group1_categories:
        if cat in group1_data:
            rates = group1_data[cat]
            print(f"  {cat.replace('_', ' ').title():25s}: n={len(rates):4d}, median={np.median(rates):.2e}")

    print("\nGROUP 2:")
    for cat in group2_categories:
        if cat in group2_data:
            rates = group2_data[cat]
            print(f"  {cat.replace('_', ' ').title():25s}: n={len(rates):4d}, median={np.median(rates):.2e}")


if __name__ == '__main__':
    main()

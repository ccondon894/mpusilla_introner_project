#!/usr/bin/env python3
"""
Phase 2 Results Visualization Script
Introner Functional Enrichment Analysis

Generates visualizations for Phase 2 results:
- Enrichment heatmap across GO categories and introner types
- Volcano plots showing effect size vs significance
- Top enriched categories bar chart

Usage: python plot_phase2_results.py <phase2_results.tsv> <output_dir>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# Set consistent style
plt.style.use('default')
sns.set_palette("husl")

def load_phase2_data(phase2_file):
    """Load and process Phase 2 enrichment results."""
    print(f"Loading Phase 2 data from {phase2_file}...")
    
    df = pd.read_csv(phase2_file, sep='\t')
    
    # Clean up data
    df = df[df['genes_with_introner_and_go'] > 0]  # Remove zero-overlap results
    
    print(f"  Loaded {len(df)} enrichment tests")
    print(f"  {df['significant'].sum()} significant results (FDR < 0.05)")
    print(f"  {len(df['go_category'].unique())} unique GO categories")
    print(f"  {len(df['introner_type'].unique())} introner types")
    
    return df

def plot_enrichment_heatmap(df, output_dir):
    """
    Plot 1: Enrichment Heatmap
    Heatmap showing FDR values (colors) and enrichment/depletion proportions (text)
    across GO categories (rows) and introner types (columns).
    """
    print("Creating enrichment heatmap...")

    # Get GO categories that have at least one significant result (FDR < 0.05)
    significant_go_categories = df[df['fdr'] < 0.05]['go_category'].unique()

    if len(significant_go_categories) == 0:
        print("  No significant enrichments to plot")
        return

    # Filter to only those GO categories, but show ALL FDR values for them
    filtered_df = df[df['go_category'].isin(significant_go_categories)]

    # Get all unique introner types to ensure complete coverage
    all_introner_types = df['introner_type'].unique()

    # Create pivot tables for FDR values and fold enrichment
    fdr_data = filtered_df.pivot_table(
        index='go_category',
        columns='introner_type',
        values='fdr'
    )

    # Use fold_enrichment for colors (biological effect size)
    fold_enrichment_data = filtered_df.pivot_table(
        index='go_category',
        columns='introner_type',
        values='fold_enrichment'
    )

    # Reindex to include all introner types as columns
    fdr_data = fdr_data.reindex(columns=all_introner_types)
    fold_enrichment_data = fold_enrichment_data.reindex(columns=all_introner_types)

    # Determine color scale range: 0 to max fold enrichment
    max_fold_enrichment = filtered_df['fold_enrichment'].max()
    print(f"  Color scale: 0 to {max_fold_enrichment:.2f} (fold enrichment)")

    # Fill missing fold enrichment values with 0 (gray color)
    fold_enrichment_filled = fold_enrichment_data.fillna(0)

    if len(fold_enrichment_filled) == 0:
        print("  No significant enrichments to plot")
        return

    # Create custom annotation matrix with FDR values or significance stars
    annot_data = fdr_data.copy()
    for i in range(len(annot_data)):
        for j in range(len(annot_data.columns)):
            fdr_val = fdr_data.iloc[i, j]
            if pd.isna(fdr_val):
                annot_data.iloc[i, j] = "NA"
            else:
                # Convert FDR to significance stars
                if fdr_val < 0.001:
                    annot_data.iloc[i, j] = "***"
                elif fdr_val < 0.01:
                    annot_data.iloc[i, j] = "**"
                elif fdr_val < 0.05:
                    annot_data.iloc[i, j] = "*"
                else:
                    annot_data.iloc[i, j] = f"{fdr_val:.2f}"

    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(8, len(fold_enrichment_filled) * 0.4)))

    # Create custom colormap: gray for 0, then viridis for positive values
    from matplotlib.colors import ListedColormap
    import matplotlib.cm as cm

    # Create colormap that handles 0 values as gray
    viridis = cm.get_cmap('RdBu')
    colors = ['#808080']  # Gray for 0/NA values
    colors.extend([viridis(i/256) for i in range(256)])
    custom_cmap = ListedColormap(colors)

    # Heatmap with fold enrichment colors and FDR annotations
    sns.heatmap(fold_enrichment_filled,
                annot=annot_data,
                fmt='',  # Use string format since we're providing custom annotations
                cmap=custom_cmap,
                cbar_kws={'label': 'Fold Enrichment'},
                linewidths=0.5,
                vmin=0, vmax=max_fold_enrichment,
                ax=ax)

    ax.set_xlabel('Introner-Gene Category', fontsize=12, fontweight='bold')
    ax.set_ylabel('GO Category', fontsize=12, fontweight='bold')

    # Add legend for significance stars - positioned outside the plot area
    legend_text = "Significance: *** FDR<0.001, ** FDR<0.01, * FDR<0.05"
    ax.text(-0.15, 1.01, legend_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=0, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase2_enrichment_heatmap.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    return fold_enrichment_filled

def plot_volcano_plots(df, output_dir):
    """
    Plot 2: Volcano Plots
    5-panel subplot showing effect size vs significance for each introner type.
    """
    print("Creating volcano plots...")
    
    introner_types = df['introner_type'].unique()
    
    # Create 2x3 grid to accommodate 5 plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    colors = {'Enriched': '#d62728', 'Depleted': '#1f77b4', 'Not Significant': '#7f7f7f'}
    
    for i, introner_type in enumerate(introner_types):
        if i >= 6:  # Only plot first 6 types (we have 5)
            break
            
        ax = axes[i]
        subset = df[df['introner_type'] == introner_type].copy()
        
        if len(subset) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{introner_type}\n(No data)', fontsize=12, fontweight='bold')
            continue
        
        # Calculate log2 odds ratio and -log10 p-value
        subset['log2_or'] = np.log2(np.clip(subset['odds_ratio'], 0.01, 100))
        subset['neg_log10_p'] = -np.log10(np.clip(subset['p_value'], 1e-300, 1))
        
        # Create scatter plot with colors by enrichment type
        for enrichment_type, color in colors.items():
            mask = subset['enrichment_type'] == enrichment_type
            if mask.any():
                ax.scatter(subset[mask]['log2_or'], 
                          subset[mask]['neg_log10_p'],
                          c=color, alpha=0.7, s=50, 
                          label=f'{enrichment_type} ({mask.sum()})')
        
        # Add significance threshold line
        sig_threshold = -np.log10(0.05)
        ax.axhline(y=sig_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
        
        # Labels and title
        ax.set_xlabel('log2(Odds Ratio)', fontsize=11, fontweight='bold')
        ax.set_ylabel('-log10(p-value)', fontsize=11, fontweight='bold')
        ax.set_title(f'{introner_type}\n({len(subset)} tests)', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
        
        # Annotate most significant points
        if subset['significant'].any():
            sig_subset = subset[subset['significant']].nlargest(3, 'neg_log10_p')
            for _, row in sig_subset.iterrows():
                ax.annotate(row['go_category'][:15] + '...', 
                           (row['log2_or'], row['neg_log10_p']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8)
    
    # Hide unused subplots
    for i in range(len(introner_types), 6):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase2_volcano_plots.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_enrichments(df, output_dir):
    """
    Plot 3: Top Enriched Categories
    Horizontal bar chart showing the most significant enrichments across all groups.
    """
    print("Creating top enrichments plot...")
    
    # Get significant results and sort by FDR
    significant = df[df['significant']].copy()
    
    if len(significant) == 0:
        print("  No significant enrichments to plot")
        # Create empty plot with message
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, 'No significant enrichments found\n(FDR < 0.05)', 
                ha='center', va='center', transform=ax.transAxes, 
                fontsize=16, fontweight='bold')
        ax.set_title('Top Enriched GO Categories', fontsize=14, fontweight='bold')
        plt.savefig(os.path.join(output_dir, 'phase2_top_enriched.pdf'), 
                    dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Sort by FDR and take top 15
    top_results = significant.nsmallest(15, 'fdr')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(8, len(top_results) * 0.4)))
    
    # Color by introner type
    type_colors = {
        'Group1-Fixed': '#2E86AB',
        'Group1-Polymorphic': '#A23B72', 
        'Group2-Fixed': '#F18F01',
        'Group2-Polymorphic': '#C73E1D',
        'All-Introners': '#7B68EE'  # Medium slate blue for the new category
    }
    
    colors = [type_colors.get(itype, '#7f7f7f') for itype in top_results['introner_type']]
    
    # Create horizontal bar plot
    y_pos = np.arange(len(top_results))
    bars = ax.barh(y_pos, -np.log10(top_results['fdr']), color=colors, alpha=0.8, edgecolor='black')
    
    # Customize plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row['go_category']}\n({row['introner_type']})" 
                       for _, row in top_results.iterrows()], fontsize=10)
    ax.set_xlabel('-log10(FDR)', fontsize=12, fontweight='bold')
    ax.set_title('Top Enriched GO Categories', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add FDR values as text
    for i, (bar, fdr) in enumerate(zip(bars, top_results['fdr'])):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                f'FDR={fdr:.2e}', ha='left', va='center', fontsize=8)
    
    # Add legend for introner types
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.8, edgecolor='black') 
                      for color in type_colors.values()]
    ax.legend(legend_elements, type_colors.keys(), 
              loc='lower right', fontsize=10, title='Introner Type')
    
    # Invert y-axis to show most significant at top
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase2_top_enriched.pdf'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return top_results

def generate_summary_stats(df, output_dir):
    """Generate summary statistics text file."""
    print("Generating summary statistics...")
    
    stats_file = os.path.join(output_dir, 'phase2_summary_stats.txt')
    
    with open(stats_file, 'w') as f:
        f.write("=== PHASE 2 ENRICHMENT ANALYSIS SUMMARY ===\n\n")
        
        f.write(f"Total enrichment tests: {len(df)}\n")
        f.write(f"Significant results (FDR < 0.05): {df['significant'].sum()}\n")
        f.write(f"GO categories tested: {len(df['go_category'].unique())}\n")
        f.write(f"Introner types analyzed: {len(df['introner_type'].unique())}\n\n")
        
        # Results by introner type
        f.write("Results by introner type:\n")
        for introner_type in df['introner_type'].unique():
            subset = df[df['introner_type'] == introner_type]
            sig_count = subset['significant'].sum()
            enriched = (subset['enrichment_type'] == 'Enriched').sum()
            depleted = (subset['enrichment_type'] == 'Depleted').sum()
            
            f.write(f"  {introner_type}:\n")
            f.write(f"    Total tests: {len(subset)}\n")
            f.write(f"    Significant: {sig_count}\n")
            f.write(f"    Enriched: {enriched}\n")
            f.write(f"    Depleted: {depleted}\n\n")
        
        # Top significant results
        if df['significant'].any():
            f.write("Top 5 most significant enrichments:\n")
            top5 = df[df['significant']].nsmallest(5, 'fdr')
            for i, (_, row) in enumerate(top5.iterrows(), 1):
                f.write(f"  {i}. {row['go_category']} ({row['introner_type']})\n")
                f.write(f"     OR={row['odds_ratio']:.2f}, FDR={row['fdr']:.2e}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_phase2_results.py <phase2_results.tsv> <output_dir>")
        print("\nExample:")
        print("python plot_phase2_results.py results/phase2_enrichment_results.tsv results/plots/")
        sys.exit(1)
    
    phase2_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== PHASE 2 RESULTS VISUALIZATION ===")
    print(f"Input file: {phase2_file}")
    print(f"Output directory: {output_dir}")
    print()
    
    try:
        # Load data
        df = load_phase2_data(phase2_file)
        
        # Generate plots
        heatmap_data = plot_enrichment_heatmap(df, output_dir)
        plot_volcano_plots(df, output_dir)
        top_results = plot_top_enrichments(df, output_dir)
        
        # Generate summary statistics
        generate_summary_stats(df, output_dir)
        
        print("\n=== VISUALIZATION SUMMARY ===")
        print(f"Generated 3 plots in {output_dir}")
        print(f"  • phase2_enrichment_heatmap.pdf")
        print(f"  • phase2_volcano_plots.pdf")
        print(f"  • phase2_top_enriched.pdf")
        print(f"  • phase2_summary_stats.txt")
        
        print(f"\nKey findings:")
        print(f"  • {df['significant'].sum()}/{len(df)} tests were significant (FDR < 0.05)")
        if df['significant'].any():
            enriched_count = (df['enrichment_type'] == 'Enriched').sum()
            depleted_count = (df['enrichment_type'] == 'Depleted').sum()
            print(f"  • {enriched_count} enriched, {depleted_count} depleted categories")
            
            if heatmap_data is not None and len(heatmap_data) > 0:
                print(f"  • {len(heatmap_data)} GO categories with significant enrichments")
        
        print("\n✓ Phase 2 visualization completed successfully!")
        
    except Exception as e:
        print(f"✗ Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
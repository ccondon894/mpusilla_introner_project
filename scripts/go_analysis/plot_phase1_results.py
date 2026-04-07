#!/usr/bin/env python3
"""
Phase 1 Results Visualization Script
Introner Insertion Pattern Analysis

Generates visualizations for Phase 1 results:
- Introner distribution overview 
- Group comparison and overlap analysis
- GO annotation coverage assessment

Usage: python plot_phase1_results.py <phase1_results.tsv> <output_dir>
"""

import argparse
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

def load_phase1_data(phase1_file):
    """Load and process Phase 1 results."""
    print(f"Loading Phase 1 data from {phase1_file}...")
    
    df = pd.read_csv(phase1_file, sep='\t')
    
    # Parse GO terms back from string format
    df['go_terms'] = df['go_terms_str'].apply(
        lambda x: x.split(';') if x and str(x).strip() and str(x) != 'nan' else []
    )
    df['n_go_terms'] = df['go_terms'].apply(len)
    df['has_go_terms'] = df['n_go_terms'] > 0
    
    print(f"  Loaded {len(df)} genes")
    print(f"  {df['has_go_terms'].sum()} genes have GO annotations")
    
    return df

def plot_introner_distribution(df, output_dir):
    """
    Plot 1: Introner Distribution Overview
    Stacked bar chart + pie chart showing introner prevalence and group distributions.
    """
    print("Creating introner distribution plot...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Classify genes by introner presence
    df['has_group1'] = (df['group1_fixed'] > 0) | (df['group1_polymorphic'] > 0)
    df['has_group2'] = (df['group2_fixed'] > 0) | (df['group2_polymorphic'] > 0)
    df['has_any_introner'] = df['has_group1'] | df['has_group2']
    df['has_both_groups'] = df['has_group1'] & df['has_group2']
    
    # Left panel: Stacked bar chart of introner types
    introner_counts = {
        'Group1-Fixed': (df['group1_fixed'] > 0).sum(),
        'Group1-Polymorphic': (df['group1_polymorphic'] > 0).sum(),
        'Group2-Fixed': (df['group2_fixed'] > 0).sum(),
        'Group2-Polymorphic': (df['group2_polymorphic'] > 0).sum()
    }
    
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    bars = ax1.bar(range(len(introner_counts)), list(introner_counts.values()), 
                   color=colors, alpha=0.8, edgecolor='black', linewidth=0.8)
    
    ax1.set_xlabel('Introner Type', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax1.set_title('Introner Distribution by Type', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(introner_counts)))
    ax1.set_xticklabels(list(introner_counts.keys()), rotation=45, ha='right')
    
    # Add value labels on bars
    # for bar, count in zip(bars, introner_counts.values()):
    #     height = bar.get_height()
    #     ax1.text(bar.get_x() + bar.get_width()/2., height + 20,
    #             f'{count}', ha='center', va='bottom', fontweight='bold')
    
    ax1.grid(axis='y', alpha=0.3)
    
    # Right panel: Pie chart of overall introner presence
    # Total M. pusilla genes in genome (from annotation file)
    TOTAL_GENOME_GENES = 10660
    genes_with_introners = df['has_any_introner'].sum()
    genes_without_introners = TOTAL_GENOME_GENES - genes_with_introners
    
    introner_presence = {
        'Genes with Introners': genes_with_introners,
        'Genes without Introners': genes_without_introners
    }
    
    colors_pie = ['#2E86AB', '#E0E0E0']
    wedges, texts, autotexts = ax2.pie(list(introner_presence.values()), 
                                      labels=list(introner_presence.keys()),
                                      colors=colors_pie, autopct='%1.1f%%',
                                      startangle=90, textprops={'fontsize': 11})
    
    ax2.set_title('Overall Introner Prevalence', fontsize=14, fontweight='bold')
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase1_overview.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return introner_counts, introner_presence

def plot_group_comparison(df, output_dir):
    """
    Plot 2: Group Comparison
    Venn diagram + bar chart showing overlap between Group1 and Group2.
    """
    print("Creating group comparison plot...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: Simple overlap visualization using bar chart
    group1_genes = set(df[df['has_group1']]['gene_id'])
    group2_genes = set(df[df['has_group2']]['gene_id'])
    
    overlap = len(group1_genes & group2_genes)
    group1_only = len(group1_genes - group2_genes)
    group2_only = len(group2_genes - group1_genes)
    
    categories = ['Group1\nOnly', 'Overlap', 'Group2\nOnly']
    values = [group1_only, overlap, group2_only]
    colors = ['#2E86AB', '#A23B72', '#F18F01']
    
    bars = ax1.bar(categories, values, color=colors, alpha=0.8, edgecolor='black')
    ax1.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax1.set_title('Group1 vs Group2 Gene Overlap', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, value in zip(bars, values):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 10,
                f'{value}', ha='center', va='bottom', fontweight='bold')
    
    # Right panel: Fixed vs Polymorphic comparison
    group_data = {
        'Group1': [
            (df['group1_fixed'] > 0).sum(),
            (df['group1_polymorphic'] > 0).sum()
        ],
        'Group2': [
            (df['group2_fixed'] > 0).sum(),
            (df['group2_polymorphic'] > 0).sum()
        ]
    }
    
    x = np.arange(2)
    width = 0.35
    
    bars1 = ax2.bar(x - width/2, group_data['Group1'], width, 
                    label='Group1', color='#2E86AB', alpha=0.8, edgecolor='black')
    bars2 = ax2.bar(x + width/2, group_data['Group2'], width,
                    label='Group2', color='#F18F01', alpha=0.8, edgecolor='black')
    
    ax2.set_xlabel('Introner Pattern', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax2.set_title('Fixed vs Polymorphic Patterns', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Fixed', 'Polymorphic'])
    ax2.legend(fontsize=11)
    ax2.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 10,
                    f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase1_venn.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return len(group1_genes), len(group2_genes), len(group1_genes & group2_genes)

def plot_go_coverage(df, output_dir, go_json=None):
    """
    Plot 3: GO Annotation Coverage
    Bar chart + histogram showing GO annotation coverage and distribution.
    """
    print("Creating GO annotation coverage plot...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Load full genome GO data to calculate coverage for genes without introners
    try:
        import json
        if go_json is None:
            raise FileNotFoundError("No GO JSON file provided")
        with open(go_json, 'r') as f:
            full_go_data = json.load(f)

        # Create set of genes with introners from our dataset
        genes_with_introners = set(df['gene_id'])

        # Find genes without introners (genes in full GO data but not in our introner dataset)
        all_genes = set(full_go_data.keys())
        genes_without_introners = all_genes - genes_with_introners

        # Calculate GO coverage for genes without introners
        # A gene has GO terms if its array is non-empty
        genes_without_introners_with_go = sum(1 for gene in genes_without_introners
                                             if len(full_go_data[gene]) > 0)
        no_introners_coverage = (genes_without_introners_with_go / len(genes_without_introners) * 100) if genes_without_introners else 0

        print(f"  Found {len(genes_without_introners)} genes without introners")
        print(f"  {genes_without_introners_with_go} of them have GO annotations ({no_introners_coverage:.1f}%)")

    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"  Warning: Could not load full GO data ({e}). Excluding 'No Introners' category.")
        no_introners_coverage = None

    # Left panel: GO coverage by introner presence (simplified)
    go_coverage = {}

    # Calculate coverage for genes with introners (all genes in our dataset have introners)
    genes_with_introners_coverage = df['has_go_terms'].mean() * 100
    go_coverage['Genes with Introners'] = genes_with_introners_coverage

    # Add coverage for genes without introners if available
    if no_introners_coverage is not None:
        go_coverage['Genes without Introners'] = no_introners_coverage
    
    # Use appropriate colors for the simplified categories
    if no_introners_coverage is not None:
        colors = ['#2E86AB', '#E0E0E0']  # Blue for with introners, gray for without
    else:
        colors = ['#2E86AB']  # Just blue for with introners

    bars = ax1.bar(range(len(go_coverage)), list(go_coverage.values()),
                   color=colors, alpha=0.8, edgecolor='black', linewidth=0.8)
    
    ax1.set_xlabel('Gene Category', fontsize=12, fontweight='bold')
    ax1.set_ylabel('GO Annotation Coverage (%)', fontsize=12, fontweight='bold')
    ax1.set_title('GO Annotation Coverage by Category', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(go_coverage)))
    ax1.set_xticklabels(list(go_coverage.keys()), rotation=45, ha='right')
    ax1.set_ylim(0, 100)
    
    # Add percentage labels on bars
    # for bar, coverage in zip(bars, go_coverage.values()):
    #     height = bar.get_height()
    #     ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
    #             f'{coverage:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    ax1.grid(axis='y', alpha=0.3)
    
    # Right panel: Histogram of GO terms per gene
    go_annotated = df[df['has_go_terms']]
    if len(go_annotated) > 0:
        ax2.hist(go_annotated['n_go_terms'], bins=20, color='#2E86AB', 
                alpha=0.7, edgecolor='black', linewidth=0.8)
        
        ax2.set_xlabel('Number of GO Terms per Gene', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
        ax2.set_title('Distribution of GO Terms per Gene', fontsize=14, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add statistics as text
        mean_go = go_annotated['n_go_terms'].mean()
        median_go = go_annotated['n_go_terms'].median()
        ax2.axvline(mean_go, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_go:.1f}')
        ax2.axvline(median_go, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_go:.1f}')
        ax2.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'phase1_go_coverage.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return go_coverage

def generate_summary_stats(df, output_dir):
    """Generate summary statistics text file."""
    print("Generating summary statistics...")
    
    stats_file = os.path.join(output_dir, 'phase1_summary_stats.txt')
    
    # Total M. pusilla genes in genome
    TOTAL_GENOME_GENES = 10660
    genes_with_introners = df['has_any_introner'].sum()
    genes_without_introners = TOTAL_GENOME_GENES - genes_with_introners
    
    with open(stats_file, 'w') as f:
        f.write("=== PHASE 1 SUMMARY STATISTICS ===\n\n")
        
        f.write(f"Total M. pusilla genes in genome: {TOTAL_GENOME_GENES}\n")
        f.write(f"Genes analyzed (with introners): {len(df)}\n")
        f.write(f"Genes with any introners: {genes_with_introners} ({genes_with_introners/TOTAL_GENOME_GENES*100:.1f}% of genome)\n")
        f.write(f"Genes without introners: {genes_without_introners} ({genes_without_introners/TOTAL_GENOME_GENES*100:.1f}% of genome)\n\n")
        f.write(f"Genes with Group1 introners: {df['has_group1'].sum()}\n")
        f.write(f"Genes with Group2 introners: {df['has_group2'].sum()}\n")
        f.write(f"Genes with both groups: {df['has_both_groups'].sum()}\n\n")
        
        f.write("Introner type breakdown:\n")
        f.write(f"  Group1-Fixed: {(df['group1_fixed'] > 0).sum()} genes\n")
        f.write(f"  Group1-Polymorphic: {(df['group1_polymorphic'] > 0).sum()} genes\n")
        f.write(f"  Group2-Fixed: {(df['group2_fixed'] > 0).sum()} genes\n")
        f.write(f"  Group2-Polymorphic: {(df['group2_polymorphic'] > 0).sum()} genes\n\n")
        
        f.write("GO annotation coverage:\n")
        f.write(f"  Total genes with GO terms: {df['has_go_terms'].sum()} ({df['has_go_terms'].mean()*100:.1f}%)\n")
        
        if df['has_go_terms'].any():
            go_annotated = df[df['has_go_terms']]
            f.write(f"  Average GO terms per annotated gene: {go_annotated['n_go_terms'].mean():.1f}\n")
            f.write(f"  Median GO terms per annotated gene: {go_annotated['n_go_terms'].median():.1f}\n")
            f.write(f"  Max GO terms for single gene: {go_annotated['n_go_terms'].max()}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Phase 1 Results Visualization - Introner Insertion Pattern Analysis'
    )
    parser.add_argument('phase1_file', help='Phase 1 gene classification TSV')
    parser.add_argument('output_dir', help='Output directory for plots')
    parser.add_argument('--go_json', default=None,
                        help='Path to gene2go JSON file (optional, for coverage comparison)')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    print("=== PHASE 1 RESULTS VISUALIZATION ===")
    print(f"Input file: {args.phase1_file}")
    print(f"Output directory: {args.output_dir}")
    if args.go_json:
        print(f"GO JSON: {args.go_json}")
    print()

    try:
        # Load data
        df = load_phase1_data(args.phase1_file)

        # Generate plots
        introner_counts, introner_presence = plot_introner_distribution(df, args.output_dir)
        group1_size, group2_size, overlap_size = plot_group_comparison(df, args.output_dir)
        go_coverage = plot_go_coverage(df, args.output_dir, go_json=args.go_json)
        
        # Generate summary statistics
        generate_summary_stats(df, args.output_dir)

        print("\n=== VISUALIZATION SUMMARY ===")
        print(f"Generated 3 plots in {args.output_dir}")
        print(f"  • phase1_overview.pdf")
        print(f"  • phase1_venn.pdf")
        print(f"  • phase1_go_coverage.pdf")
        print(f"  • phase1_summary_stats.txt")
        
        # Calculate genome-wide statistics
        TOTAL_GENOME_GENES = 10660
        genes_with_introners = df['has_any_introner'].sum()
        
        print(f"\nKey findings:")
        print(f"  • {genes_with_introners}/{TOTAL_GENOME_GENES} genes contain introners ({genes_with_introners/TOTAL_GENOME_GENES*100:.1f}% of genome)")
        print(f"  • {overlap_size} genes have both Group1 and Group2 introners")
        print(f"  • {df['has_go_terms'].mean()*100:.1f}% GO annotation coverage among analyzed genes")
        
        print("\n✓ Phase 1 visualization completed successfully!")
        
    except Exception as e:
        print(f"✗ Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
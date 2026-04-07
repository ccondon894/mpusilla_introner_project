#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys

def calculate_dn_ds_ratios():
    """
    Calculate dN/dS ratios for each ortholog group comparing introner-present vs introner-absent samples
    """
    print("=== Phase 2: dN/dS Ratio Analysis ===")
    
    # Load the final dataset
    print("Loading selection dataset...")
    df = pd.read_csv('/scratch1/chris/introner_vis/final_selection_dataset.tsv', sep='\t')
    print(f"Loaded {len(df)} sample-variant records")
    
    # Filter for coding variants only (synonymous and missense/nonsense)
    coding_variants = df[df['mutation_type'].isin(['synonymous', 'missense', 'nonsense'])].copy()
    print(f"Coding variants: {len(coding_variants)}")
    
    # Calculate mutation counts per ortholog group and introner status
    print("Calculating mutation counts per ortholog group...")
    
    results = []
    
    for ortholog_id in coding_variants['ortholog_id'].unique():
        ortholog_data = coding_variants[coding_variants['ortholog_id'] == ortholog_id]
        
        # Split by introner status
        present_data = ortholog_data[ortholog_data['introner_status'] == 'present']
        absent_data = ortholog_data[ortholog_data['introner_status'] == 'absent']
        
        # Only analyze if we have data for both states
        if len(present_data) == 0 or len(absent_data) == 0:
            continue
            
        # Calculate counts for introner-present samples
        present_synonymous = len(present_data[present_data['mutation_type'] == 'synonymous'])
        present_nonsynonymous = len(present_data[present_data['mutation_type'].isin(['missense', 'nonsense'])])
        present_variants = len(present_data[present_data['has_variant'] == True])
        present_total = len(present_data)
        
        # Calculate counts for introner-absent samples  
        absent_synonymous = len(absent_data[absent_data['mutation_type'] == 'synonymous'])
        absent_nonsynonymous = len(absent_data[absent_data['mutation_type'].isin(['missense', 'nonsense'])])
        absent_variants = len(absent_data[absent_data['has_variant'] == True])
        absent_total = len(absent_data)
        
        # Calculate dN/dS ratios (add pseudocount to avoid division by zero)
        present_dn_ds = (present_nonsynonymous + 0.1) / (present_synonymous + 0.1)
        absent_dn_ds = (absent_nonsynonymous + 0.1) / (absent_synonymous + 0.1)
        
        # Calculate variant rates
        present_variant_rate = present_variants / present_total if present_total > 0 else 0
        absent_variant_rate = absent_variants / absent_total if absent_total > 0 else 0
        
        results.append({
            'ortholog_id': ortholog_id,
            'present_synonymous': present_synonymous,
            'present_nonsynonymous': present_nonsynonymous,
            'present_total': present_total,
            'present_variants': present_variants,
            'present_dn_ds': present_dn_ds,
            'present_variant_rate': present_variant_rate,
            'absent_synonymous': absent_synonymous,
            'absent_nonsynonymous': absent_nonsynonymous,
            'absent_total': absent_total,
            'absent_variants': absent_variants,
            'absent_dn_ds': absent_dn_ds,
            'absent_variant_rate': absent_variant_rate,
            'dn_ds_ratio': absent_dn_ds / present_dn_ds if present_dn_ds > 0 else np.nan
        })
    
    results_df = pd.DataFrame(results)
    
    print(f"\nAnalyzed {len(results_df)} ortholog groups with both introner states")
    
    # Summary statistics
    print("\n=== dN/dS Ratio Summary ===")
    print(f"Mean dN/dS (introner present): {results_df['present_dn_ds'].mean():.4f}")
    print(f"Mean dN/dS (introner absent): {results_df['absent_dn_ds'].mean():.4f}")
    print(f"Median dN/dS (introner present): {results_df['present_dn_ds'].median():.4f}")
    print(f"Median dN/dS (introner absent): {results_df['absent_dn_ds'].median():.4f}")
    
    # Statistical test
    present_dn_ds = results_df['present_dn_ds'].dropna()
    absent_dn_ds = results_df['absent_dn_ds'].dropna()
    
    # Wilcoxon signed-rank test (paired)
    valid_pairs = results_df.dropna(subset=['present_dn_ds', 'absent_dn_ds'])
    if len(valid_pairs) > 5:
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(valid_pairs['present_dn_ds'], valid_pairs['absent_dn_ds'])
        print(f"\nWilcoxon signed-rank test:")
        print(f"Statistic: {wilcoxon_stat}")
        print(f"P-value: {wilcoxon_p:.6f}")
        
        # Effect size (Cohen's d for paired data)
        differences = valid_pairs['absent_dn_ds'] - valid_pairs['present_dn_ds']
        effect_size = differences.mean() / differences.std()
        print(f"Effect size (Cohen's d): {effect_size:.4f}")
    
    # Mann-Whitney U test (unpaired)
    if len(present_dn_ds) > 0 and len(absent_dn_ds) > 0:
        mw_stat, mw_p = stats.mannwhitneyu(present_dn_ds, absent_dn_ds, alternative='two-sided')
        print(f"\nMann-Whitney U test:")
        print(f"Statistic: {mw_stat}")
        print(f"P-value: {mw_p:.6f}")
    
    # Save results
    output_file = '/scratch1/chris/introner_vis/dn_ds_analysis_results.tsv'
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Create visualizations
    create_dn_ds_plots(results_df)
    
    return results_df

def create_dn_ds_plots(results_df):
    """
    Create visualizations for dN/dS analysis
    """
    print("\nCreating dN/dS visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('dN/dS Ratio Analysis: Introner Present vs Absent', fontsize=16)
    
    # 1. Box plot comparison
    ax1 = axes[0, 0]
    data_for_box = []
    labels_for_box = []
    
    present_dn_ds = results_df['present_dn_ds'].dropna()
    absent_dn_ds = results_df['absent_dn_ds'].dropna()
    
    data_for_box.extend([present_dn_ds, absent_dn_ds])
    labels_for_box.extend(['Introner Present', 'Introner Absent'])
    
    bp = ax1.boxplot(data_for_box, labels=labels_for_box, patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    
    ax1.set_ylabel('dN/dS Ratio')
    ax1.set_title('dN/dS Distribution by Introner Status')
    ax1.grid(True, alpha=0.3)
    
    # 2. Violin plot
    ax2 = axes[0, 1]
    plot_data = pd.DataFrame({
        'dN/dS': list(present_dn_ds) + list(absent_dn_ds),
        'Introner Status': ['Present'] * len(present_dn_ds) + ['Absent'] * len(absent_dn_ds)
    })
    
    sns.violinplot(data=plot_data, x='Introner Status', y='dN/dS', ax=ax2)
    ax2.set_title('dN/dS Distribution (Violin Plot)')
    ax2.grid(True, alpha=0.3)
    
    # 3. Paired comparison plot
    ax3 = axes[1, 0]
    valid_pairs = results_df.dropna(subset=['present_dn_ds', 'absent_dn_ds'])
    
    if len(valid_pairs) > 0:
        ax3.scatter(valid_pairs['present_dn_ds'], valid_pairs['absent_dn_ds'], alpha=0.6)
        
        # Add diagonal line
        min_val = min(valid_pairs['present_dn_ds'].min(), valid_pairs['absent_dn_ds'].min())
        max_val = max(valid_pairs['present_dn_ds'].max(), valid_pairs['absent_dn_ds'].max())
        ax3.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7, label='No difference')
        
        ax3.set_xlabel('dN/dS (Introner Present)')
        ax3.set_ylabel('dN/dS (Introner Absent)')
        ax3.set_title('Paired dN/dS Comparison')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 4. Histogram of differences
    ax4 = axes[1, 1]
    if len(valid_pairs) > 0:
        differences = valid_pairs['absent_dn_ds'] - valid_pairs['present_dn_ds']
        ax4.hist(differences, bins=20, alpha=0.7, color='purple', edgecolor='black')
        ax4.axvline(x=0, color='red', linestyle='--', label='No difference')
        ax4.axvline(x=differences.mean(), color='orange', linestyle='-', label=f'Mean diff: {differences.mean():.4f}')
        ax4.set_xlabel('dN/dS Difference (Absent - Present)')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Distribution of dN/dS Differences')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = '/scratch1/chris/introner_vis/dn_ds_analysis_plots.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plots saved to: {plot_file}")
    
    plt.close()

if __name__ == "__main__":
    result = calculate_dn_ds_ratios()
    if result is not None:
        print(f"\n✓ dN/dS analysis completed successfully")
    else:
        print("\n✗ dN/dS analysis failed")
        sys.exit(1)
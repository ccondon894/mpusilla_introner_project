#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create box plots for all-samples diversity analysis')
    parser.add_argument('--input', help='Input file with all-samples diversity metrics')
    parser.add_argument('--output', help='Output file for box plot')
    parser.add_argument('--flank_length', help='Flanking sequence length for plot titles')
    return parser.parse_args()

def extract_diversity_categories(df):
    """Extract diversity values for each category based on fixation patterns"""
    
    categories_data = {
        'pi_group1_fixed': [],
        'pi_group1_absent': [],
        'pi_group2_fixed': [],
        'pi_group2_absent': [],
        'Dxy_group1_present_group2_absent': [],
        'Dxy_group1_absent_group2_present': [],
        'Dxy_group1_present_group2_present': []
    }
    
    for _, row in df.iterrows():
        category = row['category']
        
        # Extract pi_group1 values based on fixation state
        if category in ['group1_fixed_group2_absent', 'group1_fixed_group2_fixed']:
            # Group1 is fixed (present), so we have diversity within Group1
            if pd.notna(row['pi_group1']):
                categories_data['pi_group1_fixed'].append(row['pi_group1'])
        elif category == 'group1_absent_group2_fixed':
            # Group1 is absent, so we have diversity within Group1 absent lineages
            if pd.notna(row['pi_group1']):
                categories_data['pi_group1_absent'].append(row['pi_group1'])
        
        # Extract pi_group2 values based on fixation state
        if category in ['group1_absent_group2_fixed', 'group1_fixed_group2_fixed']:
            # Group2 is fixed (present), so we have diversity within Group2
            if pd.notna(row['pi_group2']):
                categories_data['pi_group2_fixed'].append(row['pi_group2'])
        elif category == 'group1_fixed_group2_absent':
            # Group2 is absent, so we have diversity within Group2 absent lineages
            if pd.notna(row['pi_group2']):
                categories_data['pi_group2_absent'].append(row['pi_group2'])
        
        # Extract Dxy values based on specific category combinations
        if pd.notna(row['dxy_group1_group2']):
            if category == 'group1_fixed_group2_absent':
                categories_data['Dxy_group1_present_group2_absent'].append(row['dxy_group1_group2'])
            elif category == 'group1_absent_group2_fixed':
                categories_data['Dxy_group1_absent_group2_present'].append(row['dxy_group1_group2'])
            elif category == 'group1_fixed_group2_fixed':
                categories_data['Dxy_group1_present_group2_present'].append(row['dxy_group1_group2'])
    
    return categories_data

def write_summary_statistics(categories_data, flank_length):
    """Write detailed summary statistics to a file"""
    
    # Create the summary file path
    summary_dir = "all_samples_evolution_analysis/diversity_metrics"
    os.makedirs(summary_dir, exist_ok=True)
    summary_file = os.path.join(summary_dir, f"box_plot_summary_stats_{flank_length}bp.txt")
    
    with open(summary_file, 'w') as f:
        f.write(f"All-Samples Diversity Analysis - Box Plot Summary Statistics\n")
        f.write(f"Flanking sequence length: {flank_length}bp\n")
        f.write("="*80 + "\n\n")
        
        # Group categories by type
        pi_categories = ['pi_group1_fixed', 'pi_group1_absent', 'pi_group2_fixed', 'pi_group2_absent']
        dxy_categories = ['Dxy_group1_present_group2_absent', 'Dxy_group1_absent_group2_present', 'Dxy_group1_present_group2_present']
        
        # Write pi statistics
        f.write("NUCLEOTIDE DIVERSITY (π) STATISTICS\n")
        f.write("-" * 40 + "\n")
        for category in pi_categories:
            values = categories_data[category]
            display_name = category.replace('_', ' ').replace('pi ', 'π ')
            
            if len(values) > 0:
                mean_val = np.mean(values)
                std_val = np.std(values)
                median_val = np.median(values)
                q25 = np.percentile(values, 25)
                q75 = np.percentile(values, 75)
                min_val = np.min(values)
                max_val = np.max(values)
                
                f.write(f"\n{display_name}:\n")
                f.write(f"  Sample size (n): {len(values)}\n")
                f.write(f"  Mean: {mean_val:.6f}\n")
                f.write(f"  Standard deviation: {std_val:.6f}\n")
                f.write(f"  Median: {median_val:.6f}\n")
                f.write(f"  25th percentile: {q25:.6f}\n")
                f.write(f"  75th percentile: {q75:.6f}\n")
                f.write(f"  Minimum: {min_val:.6f}\n")
                f.write(f"  Maximum: {max_val:.6f}\n")
            else:
                f.write(f"\n{display_name}:\n")
                f.write(f"  Sample size (n): 0 (no data)\n")
        
        # Write dxy statistics
        f.write(f"\n\n" + "="*80 + "\n")
        f.write("DIVERGENCE (dXY) STATISTICS\n")
        f.write("-" * 40 + "\n")
        for category in dxy_categories:
            values = categories_data[category]
            display_name = category.replace('_', ' ').replace('Dxy', 'dXY')
            
            if len(values) > 0:
                mean_val = np.mean(values)
                std_val = np.std(values)
                median_val = np.median(values)
                q25 = np.percentile(values, 25)
                q75 = np.percentile(values, 75)
                min_val = np.min(values)
                max_val = np.max(values)
                
                f.write(f"\n{display_name}:\n")
                f.write(f"  Sample size (n): {len(values)}\n")
                f.write(f"  Mean: {mean_val:.6f}\n")
                f.write(f"  Standard deviation: {std_val:.6f}\n")
                f.write(f"  Median: {median_val:.6f}\n")
                f.write(f"  25th percentile: {q25:.6f}\n")
                f.write(f"  75th percentile: {q75:.6f}\n")
                f.write(f"  Minimum: {min_val:.6f}\n")
                f.write(f"  Maximum: {max_val:.6f}\n")
            else:
                f.write(f"\n{display_name}:\n")
                f.write(f"  Sample size (n): 0 (no data)\n")
    
    print(f"Created summary statistics file: {summary_file}")
    return summary_file

def plot_all_samples_box_plots(df, output_file, flank_length):
    """Create comprehensive box plots for all diversity categories with separate subplots for pi and dxy"""
    
    # Extract diversity values for each category
    categories_data = extract_diversity_categories(df)
    
    # Define colors for each category type
    color_map = {
        'pi_group1_fixed': '#3b528b',
        'pi_group1_absent': '#21918c',
        'pi_group2_fixed': '#B91C1C',
        'pi_group2_absent': '#DC2626',
        'Dxy_group1_present_group2_absent': '#3b528b',
        'Dxy_group1_absent_group2_present': '#B91C1C',
        'Dxy_group1_present_group2_present': '#5ec962'
    }
    
    # Separate pi and dxy categories
    pi_categories = ['pi_group1_fixed', 'pi_group1_absent', 'pi_group2_fixed', 'pi_group2_absent']
    dxy_categories = ['Dxy_group1_present_group2_absent', 'Dxy_group1_absent_group2_present', 'Dxy_group1_present_group2_present']
    
    # Prepare data for pi subplot
    pi_data = []
    pi_labels = []
    pi_colors = []
    for category in pi_categories:
        values = categories_data[category]
        if len(values) > 0:
            pi_data.append(values)
            # Clean label: remove 'pi_' prefix and replace underscores with spaces
            clean_label = category.replace('pi_', '').replace('_', ' ')
            pi_labels.append(clean_label)
            pi_colors.append(color_map[category])
    
    pi_labels = ['Group1 fixed', 'Group1 absent', 'Group2 fixed', 'Group2 absent']
    # Prepare data for dxy subplot
    dxy_data = []
    dxy_labels = []
    dxy_colors = []
    for category in dxy_categories:
        values = categories_data[category]
        if len(values) > 0:
            dxy_data.append(values)
            # Clean label: remove 'Dxy_' prefix and replace underscores with spaces
            clean_label = category.replace('Dxy_', '').replace('_', ' ')
            dxy_labels.append(clean_label)
            dxy_colors.append(color_map[category])
    
    dxy_labels = ["Group1 present &\nGroup2 absent", "Group1 absent &\nGroup2 present", "Group1 present &\nGroup2 present"]
    if not pi_data and not dxy_data:
        print("Warning: No data available for plotting")
        return
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Subplot 1: Pi swarm plots
    if pi_data:
        # Create swarm plot data in a format suitable for seaborn
        pi_plot_data = []
        for i, (data, label) in enumerate(zip(pi_data, pi_labels)):
            for value in data:
                pi_plot_data.append({'category': label, 'value': value})

        pi_df = pd.DataFrame(pi_plot_data)

        # Create jitter plots with custom colors
        for i, (label, color) in enumerate(zip(pi_labels, pi_colors)):
            data_subset = pi_df[pi_df['category'] == label]
            sns.stripplot(x='category', y='value', data=data_subset,
                         color=color, edgecolor='black', linewidth=0.8,
                         size=5, ax=ax1, zorder=2, jitter=True, alpha=0.6)

        ax1.set_xticks(range(len(pi_labels)))
        ax1.set_xticklabels(pi_labels, rotation=45, ha='right', fontsize=15)

        # Plot medians as black squares
        pi_medians = [np.median(data) for data in pi_data]
        pi_positions = range(len(pi_data))
        ax1.scatter(pi_positions, pi_medians, color='black', marker='s', s=120,
                   label='Median', zorder=4, edgecolor='white', linewidth=1)

        # Calculate and plot means as diamonds
        pi_means = [np.mean(data) for data in pi_data]
        ax1.scatter(pi_positions, pi_means, color='black', marker='D', s=100,
                   label='Mean', zorder=4, edgecolor='white', linewidth=1)

        # Formatting for pi subplot
        ax1.set_xlabel('')
        ax1.set_ylabel('Nucleotide Diversity (π)', fontsize=16)
        ax1.tick_params(labelsize=15)
        ax1.set_ylim(bottom=0, top=0.1)
        ax1.grid(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        
    else:
        ax1.text(0.5, 0.5, 'No π data available', ha='center', va='center', transform=ax1.transAxes)
    
    # Subplot 2: Dxy swarm plots
    if dxy_data:
        # Create swarm plot data in a format suitable for seaborn
        dxy_plot_data = []
        for i, (data, label) in enumerate(zip(dxy_data, dxy_labels)):
            for value in data:
                dxy_plot_data.append({'category': label, 'value': value})

        dxy_df = pd.DataFrame(dxy_plot_data)

        # Create jitter plots with custom colors
        for i, (label, color) in enumerate(zip(dxy_labels, dxy_colors)):
            data_subset = dxy_df[dxy_df['category'] == label]
            sns.stripplot(x='category', y='value', data=data_subset,
                         color=color, edgecolor='black', linewidth=0.8,
                         size=5, ax=ax2, zorder=2, jitter=True, alpha=0.6)

        ax2.set_xticks(range(len(dxy_labels)))
        ax2.set_xticklabels(dxy_labels, rotation=45, ha='right', fontsize=15)

        # Plot medians as black squares
        dxy_medians = [np.median(data) for data in dxy_data]
        dxy_positions = range(len(dxy_data))
        ax2.scatter(dxy_positions, dxy_medians, color='black', marker='s', s=120,
                   label='Median', zorder=4, edgecolor='white', linewidth=1)

        # Calculate and plot means as diamonds
        dxy_means = [np.mean(data) for data in dxy_data]
        ax2.scatter(dxy_positions, dxy_means, color='black', marker='D', s=100,
                   label='Mean', zorder=4, edgecolor='white', linewidth=1)

        # Formatting for dxy subplot
        ax2.set_xlabel('')
        ax2.set_ylabel('Divergence (dXY)', fontsize=16)
        ax2.tick_params(labelsize=15)
        ax2.legend(loc='upper right', fontsize=10, title_fontsize=11,
                  frameon=True, edgecolor='black')
        ax2.set_ylim(bottom=0, top=0.4)
        ax2.grid(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        
    else:
        ax2.text(0.5, 0.5, 'No dXY data available', ha='center', va='center', transform=ax2.transAxes)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created all-samples box plot: {output_file}")
    
    # Print summary statistics
    print("\nSummary statistics:")
    for category, values in categories_data.items():
        if len(values) > 0:
            mean_val = np.mean(values)
            std_val = np.std(values)
            median_val = np.median(values)
            print(f"{category}: n={len(values)}, mean={mean_val:.6f}, std={std_val:.6f}, median={median_val:.6f}")
        else:
            print(f"{category}: n=0 (no data)")

def main():
    args = parse_arguments()
    
    # Load diversity metrics
    print(f"Loading all-samples diversity metrics from {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    print(f"Loaded {len(df)} ortholog groups")
    
    # Print category distribution
    if 'category' in df.columns:
        print("\nCategory distribution:")
        category_counts = df['category'].value_counts()
        for category, count in category_counts.items():
            print(f"  {category}: {count} ortholog groups")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Create the box plots and get categories data for summary
    categories_data = extract_diversity_categories(df)
    plot_all_samples_box_plots(df, args.output, args.flank_length)
    
    # Write summary statistics file
    write_summary_statistics(categories_data, args.flank_length)
    
    print("Done!")

if __name__ == "__main__":
    main()
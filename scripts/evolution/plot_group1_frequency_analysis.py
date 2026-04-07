#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import stats

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create plots for Group1 frequency-based diversity analysis')
    parser.add_argument('--input', help='Input file with Group1 diversity metrics')
    parser.add_argument('--output', help='Output file prefix for plots')
    parser.add_argument('--plot_type', choices=['frequency_lines', 'box_plot', 'combined_box_plot'],
                       help='Type of plot to create')
    parser.add_argument('--flank_length', help='Flanking sequence length for plot titles')
    return parser.parse_args()

def plot_frequency_line_plots(df, output_prefix, flank_length):
    """Create line plots showing both π(present) and π(absent) polarized by present and absent frequencies"""
    
    # Calculate summary statistics by frequency (all metrics together)
    freq_summary = df.groupby('frequency').agg({
        'pi_present': ['mean', 'std', 'count'],
        'pi_absent': ['mean', 'std', 'count'],
        'pi_between': ['mean', 'std', 'count']
    }).reset_index()
    
    # Flatten column names
    freq_summary.columns = ['frequency'] + ['_'.join(col).strip('_') for col in freq_summary.columns.values[1:]]
    
    # Calculate Standard Error of the Mean (SEM) for error bars
    freq_summary['pi_present_sem'] = freq_summary['pi_present_std'] / np.sqrt(freq_summary['pi_present_count'])
    freq_summary['pi_absent_sem'] = freq_summary['pi_absent_std'] / np.sqrt(freq_summary['pi_absent_count'])
    freq_summary['pi_between_sem'] = freq_summary['pi_between_std'] / np.sqrt(freq_summary['pi_between_count'])
    
    # Plot 1: Present frequency polarization (π(present), π(absent), π(between) vs present frequency)
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6))
    
    frequencies = freq_summary['frequency']
    
    # Pi present
    present_mask = freq_summary['pi_present_count'] > 0
    if present_mask.any():
        ax1.errorbar(frequencies[present_mask], freq_summary['pi_present_mean'][present_mask], 
                    yerr=freq_summary['pi_present_sem'][present_mask],
                    marker='o', label='π(present)', linewidth=2, markersize=6, color='blue')
    
    # Pi absent 
    absent_mask = freq_summary['pi_absent_count'] > 0
    if absent_mask.any():
        ax1.errorbar(frequencies[absent_mask], freq_summary['pi_absent_mean'][absent_mask],
                    yerr=freq_summary['pi_absent_sem'][absent_mask], 
                    marker='s', label='π(absent)', linewidth=2, markersize=6, color='red')
    
    # Pi between
    between_mask = freq_summary['pi_between_count'] > 0
    if between_mask.any():
        ax1.errorbar(frequencies[between_mask], freq_summary['pi_between_mean'][between_mask],
                    yerr=freq_summary['pi_between_sem'][between_mask],
                    marker='^', label='π(between) = dXY', linewidth=2, markersize=6, color='green')
    
    ax1.set_xlabel('Allele Frequency (# present out of 11)')
    ax1.set_ylabel('Nucleotide Diversity (π)')
    ax1.set_title(f'Frequency Spectrum - Present Polarization ({flank_length}bp flanks)\nError bars show ± SEM')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(range(1, 11))
    ax1.set_ylim(bottom=0)
    
    plt.tight_layout()
    present_output = os.path.join(output_prefix, f"frequency_line_plots_{flank_length}bp_present_spectrum.png")
    plt.savefig(present_output, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Absent frequency polarization (same data, x-axis shows absent frequency)
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
    
    # Convert present frequencies to absent frequencies for x-axis
    absent_frequencies = 11 - frequencies
    
    # Pi present (same data, different x-axis)
    if present_mask.any():
        ax2.errorbar(absent_frequencies[present_mask], freq_summary['pi_present_mean'][present_mask], 
                    yerr=freq_summary['pi_present_sem'][present_mask],
                    marker='o', label='π(present)', linewidth=2, markersize=6, color='blue')
    
    # Pi absent (same data, different x-axis)
    if absent_mask.any():
        ax2.errorbar(absent_frequencies[absent_mask], freq_summary['pi_absent_mean'][absent_mask],
                    yerr=freq_summary['pi_absent_sem'][absent_mask], 
                    marker='s', label='π(absent)', linewidth=2, markersize=6, color='red')
    
    # Pi between (same data, different x-axis)
    if between_mask.any():
        ax2.errorbar(absent_frequencies[between_mask], freq_summary['pi_between_mean'][between_mask],
                    yerr=freq_summary['pi_between_sem'][between_mask],
                    marker='^', label='π(between) = dXY', linewidth=2, markersize=6, color='green')
    
    ax2.set_xlabel('Allele Frequency (# absent out of 11)')
    ax2.set_ylabel('Nucleotide Diversity (π)')
    ax2.set_title(f'Frequency Spectrum - Absent Polarization ({flank_length}bp flanks)\nError bars show ± SEM')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(range(1, 11))
    ax2.set_ylim(bottom=0)
    
    plt.tight_layout()
    absent_output = os.path.join(output_prefix, f"frequency_line_plots_{flank_length}bp_absent_spectrum.png")
    plt.savefig(absent_output, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created present frequency spectrum plot: {present_output}")
    print(f"Created absent frequency spectrum plot: {absent_output}")

def plot_box_plots(df, output_prefix, flank_length):
    """Create violin plots for present and absent allele frequency spectrums"""

    # Generate viridis colors for frequency bins
    viridis_colors = cm.viridis(np.linspace(0.2, 0.8, 10))

    # Filter data for valid pi values
    present_data = df[df['pi_present'].notna()].copy()
    absent_data = df[df['pi_absent'].notna()].copy()

    # Add absent frequency column for absent spectrum
    absent_data['absent_frequency'] = 11 - absent_data['frequency']

    # Plot 1: Present AFS violin plot
    fig1, ax1 = plt.subplots(1, 1, figsize=(14, 8))

    # Create violin plot for present frequencies
    if not present_data.empty:
        frequencies = sorted(present_data['frequency'].unique())
        box_data = [present_data[present_data['frequency'] == freq]['pi_present'].values
                   for freq in frequencies]

        # Create violin plots with custom colors
        box_positions = range(len(frequencies))
        for i, (pos, freq, data) in enumerate(zip(box_positions, frequencies, box_data)):
            parts = ax1.violinplot([data], positions=[pos], widths=0.7,
                                   showmeans=False, showmedians=False, showextrema=False)

            # Color the violin
            for pc in parts['bodies']:
                pc.set_facecolor(viridis_colors[i % len(viridis_colors)])
                pc.set_edgecolor('black')
                pc.set_alpha(0.7)
                pc.set_linewidth(1)

        # Calculate and plot medians as black squares
        medians = [present_data[present_data['frequency'] == freq]['pi_present'].median()
                  for freq in frequencies]
        ax1.scatter(box_positions, medians, color='black', marker='s', s=100,
                   label='Median', zorder=4, edgecolor='white', linewidth=0.5)

        # Calculate and plot means as black diamonds
        means = [present_data[present_data['frequency'] == freq]['pi_present'].mean()
                for freq in frequencies]
        ax1.scatter(box_positions, means, color='black', marker='D', s=80,
                   label='Mean', zorder=4, edgecolor='white', linewidth=0.5)

    ax1.set_xlabel('Allele Frequency (# present out of 11)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('π(present)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Present Allele Frequency Spectrum ({flank_length}bp flanks)',
                 fontsize=13, fontweight='bold')
    ax1.legend(loc='upper right', frameon=True, fancybox=True)
    ax1.set_xticks(box_positions)
    ax1.set_xticklabels(frequencies)
    ax1.set_ylim(bottom=0, top=0.04)

    plt.tight_layout()
    present_box_output = os.path.join(output_prefix, f"frequency_box_plots_{flank_length}bp_present_spectrum.png")
    plt.savefig(present_box_output, dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Absent AFS violin plot
    fig2, ax2 = plt.subplots(1, 1, figsize=(14, 8))

    # Create violin plot for absent frequencies
    if not absent_data.empty:
        absent_frequencies = sorted(absent_data['absent_frequency'].unique())
        box_data_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].values
                          for freq in absent_frequencies]

        # Create violin plots with custom colors
        box_positions_absent = range(len(absent_frequencies))
        for i, (pos, freq, data) in enumerate(zip(box_positions_absent, absent_frequencies, box_data_absent)):
            parts = ax2.violinplot([data], positions=[pos], widths=0.7,
                                   showmeans=False, showmedians=False, showextrema=False)

            # Color the violin
            for pc in parts['bodies']:
                pc.set_facecolor(viridis_colors[i % len(viridis_colors)])
                pc.set_edgecolor('black')
                pc.set_alpha(0.7)
                pc.set_linewidth(1)

        # Calculate and plot medians as black squares
        medians_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].median()
                         for freq in absent_frequencies]
        ax2.scatter(box_positions_absent, medians_absent, color='black', marker='s', s=100,
                   label='Median', zorder=4, edgecolor='white', linewidth=0.5)

        # Calculate and plot means as black diamonds
        means_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].mean()
                       for freq in absent_frequencies]
        ax2.scatter(box_positions_absent, means_absent, color='black', marker='D', s=80,
                   label='Mean', zorder=4, edgecolor='white', linewidth=0.5)

    ax2.set_xlabel('Allele Frequency (# absent out of 11)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('π(absent)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Absent Allele Frequency Spectrum ({flank_length}bp flanks)',
                 fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right', frameon=True, fancybox=True)
    ax2.set_xticks(box_positions_absent)
    ax2.set_xticklabels(absent_frequencies)
    ax2.set_ylim(bottom=0, top=0.04)

    plt.tight_layout()
    absent_box_output = os.path.join(output_prefix, f"frequency_box_plots_{flank_length}bp_absent_spectrum.png")
    plt.savefig(absent_box_output, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Created present AFS violin plot: {present_box_output}")
    print(f"Created absent AFS violin plot: {absent_box_output}")

def plot_combined_box_plots(df, output_prefix, flank_length):
    """Create combined violin plots with present and absent subplots side by side"""

    # Generate viridis colors for frequency bins
    viridis_colors = cm.viridis(np.linspace(0.2, 0.8, 10))

    # Filter data for valid pi values
    present_data = df[df['pi_present'].notna()].copy()
    absent_data = df[df['pi_absent'].notna()].copy()

    # Add absent frequency column for absent spectrum
    absent_data['absent_frequency'] = 11 - absent_data['frequency']

    # Create combined figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Left subplot: Present AFS violin plot
    if not present_data.empty:
        frequencies = sorted(present_data['frequency'].unique())
        box_data = [present_data[present_data['frequency'] == freq]['pi_present'].values
                   for freq in frequencies]

        # Create violin plots with custom colors
        box_positions = range(len(frequencies))
        for i, (pos, freq, data) in enumerate(zip(box_positions, frequencies, box_data)):
            parts = ax1.violinplot([data], positions=[pos], widths=0.7,
                                   showmeans=False, showmedians=False, showextrema=False)

            # Color the violin
            for pc in parts['bodies']:
                pc.set_facecolor(viridis_colors[i % len(viridis_colors)])
                pc.set_edgecolor('black')
                pc.set_alpha(0.7)
                pc.set_linewidth(1)

        # Calculate and plot medians as black squares
        medians = [present_data[present_data['frequency'] == freq]['pi_present'].median()
                  for freq in frequencies]
        ax1.scatter(box_positions, medians, color='black', marker='s', s=100,
                   label='Median', zorder=4, edgecolor='white', linewidth=0.5)

        # Calculate and plot means as black diamonds
        means = [present_data[present_data['frequency'] == freq]['pi_present'].mean()
                for freq in frequencies]
        ax1.scatter(box_positions, means, color='black', marker='D', s=80,
                   label='Mean', zorder=4, edgecolor='white', linewidth=0.5)

    ax1.set_xlabel('Allele Frequency (# present out of 11)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('π(present)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Present Allele Frequency Spectrum\n({flank_length}bp flanks)', fontsize=13, fontweight='bold')
    ax1.legend(loc='upper right', frameon=True, fancybox=True)
    ax1.set_xticks(box_positions)
    ax1.set_xticklabels(frequencies)
    ax1.set_ylim(bottom=0, top=0.04)

    # Right subplot: Absent AFS violin plot
    if not absent_data.empty:
        absent_frequencies = sorted(absent_data['absent_frequency'].unique())
        box_data_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].values
                          for freq in absent_frequencies]

        # Create violin plots with custom colors
        box_positions_absent = range(len(absent_frequencies))
        for i, (pos, freq, data) in enumerate(zip(box_positions_absent, absent_frequencies, box_data_absent)):
            parts = ax2.violinplot([data], positions=[pos], widths=0.7,
                                   showmeans=False, showmedians=False, showextrema=False)

            # Color the violin
            for pc in parts['bodies']:
                pc.set_facecolor(viridis_colors[i % len(viridis_colors)])
                pc.set_edgecolor('black')
                pc.set_alpha(0.7)
                pc.set_linewidth(1)

        # Calculate and plot medians as black squares
        medians_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].median()
                         for freq in absent_frequencies]
        ax2.scatter(box_positions_absent, medians_absent, color='black', marker='s', s=100,
                   label='Median', zorder=4, edgecolor='white', linewidth=0.5)

        # Calculate and plot means as black diamonds
        means_absent = [absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].mean()
                       for freq in absent_frequencies]
        ax2.scatter(box_positions_absent, means_absent, color='black', marker='D', s=80,
                   label='Mean', zorder=4, edgecolor='white', linewidth=0.5)

    ax2.set_xlabel('Allele Frequency (# absent out of 11)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('π(absent)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Absent Allele Frequency Spectrum\n({flank_length}bp flanks)', fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right', frameon=True, fancybox=True)
    ax2.set_xticks(box_positions_absent)
    ax2.set_xticklabels(absent_frequencies)
    ax2.set_ylim(bottom=0, top=0.04)

    plt.tight_layout()
    combined_output = os.path.join(output_prefix, f"frequency_box_plots_{flank_length}bp_combined_spectrum.png")
    plt.savefig(combined_output, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Created combined AFS violin plot: {combined_output}")
    return combined_output

def write_summary_statistics(df, flank_length):
    """Write detailed summary statistics for Group1 frequency analysis"""
    
    # Create the summary file path
    summary_dir = "group1_evolution_analysis/diversity_metrics"
    os.makedirs(summary_dir, exist_ok=True)
    summary_file = os.path.join(summary_dir, f"box_plot_summary_stats_{flank_length}bp.txt")
    
    with open(summary_file, 'w') as f:
        f.write(f"Group1 Diversity Analysis - Box Plot Summary Statistics\n")
        f.write(f"Flanking sequence length: {flank_length}bp\n")
        f.write("="*80 + "\n\n")
        
        # Process Present AFS data
        present_data = df[df['pi_present'].notna()].copy()
        
        f.write("PRESENT ALLELE FREQUENCY SPECTRUM - π(present) STATISTICS\n")
        f.write("-" * 60 + "\n")
        
        if not present_data.empty:
            frequencies = sorted(present_data['frequency'].unique())
            for freq in frequencies:
                freq_data = present_data[present_data['frequency'] == freq]['pi_present'].values
                
                if len(freq_data) > 0:
                    mean_val = np.mean(freq_data)
                    std_val = np.std(freq_data)
                    median_val = np.median(freq_data)
                    q25 = np.percentile(freq_data, 25)
                    q75 = np.percentile(freq_data, 75)
                    min_val = np.min(freq_data)
                    max_val = np.max(freq_data)
                    
                    f.write(f"\nFrequency {freq} (# present out of 11):\n")
                    f.write(f"  Sample size (n): {len(freq_data)}\n")
                    f.write(f"  Mean: {mean_val:.6f}\n")
                    f.write(f"  Standard deviation: {std_val:.6f}\n")
                    f.write(f"  Median: {median_val:.6f}\n")
                    f.write(f"  25th percentile: {q25:.6f}\n")
                    f.write(f"  75th percentile: {q75:.6f}\n")
                    f.write(f"  Minimum: {min_val:.6f}\n")
                    f.write(f"  Maximum: {max_val:.6f}\n")
        else:
            f.write("\nNo present AFS data available\n")
        
        # Process Absent AFS data  
        absent_data = df[df['pi_absent'].notna()].copy()
        absent_data['absent_frequency'] = 11 - absent_data['frequency']
        
        f.write(f"\n\n" + "="*80 + "\n")
        f.write("ABSENT ALLELE FREQUENCY SPECTRUM - π(absent) STATISTICS\n")
        f.write("-" * 59 + "\n")
        
        if not absent_data.empty:
            absent_frequencies = sorted(absent_data['absent_frequency'].unique())
            for freq in absent_frequencies:
                freq_data = absent_data[absent_data['absent_frequency'] == freq]['pi_absent'].values
                
                if len(freq_data) > 0:
                    mean_val = np.mean(freq_data)
                    std_val = np.std(freq_data)
                    median_val = np.median(freq_data)
                    q25 = np.percentile(freq_data, 25)
                    q75 = np.percentile(freq_data, 75)
                    min_val = np.min(freq_data)
                    max_val = np.max(freq_data)
                    
                    f.write(f"\nFrequency {freq} (# absent out of 11):\n")
                    f.write(f"  Sample size (n): {len(freq_data)}\n")
                    f.write(f"  Mean: {mean_val:.6f}\n")
                    f.write(f"  Standard deviation: {std_val:.6f}\n")
                    f.write(f"  Median: {median_val:.6f}\n")
                    f.write(f"  25th percentile: {q25:.6f}\n")
                    f.write(f"  75th percentile: {q75:.6f}\n")
                    f.write(f"  Minimum: {min_val:.6f}\n")
                    f.write(f"  Maximum: {max_val:.6f}\n")
        else:
            f.write("\nNo absent AFS data available\n")
        
        # Process π(between) data by frequency
        between_data = df[df['pi_between'].notna()].copy()
        
        f.write(f"\n\n" + "="*80 + "\n")
        f.write("π(between) STATISTICS BY FREQUENCY\n")
        f.write("-" * 34 + "\n")
        
        if not between_data.empty:
            frequencies = sorted(between_data['frequency'].unique())
            for freq in frequencies:
                freq_data = between_data[between_data['frequency'] == freq]['pi_between'].values
                
                if len(freq_data) > 0:
                    mean_val = np.mean(freq_data)
                    std_val = np.std(freq_data)
                    median_val = np.median(freq_data)
                    q25 = np.percentile(freq_data, 25)
                    q75 = np.percentile(freq_data, 75)
                    min_val = np.min(freq_data)
                    max_val = np.max(freq_data)
                    
                    f.write(f"\nFrequency {freq}:\n")
                    f.write(f"  Sample size (n): {len(freq_data)}\n")
                    f.write(f"  Mean: {mean_val:.6f}\n")
                    f.write(f"  Standard deviation: {std_val:.6f}\n")
                    f.write(f"  Median: {median_val:.6f}\n")
                    f.write(f"  25th percentile: {q25:.6f}\n")
                    f.write(f"  75th percentile: {q75:.6f}\n")
                    f.write(f"  Minimum: {min_val:.6f}\n")
                    f.write(f"  Maximum: {max_val:.6f}\n")
        else:
            f.write("\nNo π(between) data available\n")
    
    print(f"Created summary statistics file: {summary_file}")
    return summary_file

def main():
    args = parse_arguments()
    
    # Load diversity metrics
    print(f"Loading diversity metrics from {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    print(f"Loaded {len(df)} ortholog groups")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Create the appropriate plots
    if args.plot_type == 'frequency_lines':
        plot_frequency_line_plots(df, args.output, args.flank_length)
    elif args.plot_type == 'box_plot':
        plot_box_plots(df, args.output, args.flank_length)
        # Generate summary statistics for box plots
        write_summary_statistics(df, args.flank_length)
    elif args.plot_type == 'combined_box_plot':
        plot_combined_box_plots(df, args.output, args.flank_length)
        # Generate summary statistics for box plots
        write_summary_statistics(df, args.flank_length)
    else:
        print(f"Unknown plot type: {args.plot_type}")
    
    print("Done!")

if __name__ == "__main__":
    main()
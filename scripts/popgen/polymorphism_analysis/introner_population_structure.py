#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform, jaccard
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import math
import os
import sys

def create_output_directories():
    """Create directory structure for analysis outputs"""
    directories = [
        'introner_analysis',
        'introner_analysis/data_exploration',
        'introner_analysis/clustering_by_frequency',
        'introner_analysis/comparative_analysis',
        'introner_analysis/demographic_integration',
        'introner_analysis/summary'
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        
    # Create frequency-specific directories
    for freq in range(2, 11):
        freq_dir = f'introner_analysis/clustering_by_frequency/freq_{freq}'
        os.makedirs(freq_dir, exist_ok=True)
    
    print("Created analysis directory structure")
    
def load_and_prepare_data():
    """
    Load polymorphic groups data and create binary presence/absence matrix
    """
    print("=== Phase 1: Data Preparation and Exploration ===")
    
    # Load the polymorphic groups data
    print("Loading polymorphic groups data...")
    polymorphic_df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t')
    print(f"Loaded {len(polymorphic_df)} polymorphic introner groups")
    
    # Define Group1 samples
    group1_samples = ['CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
                      'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    print(f"Analyzing {len(group1_samples)} Group1 samples:")
    print(f"  {group1_samples}")
    
    # Create binary presence/absence matrix
    print("Creating binary presence/absence matrix...")
    binary_matrix = []
    introner_ids = []
    
    for idx, row in polymorphic_df.iterrows():
        ortholog_id = row['ortholog_id']
        present_samples = eval(row['present_samples']) if isinstance(row['present_samples'], str) else row['present_samples']
        
        # Create binary vector for this introner
        binary_row = []
        for sample in group1_samples:
            if sample in present_samples:
                binary_row.append(1)
            else:
                binary_row.append(0)
        
        binary_matrix.append(binary_row)
        introner_ids.append(ortholog_id)
    
    # Convert to DataFrame
    binary_df = pd.DataFrame(binary_matrix, 
                           columns=group1_samples, 
                           index=introner_ids)
    
    print(f"Created binary matrix: {binary_df.shape[0]} introners × {binary_df.shape[1]} samples")
    
    # Calculate frequencies
    frequencies = binary_df.sum(axis=1)
    polymorphic_df['frequency'] = frequencies.values
    
    # Summary statistics
    print(f"\nIntroner frequency distribution:")
    freq_counts = frequencies.value_counts().sort_index()
    for freq, count in freq_counts.items():
        print(f"  Frequency {freq}: {count} introners")
    
    return binary_df, polymorphic_df, group1_samples

def calculate_combinatorial_complexity():
    """
    Calculate theoretical combinatorial complexity for each frequency
    """
    print("\nCalculating combinatorial complexity...")
    
    complexity_data = []
    max_complexity = math.comb(11, 5)  # Maximum at frequency 5
    
    for freq in range(2, 11):
        combinations = math.comb(11, freq)
        log10_combinations = math.log10(combinations)
        relative_complexity = combinations / max_complexity
        
        complexity_data.append({
            'frequency': freq,
            'samples_present': freq,
            'combinations': combinations,
            'log10_combinations': log10_combinations,
            'relative_complexity': relative_complexity
        })
    
    complexity_df = pd.DataFrame(complexity_data)
    
    print("\nCombinatorial complexity table:")
    print(complexity_df.round(3))
    
    # Save to file
    complexity_df.to_csv('introner_analysis/data_exploration/combinatorial_summary.csv', index=False)
    
    return complexity_df

def analyze_frequency_distribution(binary_df, polymorphic_df):
    """
    Analyze actual vs theoretical introner distribution across frequencies
    """
    print("\nAnalyzing frequency distribution...")
    
    # Get actual counts
    freq_counts = polymorphic_df['frequency'].value_counts().sort_index()
    
    # Create summary
    distribution_data = []
    for freq in range(2, 11):
        actual_count = freq_counts.get(freq, 0)
        theoretical_max = math.comb(11, freq)
        
        distribution_data.append({
            'frequency': freq,
            'actual_introners': actual_count,
            'theoretical_combinations': theoretical_max,
            'utilization_rate': actual_count / theoretical_max if theoretical_max > 0 else 0
        })
    
    distribution_df = pd.DataFrame(distribution_data)
    
    print("Frequency distribution analysis:")
    print(distribution_df.round(4))
    
    # Save results
    distribution_df.to_csv('introner_analysis/data_exploration/frequency_distribution.csv', index=False)
    freq_counts.to_csv('introner_analysis/data_exploration/introner_counts_by_freq.csv', header=['count'])
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot actual counts
    ax1.bar(distribution_df['frequency'], distribution_df['actual_introners'], alpha=0.7, color='steelblue')
    ax1.set_xlabel('Introner Frequency')
    ax1.set_ylabel('Number of Introners')
    ax1.set_title('Actual Introner Distribution by Frequency')
    ax1.set_xticks(range(2, 11))
    
    # Plot utilization rates
    ax2.bar(distribution_df['frequency'], distribution_df['utilization_rate'], alpha=0.7, color='darkgreen')
    ax2.set_xlabel('Introner Frequency')
    ax2.set_ylabel('Utilization Rate')
    ax2.set_title('Utilization of Theoretical Combinations')
    ax2.set_xticks(range(2, 11))
    ax2.set_ylim(0, max(distribution_df['utilization_rate']) * 1.1)
    
    plt.tight_layout()
    plt.savefig('introner_analysis/data_exploration/frequency_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return distribution_df

def filter_by_frequency_bins(binary_df, polymorphic_df):
    """
    Create frequency-specific datasets for analysis
    """
    print("\nFiltering introners by frequency bins...")
    
    frequency_datasets = {}
    
    for freq in range(2, 11):
        # Get introners at this frequency
        freq_introners = polymorphic_df[polymorphic_df['frequency'] == freq]['ortholog_id']
        freq_binary = binary_df.loc[freq_introners]
        
        frequency_datasets[freq] = {
            'binary_matrix': freq_binary,
            'introner_count': len(freq_introners),
            'introner_ids': freq_introners.tolist()
        }
        
        print(f"  Frequency {freq}: {len(freq_introners)} introners")
    
    return frequency_datasets

def main():
    """Main analysis pipeline"""
    print("Starting Introner Population Structure Analysis")
    print("=" * 50)
    
    # Create directory structure
    create_output_directories()
    
    # Phase 1: Data preparation
    binary_df, polymorphic_df, group1_samples = load_and_prepare_data()
    
    # Calculate theoretical complexity
    complexity_df = calculate_combinatorial_complexity()
    
    # Analyze frequency distribution
    distribution_df = analyze_frequency_distribution(binary_df, polymorphic_df)
    
    # Filter by frequency bins
    frequency_datasets = filter_by_frequency_bins(binary_df, polymorphic_df)
    
    print(f"\n✓ Phase 1 completed successfully")
    print(f"  - Created binary matrix: {binary_df.shape[0]} × {binary_df.shape[1]}")
    print(f"  - Analyzed {len(frequency_datasets)} frequency bins")
    print(f"  - Generated theoretical complexity analysis")
    
    return binary_df, polymorphic_df, group1_samples, frequency_datasets

if __name__ == "__main__":
    result = main()
    if result is not None:
        print(f"\n✓ Data preparation phase completed successfully")
    else:
        print("\n✗ Data preparation phase failed")
        sys.exit(1)
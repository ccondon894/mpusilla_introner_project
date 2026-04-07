#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
import os
import sys

def load_introner_clustering_results():
    """
    Load the population structure results from introner analysis
    """
    print("Loading introner-based population structure results...")
    
    # Load sample consistency matrix (shows which samples cluster together)
    consistency_df = pd.read_csv('introner_analysis/comparative_analysis/sample_relationship_matrix.csv', 
                                index_col=0)
    
    # Load cross-frequency summary
    cross_freq_df = pd.read_csv('introner_analysis/comparative_analysis/cross_frequency_summary.csv')
    
    # Load clustering results for frequency 4 (good balance of data and resolution)
    freq4_clusters = pd.read_csv('introner_analysis/clustering_by_frequency/freq_4/cluster_assignments.csv')
    
    print(f"  - Consistency matrix: {consistency_df.shape}")
    print(f"  - Using frequency 4 clustering (hierarchical k=2) as reference")
    
    return consistency_df, cross_freq_df, freq4_clusters

def load_demographic_dn_ds_data():
    """
    Load existing dN/dS analysis results if available
    """
    print("Loading demographic and dN/dS analysis data...")
    
    demographic_data = {}
    
    # Try to load dN/dS results
    try:
        dn_ds_results = pd.read_csv('/scratch1/chris/introner_vis/dn_ds_analysis_results.tsv', sep='\t')
        demographic_data['dn_ds'] = dn_ds_results
        print(f"  - Loaded dN/dS results: {len(dn_ds_results)} ortholog groups")
    except FileNotFoundError:
        print("  - dN/dS results not found, will skip dN/dS integration")
        demographic_data['dn_ds'] = None
    
    # Try to load selection dataset
    try:
        selection_data = pd.read_csv('/scratch1/chris/introner_vis/final_selection_dataset.tsv', sep='\t')
        demographic_data['selection'] = selection_data
        print(f"  - Loaded selection dataset: {len(selection_data)} records")
    except FileNotFoundError:
        print("  - Selection dataset not found, will skip selection integration")
        demographic_data['selection'] = None
    
    # Try to load polymorphic groups data
    try:
        polymorphic_data = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t')
        # Add frequency information
        frequencies = []
        for _, row in polymorphic_data.iterrows():
            present_samples = eval(row['present_samples']) if isinstance(row['present_samples'], str) else row['present_samples']
            frequencies.append(len(present_samples))
        polymorphic_data['frequency'] = frequencies
        demographic_data['polymorphic'] = polymorphic_data
        print(f"  - Loaded polymorphic groups: {len(polymorphic_data)} groups")
    except FileNotFoundError:
        print("  - Polymorphic groups data not found")
        demographic_data['polymorphic'] = None
    
    return demographic_data

def analyze_population_groups(consistency_df, freq4_clusters):
    """
    Define population groups based on introner clustering
    """
    print("\nAnalyzing population groupings...")
    
    # Use hierarchical k=2 clustering from frequency 4 as the reference
    samples = freq4_clusters['sample'].tolist()
    clusters = freq4_clusters['hierarchical_k_2'].tolist()
    
    # Create population groups
    population_groups = {}
    for sample, cluster in zip(samples, clusters):
        if cluster not in population_groups:
            population_groups[cluster] = []
        population_groups[cluster].append(sample)
    
    print(f"Population groups based on introner clustering:")
    for group_id, group_samples in population_groups.items():
        print(f"  Group {group_id}: {group_samples} (n={len(group_samples)})")
    
    # Analyze consistency within and between groups
    within_group_consistency = []
    between_group_consistency = []
    
    for group_id, group_samples in population_groups.items():
        # Within-group consistency
        for i, sample1 in enumerate(group_samples):
            for j, sample2 in enumerate(group_samples):
                if i < j:  # Avoid duplicates
                    consistency = consistency_df.loc[sample1, sample2]
                    within_group_consistency.append(consistency)
    
    # Between-group consistency
    group_ids = list(population_groups.keys())
    for i, group1 in enumerate(group_ids):
        for j, group2 in enumerate(group_ids):
            if i < j:
                for sample1 in population_groups[group1]:
                    for sample2 in population_groups[group2]:
                        consistency = consistency_df.loc[sample1, sample2]
                        between_group_consistency.append(consistency)
    
    print(f"\nConsistency analysis:")
    print(f"  Within-group mean consistency: {np.mean(within_group_consistency):.3f}")
    print(f"  Between-group mean consistency: {np.mean(between_group_consistency):.3f}")
    
    return population_groups, within_group_consistency, between_group_consistency

def analyze_frequency_spectrum_by_population(demographic_data, population_groups):
    """
    Analyze introner frequency spectrum by population groups
    """
    print("\nAnalyzing frequency spectrum by population groups...")
    
    if demographic_data['polymorphic'] is None:
        print("  Polymorphic data not available, skipping frequency analysis")
        return None
    
    polymorphic_df = demographic_data['polymorphic']
    
    # Analyze frequency distribution for each population group
    frequency_analysis = {}
    
    for group_id, group_samples in population_groups.items():
        group_frequencies = []
        
        for _, row in polymorphic_df.iterrows():
            present_samples = eval(row['present_samples']) if isinstance(row['present_samples'], str) else row['present_samples']
            
            # Count how many samples from this group have the introner
            group_present = sum(1 for sample in present_samples if sample in group_samples)
            group_freq = group_present / len(group_samples)
            group_frequencies.append(group_freq)
        
        frequency_analysis[group_id] = {
            'frequencies': group_frequencies,
            'mean_frequency': np.mean(group_frequencies),
            'std_frequency': np.std(group_frequencies)
        }
        
        print(f"  Group {group_id}: mean frequency {np.mean(group_frequencies):.3f} ± {np.std(group_frequencies):.3f}")
    
    return frequency_analysis

def analyze_dn_ds_by_population(demographic_data, population_groups):
    """
    Analyze dN/dS patterns by population groups if data is available
    """
    print("\nAnalyzing dN/dS patterns by population groups...")
    
    if demographic_data['dn_ds'] is None:
        print("  dN/dS data not available, skipping dN/dS analysis")
        return None
    
    dn_ds_df = demographic_data['dn_ds']
    
    # This is a placeholder - would need to map ortholog groups to population-specific patterns
    # The current dN/dS analysis compares introner-present vs absent across all samples
    
    print("  dN/dS analysis integration would require re-running dN/dS with population-specific groups")
    print("  Current dN/dS results are population-averaged")
    
    return {
        'note': 'dN/dS integration requires population-specific re-analysis',
        'current_results': dn_ds_df.describe()
    }

def test_population_structure_significance(population_groups, consistency_df):
    """
    Test statistical significance of population structure
    """
    print("\nTesting population structure significance...")
    
    # Perform permutation test
    n_permutations = 1000
    observed_difference = []
    
    # Calculate observed within vs between group differences
    within_scores = []
    between_scores = []
    
    group_ids = list(population_groups.keys())
    
    # Within-group scores
    for group_id, group_samples in population_groups.items():
        for i, sample1 in enumerate(group_samples):
            for j, sample2 in enumerate(group_samples):
                if i < j:
                    within_scores.append(consistency_df.loc[sample1, sample2])
    
    # Between-group scores
    for i, group1 in enumerate(group_ids):
        for j, group2 in enumerate(group_ids):
            if i < j:
                for sample1 in population_groups[group1]:
                    for sample2 in population_groups[group2]:
                        between_scores.append(consistency_df.loc[sample1, sample2])
    
    observed_diff = np.mean(within_scores) - np.mean(between_scores)
    
    # Permutation test
    permutation_diffs = []
    all_samples = [sample for group in population_groups.values() for sample in group]
    
    for perm in range(n_permutations):
        # Randomly reassign samples to groups
        shuffled_samples = np.random.permutation(all_samples)
        perm_groups = {}
        idx = 0
        for group_id, group_samples in population_groups.items():
            perm_groups[group_id] = shuffled_samples[idx:idx+len(group_samples)].tolist()
            idx += len(group_samples)
        
        # Calculate permuted within vs between scores
        perm_within = []
        perm_between = []
        
        # Within-group scores
        for group_id, group_samples in perm_groups.items():
            for i, sample1 in enumerate(group_samples):
                for j, sample2 in enumerate(group_samples):
                    if i < j:
                        perm_within.append(consistency_df.loc[sample1, sample2])
        
        # Between-group scores  
        for i, group1 in enumerate(group_ids):
            for j, group2 in enumerate(group_ids):
                if i < j:
                    for sample1 in perm_groups[group1]:
                        for sample2 in perm_groups[group2]:
                            perm_between.append(consistency_df.loc[sample1, sample2])
        
        perm_diff = np.mean(perm_within) - np.mean(perm_between)
        permutation_diffs.append(perm_diff)
    
    # Calculate p-value
    p_value = np.sum(np.array(permutation_diffs) >= observed_diff) / n_permutations
    
    print(f"Population structure significance test:")
    print(f"  Observed within-between difference: {observed_diff:.4f}")
    print(f"  Permutation test p-value: {p_value:.4f}")
    print(f"  Mean within-group consistency: {np.mean(within_scores):.4f}")
    print(f"  Mean between-group consistency: {np.mean(between_scores):.4f}")
    
    return {
        'observed_difference': observed_diff,
        'p_value': p_value,
        'within_group_mean': np.mean(within_scores),
        'between_group_mean': np.mean(between_scores),
        'permutation_differences': permutation_diffs
    }

def create_integration_visualizations(population_groups, consistency_df, 
                                    frequency_analysis, significance_test, output_dir):
    """
    Create visualizations integrating population structure with demographic patterns
    """
    print("Creating integration visualizations...")
    
    # 1. Population groups consistency heatmap
    plt.figure(figsize=(10, 8))
    
    # Reorder samples by population groups
    ordered_samples = []
    for group_id in sorted(population_groups.keys()):
        ordered_samples.extend(population_groups[group_id])
    
    ordered_consistency = consistency_df.loc[ordered_samples, ordered_samples]
    
    sns.heatmap(ordered_consistency, annot=True, fmt='.2f', cmap='YlOrRd',
                square=True, cbar_kws={'label': 'Consistency Score'})
    
    plt.title('Sample Consistency by Population Groups')
    plt.xlabel('Samples (ordered by population group)')
    plt.ylabel('Samples (ordered by population group)')
    
    # Add group boundaries
    cumulative_sizes = [0]
    for group_id in sorted(population_groups.keys()):
        cumulative_sizes.append(cumulative_sizes[-1] + len(population_groups[group_id]))
    
    for boundary in cumulative_sizes[1:-1]:
        plt.axhline(y=boundary, color='white', linewidth=2)
        plt.axvline(x=boundary, color='white', linewidth=2)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/population_groups_consistency.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Significance test visualization
    if significance_test:
        plt.figure(figsize=(10, 6))
        
        # Histogram of permutation differences
        plt.hist(significance_test['permutation_differences'], bins=50, alpha=0.7, 
                density=True, color='lightgray', label='Permuted differences')
        
        # Mark observed difference
        plt.axvline(significance_test['observed_difference'], color='red', 
                   linewidth=2, label=f'Observed difference: {significance_test["observed_difference"]:.4f}')
        
        plt.xlabel('Within-group - Between-group Consistency')
        plt.ylabel('Density')
        plt.title(f'Population Structure Significance Test\n(p = {significance_test["p_value"]:.4f})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/population_significance_test.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Frequency analysis by group
    if frequency_analysis:
        plt.figure(figsize=(12, 6))
        
        # Box plots of frequency distributions by group
        group_data = []
        group_labels = []
        
        for group_id in sorted(frequency_analysis.keys()):
            group_data.append(frequency_analysis[group_id]['frequencies'])
            group_labels.append(f'Group {group_id}')
        
        plt.boxplot(group_data, labels=group_labels)
        plt.ylabel('Introner Frequency within Group')
        plt.title('Introner Frequency Distribution by Population Group')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/frequency_by_population.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    """Main Phase 4 analysis pipeline"""
    print("Starting Phase 4: Integration with Demographic Patterns")
    print("=" * 56)
    
    output_dir = 'introner_analysis/demographic_integration'
    
    # Load introner clustering results
    consistency_df, cross_freq_df, freq4_clusters = load_introner_clustering_results()
    
    # Load demographic/dN/dS data
    demographic_data = load_demographic_dn_ds_data()
    
    # Analyze population groups
    population_groups, within_consistency, between_consistency = analyze_population_groups(
        consistency_df, freq4_clusters)
    
    # Analyze frequency spectrum by population
    frequency_analysis = analyze_frequency_spectrum_by_population(demographic_data, population_groups)
    
    # Analyze dN/dS by population
    dn_ds_analysis = analyze_dn_ds_by_population(demographic_data, population_groups)
    
    # Test significance of population structure
    significance_test = test_population_structure_significance(population_groups, consistency_df)
    
    # Create visualizations
    create_integration_visualizations(population_groups, consistency_df, 
                                    frequency_analysis, significance_test, output_dir)
    
    # Save results
    # Population groups
    population_df = pd.DataFrame([
        {'group': group_id, 'sample': sample}
        for group_id, samples in population_groups.items()
        for sample in samples
    ])
    population_df.to_csv(f'{output_dir}/population_groups.csv', index=False)
    
    # Integration summary
    integration_summary = {
        'n_population_groups': len(population_groups),
        'within_group_consistency': np.mean(within_consistency),
        'between_group_consistency': np.mean(between_consistency),
        'structure_significance_p': significance_test['p_value'] if significance_test else None,
        'frequency_analysis_available': frequency_analysis is not None,
        'dn_ds_analysis_available': dn_ds_analysis is not None
    }
    
    integration_df = pd.DataFrame([integration_summary])
    integration_df.to_csv(f'{output_dir}/integration_summary.csv', index=False)
    
    print(f"\n✓ Phase 4 completed successfully")
    print(f"  - Identified {len(population_groups)} population groups")
    print(f"  - Within-group consistency: {np.mean(within_consistency):.3f}")
    print(f"  - Between-group consistency: {np.mean(between_consistency):.3f}")
    if significance_test:
        print(f"  - Population structure p-value: {significance_test['p_value']:.4f}")
    
    return {
        'population_groups': population_groups,
        'consistency_analysis': (within_consistency, between_consistency),
        'frequency_analysis': frequency_analysis,
        'dn_ds_analysis': dn_ds_analysis,
        'significance_test': significance_test,
        'integration_summary': integration_summary
    }

if __name__ == "__main__":
    results = main()
    if results is not None:
        print(f"\n✓ Demographic integration analysis completed successfully")
    else:
        print("\n✗ Demographic integration analysis failed")
        sys.exit(1)
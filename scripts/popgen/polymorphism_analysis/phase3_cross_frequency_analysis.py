#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import os
import sys

def load_frequency_results():
    """
    Load clustering results from all frequency bins
    """
    print("Loading frequency-specific clustering results...")
    
    frequency_results = {}
    
    for freq in range(2, 11):
        freq_dir = f'introner_analysis/clustering_by_frequency/freq_{freq}'
        
        if not os.path.exists(freq_dir):
            print(f"Warning: No results found for frequency {freq}")
            continue
            
        try:
            # Load cluster assignments
            clusters_df = pd.read_csv(f'{freq_dir}/cluster_assignments.csv')
            
            # Load distance matrix
            distance_df = pd.read_csv(f'{freq_dir}/distance_matrix.csv', index_col=0)
            
            # Load metrics
            metrics_df = pd.read_csv(f'{freq_dir}/clustering_metrics.csv')
            
            # Load PCA coordinates
            pca_df = pd.read_csv(f'{freq_dir}/pca_coordinates.csv', index_col=0)
            
            frequency_results[freq] = {
                'clusters': clusters_df,
                'distance_matrix': distance_df,
                'metrics': metrics_df,
                'pca_coordinates': pca_df
            }
            
            print(f"  Loaded results for frequency {freq}")
            
        except Exception as e:
            print(f"Warning: Could not load results for frequency {freq}: {e}")
    
    return frequency_results

def calculate_distance_matrix_correlations(frequency_results):
    """
    Calculate correlations between distance matrices across frequencies
    """
    print("\nCalculating distance matrix correlations...")
    
    frequencies = sorted(frequency_results.keys())
    n_freq = len(frequencies)
    
    # Initialize correlation matrices
    pearson_corr = np.zeros((n_freq, n_freq))
    spearman_corr = np.zeros((n_freq, n_freq))
    
    for i, freq1 in enumerate(frequencies):
        for j, freq2 in enumerate(frequencies):
            if i <= j:  # Only calculate upper triangle
                dist1 = frequency_results[freq1]['distance_matrix'].values
                dist2 = frequency_results[freq2]['distance_matrix'].values
                
                # Get upper triangular elements (excluding diagonal)
                triu_indices = np.triu_indices_from(dist1, k=1)
                dist1_flat = dist1[triu_indices]
                dist2_flat = dist2[triu_indices]
                
                # Calculate correlations
                pearson_r, _ = pearsonr(dist1_flat, dist2_flat)
                spearman_r, _ = spearmanr(dist1_flat, dist2_flat)
                
                pearson_corr[i, j] = pearson_r
                pearson_corr[j, i] = pearson_r  # Mirror
                spearman_corr[i, j] = spearman_r
                spearman_corr[j, i] = spearman_r  # Mirror
    
    # Create correlation DataFrames
    freq_labels = [f'Freq_{freq}' for freq in frequencies]
    pearson_df = pd.DataFrame(pearson_corr, index=freq_labels, columns=freq_labels)
    spearman_df = pd.DataFrame(spearman_corr, index=freq_labels, columns=freq_labels)
    
    return pearson_df, spearman_df, frequencies

def calculate_clustering_concordance(frequency_results):
    """
    Calculate clustering concordance across frequencies using ARI and NMI
    """
    print("Calculating clustering concordance...")
    
    frequencies = sorted(frequency_results.keys())
    samples = frequency_results[frequencies[0]]['clusters']['sample'].tolist()
    
    # For each clustering method and k value
    clustering_methods = ['hierarchical_k_2', 'hierarchical_k_3', 'hierarchical_k_4', 
                         'kmeans_k2', 'kmeans_k3', 'kmeans_k4']
    
    concordance_results = []
    
    for method in clustering_methods:
        # Check if this method exists across frequencies
        method_exists = all(method in frequency_results[freq]['clusters'].columns 
                           for freq in frequencies)
        
        if not method_exists:
            print(f"Warning: Method {method} not available in all frequencies, skipping")
            continue
        
        # Calculate pairwise concordance
        n_freq = len(frequencies)
        ari_matrix = np.zeros((n_freq, n_freq))
        nmi_matrix = np.zeros((n_freq, n_freq))
        
        for i, freq1 in enumerate(frequencies):
            for j, freq2 in enumerate(frequencies):
                if i <= j:
                    clusters1 = frequency_results[freq1]['clusters'][method].values
                    clusters2 = frequency_results[freq2]['clusters'][method].values
                    
                    # Calculate concordance metrics
                    ari = adjusted_rand_score(clusters1, clusters2)
                    nmi = normalized_mutual_info_score(clusters1, clusters2)
                    
                    ari_matrix[i, j] = ari
                    ari_matrix[j, i] = ari
                    nmi_matrix[i, j] = nmi
                    nmi_matrix[j, i] = nmi
        
        concordance_results.append({
            'method': method,
            'ari_matrix': ari_matrix,
            'nmi_matrix': nmi_matrix,
            'mean_ari': np.mean(ari_matrix[np.triu_indices_from(ari_matrix, k=1)]),
            'mean_nmi': np.mean(nmi_matrix[np.triu_indices_from(nmi_matrix, k=1)])
        })
    
    return concordance_results, frequencies

def analyze_sample_relationship_consistency(frequency_results):
    """
    Analyze which sample pairs consistently cluster together across frequencies
    """
    print("Analyzing sample relationship consistency...")
    
    frequencies = sorted(frequency_results.keys())
    samples = frequency_results[frequencies[0]]['clusters']['sample'].tolist()
    n_samples = len(samples)
    
    # Initialize consistency matrix
    consistency_matrix = np.zeros((n_samples, n_samples))
    
    # For hierarchical clustering with k=2 (most robust)
    method = 'hierarchical_k_2'
    
    for freq in frequencies:
        if method in frequency_results[freq]['clusters'].columns:
            clusters = frequency_results[freq]['clusters'][method].values
            
            # For each pair of samples
            for i in range(n_samples):
                for j in range(i+1, n_samples):
                    # Check if they're in the same cluster
                    if clusters[i] == clusters[j]:
                        consistency_matrix[i, j] += 1
                        consistency_matrix[j, i] += 1
    
    # Normalize by number of frequencies
    n_frequencies = len(frequencies)
    consistency_matrix = consistency_matrix / n_frequencies
    
    # Convert to DataFrame
    consistency_df = pd.DataFrame(consistency_matrix, 
                                 index=samples, 
                                 columns=samples)
    
    return consistency_df

def create_correlation_heatmaps(pearson_df, spearman_df, output_dir):
    """
    Create heatmaps showing distance matrix correlations
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Pearson correlation heatmap
    sns.heatmap(pearson_df, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                square=True, ax=ax1, cbar_kws={'label': 'Pearson Correlation'})
    ax1.set_title('Distance Matrix Correlations\n(Pearson)')
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('Frequency')
    
    # Spearman correlation heatmap
    sns.heatmap(spearman_df, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                square=True, ax=ax2, cbar_kws={'label': 'Spearman Correlation'})
    ax2.set_title('Distance Matrix Correlations\n(Spearman)')
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/distance_correlations.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_concordance_summary(concordance_results, frequencies, output_dir):
    """
    Create summary plots of clustering concordance
    """
    # Prepare data for plotting
    methods = [result['method'] for result in concordance_results]
    mean_aris = [result['mean_ari'] for result in concordance_results]
    mean_nmis = [result['mean_nmi'] for result in concordance_results]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # ARI plot
    bars1 = ax1.bar(methods, mean_aris, alpha=0.7, color='steelblue')
    ax1.set_ylabel('Mean Adjusted Rand Index')
    ax1.set_title('Clustering Concordance Across Frequencies\n(Adjusted Rand Index)')
    ax1.set_xticklabels(methods, rotation=45, ha='right')
    ax1.set_ylim(0, max(mean_aris) * 1.1)
    
    # Add value labels on bars
    for bar, value in zip(bars1, mean_aris):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.3f}', ha='center', va='bottom')
    
    # NMI plot
    bars2 = ax2.bar(methods, mean_nmis, alpha=0.7, color='darkgreen')
    ax2.set_ylabel('Mean Normalized Mutual Information')
    ax2.set_title('Clustering Concordance Across Frequencies\n(Normalized Mutual Information)')
    ax2.set_xticklabels(methods, rotation=45, ha='right')
    ax2.set_ylim(0, max(mean_nmis) * 1.1)
    
    # Add value labels on bars
    for bar, value in zip(bars2, mean_nmis):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.3f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/clustering_concordance_summary.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_consistency_heatmap(consistency_df, output_dir):
    """
    Create heatmap showing sample relationship consistency
    """
    plt.figure(figsize=(10, 8))
    
    sns.heatmap(consistency_df, annot=True, fmt='.2f', cmap='YlOrRd',
                square=True, cbar_kws={'label': 'Consistency Score'})
    
    plt.title('Sample Relationship Consistency Across Frequencies\n(Proportion of frequencies where samples cluster together)')
    plt.xlabel('Samples')
    plt.ylabel('Samples')
    plt.tight_layout()
    
    plt.savefig(f'{output_dir}/sample_consistency.png', dpi=300, bbox_inches='tight')
    plt.close()

def compile_summary_statistics(frequency_results, pearson_df, spearman_df, 
                              concordance_results, consistency_df):
    """
    Compile comprehensive summary statistics
    """
    print("Compiling summary statistics...")
    
    # Collect metrics from each frequency
    summary_data = []
    
    for freq, results in frequency_results.items():
        metrics = results['metrics'].iloc[0].to_dict()
        summary_data.append(metrics)
    
    summary_df = pd.DataFrame(summary_data)
    
    # Add cross-frequency statistics
    cross_freq_stats = {
        'mean_distance_pearson_corr': np.mean(pearson_df.values[np.triu_indices_from(pearson_df.values, k=1)]),
        'mean_distance_spearman_corr': np.mean(spearman_df.values[np.triu_indices_from(spearman_df.values, k=1)]),
        'max_distance_correlation': np.max(pearson_df.values[np.triu_indices_from(pearson_df.values, k=1)]),
        'min_distance_correlation': np.min(pearson_df.values[np.triu_indices_from(pearson_df.values, k=1)])
    }
    
    # Add concordance statistics
    if concordance_results:
        best_ari = max(result['mean_ari'] for result in concordance_results)
        best_nmi = max(result['mean_nmi'] for result in concordance_results)
        cross_freq_stats['best_clustering_ari'] = best_ari
        cross_freq_stats['best_clustering_nmi'] = best_nmi
    
    # Add consistency statistics
    cross_freq_stats['mean_sample_consistency'] = np.mean(consistency_df.values[np.triu_indices_from(consistency_df.values, k=1)])
    cross_freq_stats['max_sample_consistency'] = np.max(consistency_df.values[np.triu_indices_from(consistency_df.values, k=1)])
    
    return summary_df, cross_freq_stats

def main():
    """Main Phase 3 analysis pipeline"""
    print("Starting Phase 3: Cross-Frequency Pattern Analysis")
    print("=" * 55)
    
    output_dir = 'introner_analysis/comparative_analysis'
    
    # Load frequency-specific results
    frequency_results = load_frequency_results()
    
    if not frequency_results:
        print("Error: No frequency results found. Run Phase 2 first.")
        return None
    
    print(f"Loaded results for {len(frequency_results)} frequency bins")
    
    # Calculate distance matrix correlations
    pearson_df, spearman_df, frequencies = calculate_distance_matrix_correlations(frequency_results)
    
    # Save correlation matrices
    pearson_df.to_csv(f'{output_dir}/distance_pearson_correlations.csv')
    spearman_df.to_csv(f'{output_dir}/distance_spearman_correlations.csv')
    
    # Calculate clustering concordance
    concordance_results, frequencies = calculate_clustering_concordance(frequency_results)
    
    # Save concordance results
    concordance_summary = pd.DataFrame([{
        'method': result['method'],
        'mean_ari': result['mean_ari'],
        'mean_nmi': result['mean_nmi']
    } for result in concordance_results])
    concordance_summary.to_csv(f'{output_dir}/clustering_concordance.csv', index=False)
    
    # Analyze sample relationship consistency
    consistency_df = analyze_sample_relationship_consistency(frequency_results)
    consistency_df.to_csv(f'{output_dir}/sample_relationship_matrix.csv')
    
    # Create visualizations
    create_correlation_heatmaps(pearson_df, spearman_df, output_dir)
    create_concordance_summary(concordance_results, frequencies, output_dir)
    create_consistency_heatmap(consistency_df, output_dir)
    
    # Compile comprehensive summary
    summary_df, cross_freq_stats = compile_summary_statistics(
        frequency_results, pearson_df, spearman_df, concordance_results, consistency_df)
    
    # Save summary statistics
    summary_df.to_csv(f'{output_dir}/frequency_specific_summary.csv', index=False)
    
    cross_freq_df = pd.DataFrame([cross_freq_stats])
    cross_freq_df.to_csv(f'{output_dir}/cross_frequency_summary.csv', index=False)
    
    print(f"\n✓ Phase 3 completed successfully")
    print(f"  - Analyzed correlations across {len(frequencies)} frequency bins")
    print(f"  - Calculated clustering concordance for {len(concordance_results)} methods")
    print(f"  - Generated sample relationship consistency matrix")
    print(f"  - Created comprehensive visualizations and summaries")
    
    return {
        'frequency_results': frequency_results,
        'correlations': {'pearson': pearson_df, 'spearman': spearman_df},
        'concordance': concordance_results,
        'consistency': consistency_df,
        'summary': summary_df,
        'cross_freq_stats': cross_freq_stats
    }

if __name__ == "__main__":
    results = main()
    if results is not None:
        print(f"\n✓ Cross-frequency pattern analysis completed successfully")
    else:
        print("\n✗ Cross-frequency pattern analysis failed")
        sys.exit(1)
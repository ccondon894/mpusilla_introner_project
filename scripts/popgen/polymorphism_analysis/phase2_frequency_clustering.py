#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform, jaccard
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os
import sys

def calculate_jaccard_distance_matrix(binary_matrix):
    """
    Calculate Jaccard distance matrix between samples
    """
    # Transpose to get samples as rows
    samples_matrix = binary_matrix.T
    
    # Calculate pairwise Jaccard distances
    distances = pdist(samples_matrix, metric='jaccard')
    distance_matrix = squareform(distances)
    
    return distance_matrix

def perform_hierarchical_clustering(distance_matrix, sample_names, freq, output_dir):
    """
    Perform hierarchical clustering and create dendrogram
    """
    # Perform linkage
    linkage_matrix = linkage(distance_matrix, method='ward')
    
    # Create dendrogram
    plt.figure(figsize=(10, 6))
    dendrogram(linkage_matrix, labels=sample_names, orientation='top', 
               leaf_rotation=45, leaf_font_size=10)
    plt.title(f'Hierarchical Clustering - Frequency {freq} ({len(sample_names)} samples)')
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_dir}/dendrogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Get cluster assignments for different k values
    cluster_assignments = {}
    for k in [2, 3, 4]:
        if len(sample_names) >= k:
            clusters = fcluster(linkage_matrix, k, criterion='maxclust')
            cluster_assignments[f'k_{k}'] = clusters
    
    return linkage_matrix, cluster_assignments

def perform_pca_analysis(binary_matrix, sample_names, freq, output_dir):
    """
    Perform PCA analysis and create scatter plots
    """
    # Transpose to get samples as rows
    samples_matrix = binary_matrix.T.values
    
    # Perform PCA
    pca = PCA(n_components=min(len(sample_names)-1, binary_matrix.shape[0]))
    pca_result = pca.fit_transform(samples_matrix)
    
    # Calculate explained variance
    explained_variance = pca.explained_variance_ratio_
    
    # Create PCA plot
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], 
                         alpha=0.7, s=100, c=range(len(sample_names)), 
                         cmap='tab10')
    
    # Add sample labels
    for i, sample in enumerate(sample_names):
        plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    plt.xlabel(f'PC1 ({explained_variance[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({explained_variance[1]:.1%} variance)')
    plt.title(f'PCA Analysis - Frequency {freq} ({binary_matrix.shape[0]} introners)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_dir}/pca_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save PCA results
    pca_df = pd.DataFrame(pca_result[:, :min(4, pca_result.shape[1])], 
                         index=sample_names,
                         columns=[f'PC{i+1}' for i in range(min(4, pca_result.shape[1]))])
    pca_df.to_csv(f'{output_dir}/pca_coordinates.csv')
    
    # Save explained variance
    variance_df = pd.DataFrame({
        'component': [f'PC{i+1}' for i in range(len(explained_variance))],
        'explained_variance_ratio': explained_variance,
        'cumulative_variance': np.cumsum(explained_variance)
    })
    variance_df.to_csv(f'{output_dir}/pca_explained_variance.csv', index=False)
    
    return pca_result, explained_variance

def perform_kmeans_clustering(binary_matrix, sample_names, freq, output_dir):
    """
    Perform K-means clustering with different k values
    """
    # Transpose to get samples as rows
    samples_matrix = binary_matrix.T.values
    
    kmeans_results = {}
    silhouette_scores = {}
    
    # Try different k values
    for k in [2, 3, 4]:
        if len(sample_names) >= k:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(samples_matrix)
            
            # Calculate silhouette score
            if k > 1:
                sil_score = silhouette_score(samples_matrix, cluster_labels)
                silhouette_scores[k] = sil_score
            
            kmeans_results[k] = {
                'labels': cluster_labels,
                'centers': kmeans.cluster_centers_,
                'inertia': kmeans.inertia_
            }
    
    return kmeans_results, silhouette_scores

def create_distance_heatmap(distance_matrix, sample_names, freq, output_dir):
    """
    Create distance heatmap between samples
    """
    plt.figure(figsize=(8, 6))
    
    # Create heatmap
    sns.heatmap(distance_matrix, 
                xticklabels=sample_names, 
                yticklabels=sample_names,
                annot=True, fmt='.3f', cmap='viridis',
                square=True, cbar_kws={'label': 'Jaccard Distance'})
    
    plt.title(f'Sample Distance Matrix - Frequency {freq}')
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_dir}/distance_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

def analyze_frequency_bin(freq, binary_matrix, sample_names):
    """
    Complete analysis pipeline for a single frequency bin
    """
    print(f"\n--- Analyzing Frequency {freq} ---")
    print(f"Introners: {binary_matrix.shape[0]}, Samples: {binary_matrix.shape[1]}")
    
    if binary_matrix.shape[0] == 0:
        print(f"Skipping frequency {freq}: no introners")
        return None
    
    # Create output directory
    output_dir = f'introner_analysis/clustering_by_frequency/freq_{freq}'
    
    # Calculate distance matrix
    print("Calculating Jaccard distance matrix...")
    distance_matrix = calculate_jaccard_distance_matrix(binary_matrix)
    
    # Save distance matrix
    distance_df = pd.DataFrame(distance_matrix, 
                              index=sample_names, 
                              columns=sample_names)
    distance_df.to_csv(f'{output_dir}/distance_matrix.csv')
    
    # Create distance heatmap
    create_distance_heatmap(distance_matrix, sample_names, freq, output_dir)
    
    # Hierarchical clustering
    print("Performing hierarchical clustering...")
    linkage_matrix, hierarchical_clusters = perform_hierarchical_clustering(
        distance_matrix, sample_names, freq, output_dir)
    
    # PCA analysis
    print("Performing PCA analysis...")
    pca_result, explained_variance = perform_pca_analysis(
        binary_matrix, sample_names, freq, output_dir)
    
    # K-means clustering
    print("Performing K-means clustering...")
    kmeans_results, silhouette_scores = perform_kmeans_clustering(
        binary_matrix, sample_names, freq, output_dir)
    
    # Create cluster assignment table
    cluster_data = {'sample': sample_names}
    
    # Add hierarchical clusters
    for k, clusters in hierarchical_clusters.items():
        cluster_data[f'hierarchical_{k}'] = clusters
    
    # Add k-means clusters
    for k, result in kmeans_results.items():
        cluster_data[f'kmeans_k{k}'] = result['labels'] + 1  # +1 to match hierarchical numbering
    
    cluster_df = pd.DataFrame(cluster_data)
    cluster_df.to_csv(f'{output_dir}/cluster_assignments.csv', index=False)
    
    # Save clustering metrics
    metrics = {
        'frequency': freq,
        'n_introners': binary_matrix.shape[0],
        'n_samples': binary_matrix.shape[1],
        'pca_pc1_variance': explained_variance[0] if len(explained_variance) > 0 else 0,
        'pca_pc2_variance': explained_variance[1] if len(explained_variance) > 1 else 0,
        'mean_jaccard_distance': np.mean(distance_matrix[np.triu_indices_from(distance_matrix, k=1)]),
        'std_jaccard_distance': np.std(distance_matrix[np.triu_indices_from(distance_matrix, k=1)])
    }
    
    # Add silhouette scores
    for k, score in silhouette_scores.items():
        metrics[f'silhouette_k{k}'] = score
    
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(f'{output_dir}/clustering_metrics.csv', index=False)
    
    print(f"Analysis completed for frequency {freq}")
    return {
        'frequency': freq,
        'distance_matrix': distance_matrix,
        'cluster_assignments': cluster_df,
        'metrics': metrics,
        'pca_result': pca_result,
        'explained_variance': explained_variance
    }

def load_data_from_phase1():
    """
    Load the data prepared in Phase 1
    """
    print("Loading data from Phase 1...")
    
    # Load polymorphic groups data
    polymorphic_df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t')
    
    # Define Group1 samples
    group1_samples = ['CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
                      'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    # Recreate binary matrix (same as Phase 1)
    binary_matrix = []
    introner_ids = []
    
    for idx, row in polymorphic_df.iterrows():
        ortholog_id = row['ortholog_id']
        present_samples = eval(row['present_samples']) if isinstance(row['present_samples'], str) else row['present_samples']
        
        binary_row = []
        for sample in group1_samples:
            if sample in present_samples:
                binary_row.append(1)
            else:
                binary_row.append(0)
        
        binary_matrix.append(binary_row)
        introner_ids.append(ortholog_id)
    
    binary_df = pd.DataFrame(binary_matrix, 
                           columns=group1_samples, 
                           index=introner_ids)
    
    # Calculate frequencies and add to polymorphic_df
    frequencies = binary_df.sum(axis=1)
    polymorphic_df['frequency'] = frequencies.values
    
    return binary_df, polymorphic_df, group1_samples

def main():
    """Main Phase 2 analysis pipeline"""
    print("Starting Phase 2: Frequency-Stratified Clustering Analysis")
    print("=" * 60)
    
    # Load data from Phase 1
    binary_df, polymorphic_df, group1_samples = load_data_from_phase1()
    
    # Analyze each frequency bin
    results = {}
    
    for freq in range(2, 11):
        # Get introners at this frequency
        freq_introners = polymorphic_df[polymorphic_df['frequency'] == freq]['ortholog_id']
        
        if len(freq_introners) == 0:
            print(f"Skipping frequency {freq}: no introners")
            continue
            
        freq_binary = binary_df.loc[freq_introners]
        
        # Analyze this frequency bin
        result = analyze_frequency_bin(freq, freq_binary, group1_samples)
        if result is not None:
            results[freq] = result
    
    print(f"\n✓ Phase 2 completed successfully")
    print(f"  - Analyzed {len(results)} frequency bins")
    print(f"  - Generated clustering results for each bin")
    
    return results

if __name__ == "__main__":
    results = main()
    if results is not None:
        print(f"\n✓ Frequency-stratified clustering analysis completed successfully")
    else:
        print("\n✗ Frequency-stratified clustering analysis failed")
        sys.exit(1)
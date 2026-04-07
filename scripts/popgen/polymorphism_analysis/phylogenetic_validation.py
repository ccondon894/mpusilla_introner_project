#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from ete3 import Tree
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import squareform
import os
import sys

def load_introner_population_groups():
    """
    Load the introner-based population groups identified in previous analysis
    """
    print("Loading introner-based population groups...")
    
    # Define the groups identified from introner analysis
    introner_groups = {
        'Group1': ['RCC1614', 'RCC1698', 'RCC373', 'RCC465'],
        'Group2': ['CCMP1545', 'RCC114', 'RCC2482', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    }
    
    print(f"Group1: {introner_groups['Group1']} (n={len(introner_groups['Group1'])})")
    print(f"Group2: {introner_groups['Group2']} (n={len(introner_groups['Group2'])})")
    
    # Load introner distance matrix for comparison
    try:
        introner_distances = pd.read_csv(
            'introner_analysis/comparative_analysis/sample_relationship_matrix.csv', 
            index_col=0
        )
        # Convert consistency to distance (1 - consistency)
        introner_distances = 1 - introner_distances
        print(f"Loaded introner distance matrix: {introner_distances.shape}")
    except FileNotFoundError:
        print("Warning: Introner distance matrix not found")
        introner_distances = None
    
    return introner_groups, introner_distances

def parse_and_root_tree():
    """
    Parse the phylogenetic tree and root it with outgroups
    """
    print("Parsing and rooting phylogenetic tree...")
    
    tree_file = '/scratch1/chris/introner_vis/iqtree/mpusilla.snps.4d.notMT.min4.phy.treefile'
    
    # Load tree with ete3 for better manipulation
    tree = Tree(tree_file)
    
    # Define outgroups
    outgroups = ['RCC1749', 'RCC3052']
    
    # Find outgroup nodes
    outgroup_nodes = []
    for node in tree.traverse():
        if node.is_leaf() and node.name in outgroups:
            outgroup_nodes.append(node)
    
    if len(outgroup_nodes) == 2:
        # Find common ancestor of outgroups
        outgroup_ancestor = tree.get_common_ancestor(outgroup_nodes)
        # Root tree at outgroup ancestor
        tree.set_outgroup(outgroup_ancestor)
        print(f"Tree rooted using outgroups: {outgroups}")
    else:
        print(f"Warning: Could not find both outgroups. Found: {[n.name for n in outgroup_nodes]}")
    
    # Get all leaf names
    all_samples = [leaf.name for leaf in tree.get_leaves()]
    print(f"Tree contains {len(all_samples)} samples: {all_samples}")
    
    return tree, all_samples

def extract_group1_subtree(tree, introner_groups):
    """
    Extract subtree containing only Group1 samples for focused analysis
    """
    print("Extracting Group1 subtree...")
    
    # Get all Group1 samples
    group1_samples = introner_groups['Group1'] + introner_groups['Group2']
    
    # Find Group1 leaf nodes
    group1_nodes = []
    for node in tree.traverse():
        if node.is_leaf() and node.name in group1_samples:
            group1_nodes.append(node)
    
    print(f"Found {len(group1_nodes)} Group1 samples in tree: {[n.name for n in group1_nodes]}")
    
    if len(group1_nodes) < len(group1_samples):
        missing = set(group1_samples) - set([n.name for n in group1_nodes])
        print(f"Warning: Missing samples from tree: {missing}")
    
    # Get common ancestor of all Group1 samples
    if len(group1_nodes) > 1:
        group1_ancestor = tree.get_common_ancestor(group1_nodes)
        # Create subtree
        group1_subtree = group1_ancestor.copy()
        
        # Remove non-Group1 leaves from subtree
        leaves_to_remove = []
        for leaf in group1_subtree.get_leaves():
            if leaf.name not in group1_samples:
                leaves_to_remove.append(leaf)
        
        for leaf in leaves_to_remove:
            leaf.delete()
        
        print(f"Group1 subtree contains {len(group1_subtree.get_leaves())} samples")
        return group1_subtree
    else:
        print("Error: Not enough Group1 samples found in tree")
        return None

def calculate_patristic_distances(tree, sample_list):
    """
    Calculate patristic (tree-based) distances between all sample pairs
    """
    print("Calculating patristic distances...")
    
    # Get leaf nodes for samples
    sample_nodes = {}
    for node in tree.traverse():
        if node.is_leaf() and node.name in sample_list:
            sample_nodes[node.name] = node
    
    # Calculate pairwise distances
    n_samples = len(sample_list)
    distance_matrix = np.zeros((n_samples, n_samples))
    
    for i, sample1 in enumerate(sample_list):
        for j, sample2 in enumerate(sample_list):
            if i != j and sample1 in sample_nodes and sample2 in sample_nodes:
                # Calculate patristic distance
                distance = sample_nodes[sample1].get_distance(sample_nodes[sample2])
                distance_matrix[i, j] = distance
    
    # Convert to DataFrame
    distance_df = pd.DataFrame(distance_matrix, index=sample_list, columns=sample_list)
    
    print(f"Calculated patristic distances for {len(sample_list)} samples")
    return distance_df

def test_monophyly(tree, introner_groups):
    """
    Test if introner-defined groups are monophyletic in the phylogenetic tree
    """
    print("Testing monophyly of introner groups...")
    
    monophyly_results = {}
    
    for group_name, samples in introner_groups.items():
        # Find leaf nodes for this group
        group_nodes = []
        for node in tree.traverse():
            if node.is_leaf() and node.name in samples:
                group_nodes.append(node)
        
        if len(group_nodes) > 1:
            # Check if all group members share a common ancestor
            # that doesn't include members of other groups
            common_ancestor = tree.get_common_ancestor(group_nodes)
            descendant_leaves = [leaf.name for leaf in common_ancestor.get_leaves()]
            
            # Check if descendant leaves only contain group members
            other_group_samples = []
            for other_group, other_samples in introner_groups.items():
                if other_group != group_name:
                    other_group_samples.extend(other_samples)
            
            non_group_descendants = [leaf for leaf in descendant_leaves if leaf in other_group_samples]
            
            is_monophyletic = len(non_group_descendants) == 0
            
            monophyly_results[group_name] = {
                'is_monophyletic': is_monophyletic,
                'group_size': len(samples),
                'nodes_found': len(group_nodes),
                'descendant_leaves': descendant_leaves,
                'non_group_descendants': non_group_descendants
            }
            
            print(f"{group_name}: {'Monophyletic' if is_monophyletic else 'Not monophyletic'}")
            if not is_monophyletic:
                print(f"  Non-group descendants: {non_group_descendants}")
        else:
            print(f"{group_name}: Cannot test monophyly (insufficient samples in tree)")
            monophyly_results[group_name] = {'is_monophyletic': None, 'reason': 'insufficient_samples'}
    
    return monophyly_results

def compare_distance_matrices(phylo_distances, introner_distances):
    """
    Compare phylogenetic and introner-based distance matrices
    """
    print("Comparing phylogenetic and introner distance matrices...")
    
    if introner_distances is None:
        print("Warning: Introner distances not available for comparison")
        return None
    
    # Find common samples
    common_samples = list(set(phylo_distances.index) & set(introner_distances.index))
    common_samples.sort()
    
    print(f"Common samples for comparison: {len(common_samples)}")
    
    if len(common_samples) < 3:
        print("Error: Too few common samples for meaningful comparison")
        return None
    
    # Extract distance matrices for common samples
    phylo_common = phylo_distances.loc[common_samples, common_samples]
    introner_common = introner_distances.loc[common_samples, common_samples]
    
    # Get upper triangular elements (excluding diagonal)
    triu_indices = np.triu_indices_from(phylo_common.values, k=1)
    phylo_flat = phylo_common.values[triu_indices]
    introner_flat = introner_common.values[triu_indices]
    
    # Calculate correlations
    pearson_r, pearson_p = pearsonr(phylo_flat, introner_flat)
    spearman_r, spearman_p = spearmanr(phylo_flat, introner_flat)
    
    comparison_results = {
        'common_samples': common_samples,
        'n_comparisons': len(phylo_flat),
        'pearson_correlation': pearson_r,
        'pearson_p_value': pearson_p,
        'spearman_correlation': spearman_r,
        'spearman_p_value': spearman_p,
        'phylo_distances': phylo_common,
        'introner_distances': introner_common
    }
    
    print(f"Distance matrix comparison:")
    print(f"  Pearson correlation: r = {pearson_r:.4f}, p = {pearson_p:.4f}")
    print(f"  Spearman correlation: ρ = {spearman_r:.4f}, p = {spearman_p:.4f}")
    
    return comparison_results

def visualize_tree_with_groups(tree, introner_groups, output_dir):
    """
    Create tree visualization with introner groups color-coded
    """
    print("Creating tree visualization...")
    
    # Create color mapping for groups
    group_colors = {
        'Group1': '#FF6B6B',  # Red
        'Group2': '#4ECDC4',  # Teal
        'Outgroup': '#95E1D3', # Light green
        'Other': '#CCCCCC'     # Gray
    }
    
    # Assign colors to samples
    sample_colors = {}
    for group_name, samples in introner_groups.items():
        for sample in samples:
            sample_colors[sample] = group_colors[group_name]
    
    # Add outgroup colors
    outgroups = ['RCC1749', 'RCC3052']
    for outgroup in outgroups:
        sample_colors[outgroup] = group_colors['Outgroup']
    
    try:
        # Try to create a matplotlib-based tree plot
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Convert ete3 tree to simple format for plotting
        # This is a simplified approach - for publication quality, consider using ete3's TreeStyle
        
        # Get tree as newick string and use Bio.Phylo for plotting
        from io import StringIO
        from Bio import Phylo
        
        newick_str = tree.write(format=1)
        tree_bio = Phylo.read(StringIO(newick_str), 'newick')
        
        # Plot tree
        Phylo.draw(tree_bio, axes=ax, do_show=False)
        
        # Add title and legend
        ax.set_title('Phylogenetic Tree with Introner-Based Population Groups')
        
        # Create legend
        legend_elements = []
        for group, color in group_colors.items():
            if group in ['Group1', 'Group2', 'Outgroup']:
                legend_elements.append(Rectangle((0, 0), 1, 1, facecolor=color, label=group))
        
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/phylogenetic_tree_with_groups.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Tree visualization saved to {output_dir}/phylogenetic_tree_with_groups.png")
        
    except Exception as e:
        print(f"Warning: Could not create tree visualization: {e}")
        print("Tree topology (Newick format):")
        print(tree.write(format=1))

def create_distance_comparison_plots(comparison_results, output_dir):
    """
    Create plots comparing phylogenetic and introner distances
    """
    if comparison_results is None:
        print("No comparison results available for plotting")
        return
    
    print("Creating distance comparison plots...")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Scatter plot of distances
    phylo_flat = comparison_results['phylo_distances'].values[np.triu_indices_from(comparison_results['phylo_distances'].values, k=1)]
    introner_flat = comparison_results['introner_distances'].values[np.triu_indices_from(comparison_results['introner_distances'].values, k=1)]
    
    ax1.scatter(phylo_flat, introner_flat, alpha=0.6, s=30)
    ax1.set_xlabel('Phylogenetic Distance (Patristic)')
    ax1.set_ylabel('Introner Distance (1 - Consistency)')
    ax1.set_title(f'Distance Correlation\nr = {comparison_results["pearson_correlation"]:.3f}')
    
    # Add regression line
    z = np.polyfit(phylo_flat, introner_flat, 1)
    p = np.poly1d(z)
    ax1.plot(phylo_flat, p(phylo_flat), "r--", alpha=0.8)
    ax1.grid(True, alpha=0.3)
    
    # 2. Phylogenetic distance heatmap
    sns.heatmap(comparison_results['phylo_distances'], annot=True, fmt='.4f', 
                cmap='viridis', square=True, ax=ax2, 
                cbar_kws={'label': 'Patristic Distance'})
    ax2.set_title('Phylogenetic Distances')
    
    # 3. Introner distance heatmap
    sns.heatmap(comparison_results['introner_distances'], annot=True, fmt='.3f',
                cmap='plasma', square=True, ax=ax3,
                cbar_kws={'label': 'Introner Distance'})
    ax3.set_title('Introner Distances')
    
    # 4. Residuals plot
    predicted = p(phylo_flat)
    residuals = introner_flat - predicted
    
    ax4.scatter(phylo_flat, residuals, alpha=0.6, s=30)
    ax4.axhline(y=0, color='r', linestyle='--')
    ax4.set_xlabel('Phylogenetic Distance')
    ax4.set_ylabel('Residuals')
    ax4.set_title('Residuals from Linear Fit')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/distance_comparison_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Distance comparison plots saved to {output_dir}/distance_comparison_plots.png")

def main():
    """Main phylogenetic validation pipeline"""
    print("Starting Phylogenetic Validation of Introner Population Structure")
    print("=" * 65)
    
    # Create output directory
    output_dir = 'introner_analysis/phylogenetic_validation'
    os.makedirs(output_dir, exist_ok=True)
    
    # Phase 1: Load data and parse tree
    introner_groups, introner_distances = load_introner_population_groups()
    tree, all_samples = parse_and_root_tree()
    
    # Extract Group1 subtree for focused analysis
    group1_subtree = extract_group1_subtree(tree, introner_groups)
    
    # Phase 2: Calculate distances and compare
    all_group1_samples = introner_groups['Group1'] + introner_groups['Group2']
    phylo_distances = calculate_patristic_distances(tree, all_group1_samples)
    
    # Phase 3: Test monophyly
    monophyly_results = test_monophyly(tree, introner_groups)
    
    # Phase 2 continued: Compare distance matrices
    comparison_results = compare_distance_matrices(phylo_distances, introner_distances)
    
    # Phase 5: Create visualizations
    visualize_tree_with_groups(tree, introner_groups, output_dir)
    create_distance_comparison_plots(comparison_results, output_dir)
    
    # Save results
    results_summary = {
        'total_samples_in_tree': len(all_samples),
        'group1_samples_found': len([s for s in all_group1_samples if s in all_samples]),
        'monophyly_results': monophyly_results
    }
    
    if comparison_results:
        results_summary.update({
            'distance_correlation_pearson': comparison_results['pearson_correlation'],
            'distance_correlation_spearman': comparison_results['spearman_correlation'],
            'correlation_significance_p': comparison_results['pearson_p_value']
        })
    
    # Save phylogenetic distances
    phylo_distances.to_csv(f'{output_dir}/phylogenetic_distances.csv')
    
    # Save summary
    summary_df = pd.DataFrame([results_summary])
    summary_df.to_csv(f'{output_dir}/validation_summary.csv', index=False)
    
    print(f"\n✓ Phylogenetic validation completed successfully")
    print(f"  - Tree parsed and rooted with outgroups")
    print(f"  - Patristic distances calculated for {len(all_group1_samples)} samples")
    if comparison_results:
        print(f"  - Distance correlation: r = {comparison_results['pearson_correlation']:.4f}")
    print(f"  - Monophyly test results saved")
    print(f"  - Visualizations created")
    
    return {
        'tree': tree,
        'phylo_distances': phylo_distances,
        'introner_groups': introner_groups,
        'monophyly_results': monophyly_results,
        'comparison_results': comparison_results,
        'summary': results_summary
    }

if __name__ == "__main__":
    results = main()
    if results is not None:
        print(f"\n✓ Phylogenetic validation analysis completed successfully")
    else:
        print("\n✗ Phylogenetic validation analysis failed")
        sys.exit(1)
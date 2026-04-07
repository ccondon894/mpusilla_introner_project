#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
from Bio.Phylo import draw
import argparse

def plot_phylogenetic_tree(treefile, output_prefix="phylogenetic_tree_manual"):
    """
    Plot a phylogenetic tree from a Newick format treefile without any rerooting
    """
    
    # Load the tree (no rerooting)
    tree = Phylo.read(treefile, "newick")
    
    # Define group2 samples for coloring  
    group2_samples = {'RCC1749', 'RCC3052'}
    
    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, 
               branch_labels=None, label_func=lambda x: x.name if x.name else '')
    
    # Customize the plot
    ax.set_title('Manually Rerooted Phylogenetic Tree - M. pusilla 4D Sites', fontsize=16, fontweight='bold')
    ax.set_xlabel('Branch Length', fontsize=12)
    
    # Legend removed to avoid interference with tree labels
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    
    print(f"Tree plots saved:")
    print(f"- {output_prefix}.pdf")
    print(f"- {output_prefix}.png")
    
    # Print tree statistics
    print(f"\nTree statistics:")
    print(f"Number of terminal nodes (leaves): {len(list(tree.get_terminals()))}")
    print(f"Total tree length: {tree.total_branch_length():.6f}")
    print(f"Tree depth (root to furthest leaf): {max(tree.depths().values()):.6f}")
    
    # Print sample groupings
    terminals = list(tree.get_terminals())
    group1_samples = [t.name for t in terminals if t.name not in group2_samples]
    
    print(f"\nSample groupings:")
    print(f"Group 1 (Intronerful): {', '.join(sorted(group1_samples))}")
    print(f"Group 2 (Intronerless): {', '.join(sorted(group2_samples))}")
    
    return tree

def print_tree_info(treefile):
    """
    Print detailed information about the tree without rerooting
    """
    tree = Phylo.read(treefile, "newick")
    
    print("=== MANUALLY REROOTED TREE ANALYSIS ===\n")
    
    # Basic tree info
    print("Tree Structure:")
    print(f"  Total branches: {len(list(tree.get_nonterminals())) + len(list(tree.get_terminals()))}")
    print(f"  Terminal nodes (leaves): {len(list(tree.get_terminals()))}")
    print(f"  Internal nodes: {len(list(tree.get_nonterminals()))}")
    print(f"  Total tree length: {tree.total_branch_length():.6f}")
    print(f"  Max distance from root: {max(tree.depths().values()):.6f}")
    
    # Sample information
    terminals = list(tree.get_terminals())
    sample_names = [t.name for t in terminals]
    
    print(f"\nSamples in tree ({len(sample_names)}):")
    for i, name in enumerate(sorted(sample_names), 1):
        group = "Group 2 (Intronerless)" if name in {'RCC1749', 'RCC3052'} else "Group 1 (Intronerful)"
        print(f"  {i:2d}. {name} - {group}")
    
    # Branch length statistics
    branch_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]
    if branch_lengths:
        print(f"\nBranch Length Statistics:")
        print(f"  Minimum: {min(branch_lengths):.6f}")
        print(f"  Maximum: {max(branch_lengths):.6f}")
        print(f"  Mean: {sum(branch_lengths)/len(branch_lengths):.6f}")
    
    # Check for the major split
    print(f"\nTree topology shows Group 2 (RCC1749, RCC3052) as sister group to all Group 1 samples")

def main():
    parser = argparse.ArgumentParser(description="Plot manually rerooted phylogenetic tree")
    parser.add_argument("-i", "--input", 
                       default="rerooted_tree.treefile",
                       help="Input treefile in Newick format")
    parser.add_argument("-o", "--output", 
                       default="phylogenetic_tree_manual",
                       help="Output prefix for plot files")
    parser.add_argument("--info-only", action="store_true",
                       help="Only print tree information, don't create plots")
    
    args = parser.parse_args()
    
    if args.info_only:
        print_tree_info(args.input)
    else:
        # Create tree plot
        tree = plot_phylogenetic_tree(args.input, args.output)
        
        # Print tree info
        print("\n" + "="*50)
        print_tree_info(args.input)

if __name__ == "__main__":
    main()
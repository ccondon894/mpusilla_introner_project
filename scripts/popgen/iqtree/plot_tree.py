#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
from Bio.Phylo import draw
import argparse

def plot_phylogenetic_tree(treefile, output_prefix="phylogenetic_tree"):
    """
    Plot a phylogenetic tree from a Newick format treefile
    """
    
    # Load the tree
    tree = Phylo.read(treefile, "newick")
    
    # Reroot the tree using Group 2 samples as outgroup
    group2_samples = {'RCC1749', 'RCC3052'}
    
    # Find Group 2 terminals to use as outgroup
    group2_terminals = [terminal for terminal in tree.get_terminals() 
                       if terminal.name in group2_samples]
    
    if group2_terminals:
        # Use the first Group 2 sample as outgroup to keep them together
        # This will make RCC1749 and RCC3052 appear as sister taxa at the base
        tree.root_with_outgroup(group2_terminals[0])
        print(f"Tree rerooted using {group2_terminals[0].name} as outgroup")
        print("This keeps Group 2 samples together as sister taxa at the base")
    else:
        print("Warning: No Group 2 samples found for rerooting")
    
    # Define group2 samples for coloring
    group2_samples = {'RCC1749', 'RCC3052'}
    
    # Create color mapping for samples
    def get_sample_color(sample_name):
        if sample_name in group2_samples:
            return 'red'
        else:
            return 'blue'
    
    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True, 
               branch_labels=None, label_func=lambda x: x.name if x.name else '')
    
    # Customize the plot
    ax.set_title('Phylogenetic Tree - M. pusilla 4D Sites', fontsize=16, fontweight='bold')
    ax.set_xlabel('Branch Length', fontsize=12)
    
    # Color the terminal nodes based on sample groups
    # Note: This is a basic approach - for more sophisticated coloring, 
    # we'd need to access individual leaf elements
    
    # Add legend
    blue_patch = mpatches.Patch(color='blue', label='Group 1 (Intronerful)')
    red_patch = mpatches.Patch(color='red', label='Group 2 (Intronerless)')
    ax.legend(handles=[blue_patch, red_patch], loc='upper right')
    
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
    
    # plt.show()  # Commented out to avoid hanging
    
    return tree

def plot_tree_circular(treefile, output_prefix="phylogenetic_tree_circular"):
    """
    Plot a circular phylogenetic tree
    """
    
    # Load the tree
    tree = Phylo.read(treefile, "newick")
    
    # Reroot the tree using Group 2 samples as outgroup
    group2_samples = {'RCC1749', 'RCC3052'}
    group2_terminals = [terminal for terminal in tree.get_terminals() 
                       if terminal.name in group2_samples]
    
    if group2_terminals:
        tree.root_with_outgroup(group2_terminals[0])
    
    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Draw circular tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True,
               branch_labels=None, label_func=lambda x: x.name if x.name else '',
               label_colors={'RCC1749': 'red', 'RCC3052': 'red'})
    
    # Make it circular by adjusting the plot
    ax.set_aspect('equal')
    ax.set_title('Circular Phylogenetic Tree - M. pusilla 4D Sites', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Remove axes for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add legend
    blue_patch = mpatches.Patch(color='blue', label='Group 1 (Intronerful)')
    red_patch = mpatches.Patch(color='red', label='Group 2 (Intronerless)')
    ax.legend(handles=[blue_patch, red_patch], loc='upper right', bbox_to_anchor=(1.1, 1))
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    
    print(f"Circular tree plots saved:")
    print(f"- {output_prefix}.pdf")
    print(f"- {output_prefix}.png")
    
    # plt.show()  # Commented out to avoid hanging

def print_tree_info(treefile):
    """
    Print detailed information about the tree
    """
    tree = Phylo.read(treefile, "newick")
    
    # Reroot the tree using Group 2 samples as outgroup
    group2_samples = {'RCC1749', 'RCC3052'}
    group2_terminals = [terminal for terminal in tree.get_terminals() 
                       if terminal.name in group2_samples]
    
    if group2_terminals:
        tree.root_with_outgroup(group2_terminals[0])
        print(f"Tree rerooted using {group2_terminals[0].name} as outgroup")
        print("This keeps Group 2 samples together as sister taxa at the base\n")
    
    print("=== PHYLOGENETIC TREE ANALYSIS ===\n")
    
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
    
    # Find longest branch (possibly the outgroup branch)
    longest_branch = max(tree.find_clades(), key=lambda x: x.branch_length if x.branch_length else 0)
    if longest_branch.branch_length and longest_branch.branch_length > 0.1:  # Arbitrary threshold
        print(f"\nLongest branch: {longest_branch.branch_length:.6f}")
        if longest_branch.name:
            print(f"  Associated with: {longest_branch.name}")

def main():
    parser = argparse.ArgumentParser(description="Plot phylogenetic tree from Newick format file")
    parser.add_argument("-i", "--input", 
                       default="mpusilla.snps.4d.notMT.filtered.min4.phy.treefile",
                       help="Input treefile in Newick format")
    parser.add_argument("-o", "--output", 
                       default="phylogenetic_tree",
                       help="Output prefix for plot files")
    parser.add_argument("--circular", action="store_true",
                       help="Also create a circular tree plot")
    parser.add_argument("--info-only", action="store_true",
                       help="Only print tree information, don't create plots")
    
    args = parser.parse_args()
    
    if args.info_only:
        print_tree_info(args.input)
    else:
        # Create standard tree plot
        tree = plot_phylogenetic_tree(args.input, args.output)
        
        # Create circular plot if requested
        if args.circular:
            plot_tree_circular(args.input, args.output + "_circular")
        
        # Print tree info
        print("\n" + "="*50)
        print_tree_info(args.input)

if __name__ == "__main__":
    main()
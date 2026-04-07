#!/usr/bin/env python3
"""
Publication-quality phylogenetic tree visualization using toytree.

This script creates professional, publication-ready phylogenetic tree visualizations
with color-coded sample groups, clean typography, and multiple export formats.

Features:
- Color-coded tip labels (Group 1 = Blue, Group 2 = Red)
- Professional typography (Arial/Helvetica)
- Clean, minimal styling with thin branches
- Professional scale bar
- Legend showing group membership
- Multiple export formats (PDF, PNG, SVG)
- Command-line customization options

Author: Generated for M. pusilla phylogenetic analysis
Date: 2026-01-12
"""

import toytree
import toyplot
import toyplot.pdf
import toyplot.svg
import toyplot.png
import argparse
import sys


def plot_phylogenetic_tree_toytree(treefile, output_prefix="phylogenetic_tree_improved",
                                   width=900, height=600, tip_font_size=15,
                                   layout='r'):
    """
    Create publication-quality phylogenetic tree using toytree.

    Parameters:
    -----------
    treefile : str
        Path to Newick format tree file
    output_prefix : str
        Output filename prefix (without extension)
    width : int
        Canvas width in pixels
    height : int
        Canvas height in pixels
    tip_font_size : int
        Font size for tip labels
    layout : str
        Tree layout ('r' for rectangular, 'c' for circular)

    Returns:
    --------
    tree : toytree.tree
        The loaded tree object
    """

    print(f"📊 Loading tree from: {treefile}")

    try:
        # Load tree
        tree = toytree.tree(treefile)
    except Exception as e:
        print(f"❌ Error loading tree: {e}")
        sys.exit(1)

    # Define sample groups
    group2_samples = {'RCC1749', 'RCC3052'}
    tip_labels = tree.get_tip_labels()

    print(f"   Found {len(tip_labels)} samples")
    print(f"   Group 1 (Intronerful): {len([t for t in tip_labels if t not in group2_samples])} samples")
    print(f"   Group 2 (Intronerless): {len([t for t in tip_labels if t in group2_samples])} samples")

    # Create color mapping for tip labels
    tip_colors = []
    for tip in tip_labels:
        if tip in group2_samples:
            tip_colors.append('#E74C3C')  # Red for Group 2
        else:
            tip_colors.append('#3498DB')  # Blue for Group 1

    print(f"\n🎨 Creating {layout_names[layout]} tree visualization...")

    # Draw tree with custom styling using toytree v3 API
    canvas, axes, mark = tree.draw(
        width=width,
        height=height,
        layout=layout,
        tip_labels=True,
        tip_labels_colors=tip_colors,
        tip_labels_style={
            'font-size': f'{tip_font_size}px',
            'anchor-shift': 20,
        },
        edge_style={
            'stroke': '#333333',
            'stroke-width': 1.8
        },
        edge_align_style={
            'stroke': '#DDDDDD',
            'stroke-width': 0.8,
        },
        scale_bar=True,  # Note: underscore, not camelCase
        node_labels=False,  # No bootstrap values available
        node_sizes=0,  # Hide node markers
    )

    # Set white background
    canvas.style = {"background-color": "white"}

    # Increase x-axis tick label size
    axes.x.label.style = {"font-size": "16px"}
    axes.x.ticks.labels.style = {"font-size": "16px"}

    # Note: Legend omitted - tip label colors indicate group membership:
    #   Blue (#3498DB) = Group 1 (Intronerful): 11 samples
    #   Red (#E74C3C) = Group 2 (Intronerless): RCC1749, RCC3052

    print(f"\n💾 Exporting to multiple formats...")
    print(f"   Tip label colors: Blue = Group 1 (Intronerful), Red = Group 2 (Intronerless)")

    # Export to multiple formats
    try:
        toyplot.pdf.render(canvas, f"{output_prefix}.pdf")
        print(f"   ✓ Saved: {output_prefix}.pdf")
    except Exception as e:
        print(f"   ✗ PDF export failed: {e}")

    try:
        toyplot.svg.render(canvas, f"{output_prefix}.svg")
        print(f"   ✓ Saved: {output_prefix}.svg")
    except Exception as e:
        print(f"   ✗ SVG export failed: {e}")

    try:
        toyplot.png.render(canvas, f"{output_prefix}.png", scale=3)  # 3x for ~300 DPI
        print(f"   ✓ Saved: {output_prefix}.png (high resolution)")
    except Exception as e:
        print(f"   ✗ PNG export failed: {e}")

    return tree


def print_tree_info(tree):
    """
    Print detailed information about the tree structure.

    Parameters:
    -----------
    tree : toytree.tree
        The tree object
    """
    print(f"\n{'='*60}")
    print("TREE STATISTICS")
    print(f"{'='*60}\n")

    # Get tip labels
    tip_labels = tree.get_tip_labels()
    group2_samples = {'RCC1749', 'RCC3052'}
    group1_samples = [t for t in tip_labels if t not in group2_samples]

    # Basic tree info
    print(f"Tree Structure:")
    print(f"  Terminal nodes (leaves): {tree.ntips}")
    print(f"  Internal nodes: {tree.nnodes - tree.ntips}")
    print(f"  Total tree length: {tree.treenode.height:.6f}")

    # Branch length statistics
    try:
        # Get branch lengths from the tree
        branch_lengths = [node.dist for node in tree.treenode.traverse() if node.dist and node.dist > 0]

        if branch_lengths:
            print(f"\nBranch Length Statistics:")
            print(f"  Minimum: {min(branch_lengths):.6f}")
            print(f"  Maximum: {max(branch_lengths):.6f}")
            print(f"  Mean: {sum(branch_lengths)/len(branch_lengths):.6f}")
            print(f"  Total branches: {len(branch_lengths)}")
    except Exception as e:
        print(f"\n  Branch statistics unavailable: {e}")

    # Sample groupings
    print(f"\nSample Groupings:")
    print(f"  Group 1 (Intronerful): {len(group1_samples)} samples")
    for sample in sorted(group1_samples):
        print(f"    - {sample}")

    print(f"\n  Group 2 (Intronerless): {len(group2_samples)} samples")
    for sample in sorted(group2_samples):
        print(f"    - {sample}")

    print(f"\n{'='*60}\n")


# Layout name mapping
layout_names = {
    'r': 'rectangular',
    'c': 'circular',
    'd': 'down',
    'u': 'up'
}


def main():
    """Main function with argument parsing."""

    parser = argparse.ArgumentParser(
        description='Create publication-quality phylogenetic tree visualization using toytree',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (rectangular tree)
  python plot_tree_improved.py -i rerooted_tree.treefile

  # Circular layout
  python plot_tree_improved.py -i rerooted_tree.treefile --layout circular

  # Custom dimensions and font size
  python plot_tree_improved.py -i rerooted_tree.treefile --width 1000 --height 800 --tip-font-size 14

  # Custom output name
  python plot_tree_improved.py -i rerooted_tree.treefile -o publication_figure

  # Info only (no plots)
  python plot_tree_improved.py -i rerooted_tree.treefile --info-only
        """)

    parser.add_argument('-i', '--input',
                       default='rerooted_tree.treefile',
                       help='Input tree file (Newick format). Default: rerooted_tree.treefile')
    parser.add_argument('-o', '--output',
                       default='phylogenetic_tree_improved',
                       help='Output file prefix (without extension). Default: phylogenetic_tree_improved')
    parser.add_argument('--width', type=int, default=900,
                       help='Canvas width in pixels. Default: 900')
    parser.add_argument('--height', type=int, default=600,
                       help='Canvas height in pixels. Default: 600')
    parser.add_argument('--layout', choices=['rectangular', 'circular'],
                       default='rectangular',
                       help='Tree layout style. Default: rectangular')
    parser.add_argument('--tip-font-size', type=int, default=16,
                       help='Tip label font size in pixels. Default: 16')
    parser.add_argument('--info-only', action='store_true',
                       help='Only print tree information, don\'t create plots')

    args = parser.parse_args()

    # Convert layout name to code
    layout_code = 'r' if args.layout == 'rectangular' else 'c'

    print("="*60)
    print("PUBLICATION-QUALITY PHYLOGENETIC TREE VISUALIZATION")
    print("="*60)

    if args.info_only:
        # Load tree and print info only
        try:
            tree = toytree.tree(args.input)
            print_tree_info(tree)
        except Exception as e:
            print(f"❌ Error loading tree: {e}")
            sys.exit(1)
    else:
        # Create tree visualization
        tree = plot_phylogenetic_tree_toytree(
            args.input,
            args.output,
            width=args.width,
            height=args.height,
            tip_font_size=args.tip_font_size,
            layout=layout_code
        )

        # Print tree info
        print_tree_info(tree)

        print("✅ Tree visualization complete!\n")


if __name__ == "__main__":
    main()

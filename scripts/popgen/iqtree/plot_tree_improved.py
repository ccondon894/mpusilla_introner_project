#!/usr/bin/env python3
"""
Publication-quality phylogenetic tree visualization using toytree.

Colors branches and tip labels by population (from master_figure_color_guide.tsv),
and adds group labels, legend, and a labeled scale bar.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import toytree
import toyplot.pdf
import toyplot.png
import toyplot.svg

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide

GROUP2_SAMPLES = {"RCC1749", "RCC3052"}
REFERENCE_SAMPLE = "CCMP1545"

LAYOUT_NAMES = {
    "r": "rectangular",
    "c": "circular",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create publication-quality phylogenetic tree visualization using toytree",
    )
    parser.add_argument(
        "-i",
        "--input",
        default="rerooted_tree.treefile",
        help="Input tree file (Newick format)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="phylogenetic_tree_improved",
        help="Output file prefix (without extension)",
    )
    parser.add_argument(
        "--color-guide",
        default=str(DEFAULT_GUIDE_PATH),
        help="Master figure color guide TSV",
    )
    parser.add_argument("--width", type=int, default=720, help="Canvas width in pixels")
    parser.add_argument("--height", type=int, default=380, help="Canvas height in pixels")
    parser.add_argument(
        "--layout",
        choices=["rectangular", "circular"],
        default="rectangular",
        help="Tree layout style",
    )
    parser.add_argument("--tip-font-size", type=int, default=14, help="Tip label font size")
    parser.add_argument(
        "--reference-sample",
        default=REFERENCE_SAMPLE,
        help="Reference sample highlighted in tip labels",
    )
    parser.add_argument(
        "--info-only",
        action="store_true",
        help="Only print tree information, do not create plots",
    )
    return parser.parse_args()


def load_tree(treefile):
    try:
        return toytree.tree(treefile)
    except Exception as exc:
        raise SystemExit(f"Error loading tree {treefile}: {exc}") from exc


def node_descendant_tips(tree):
    tips_by_node = {}
    for node in tree.treenode.traverse("postorder"):
        if node.is_leaf():
            tips_by_node[node.idx] = {node.name}
        else:
            tips_by_node[node.idx] = set().union(
                *(tips_by_node[child.idx] for child in node.children)
            )
    return tips_by_node


def population_edge_colors(tree, tips_by_node, group2_samples, pop1_color, pop2_color):
    colors = []
    for idx in range(tree.nnodes):
        tips = tips_by_node.get(idx, set())
        if not tips:
            colors.append("#CCCCCC")
        elif tips.issubset(group2_samples):
            colors.append(pop2_color)
        else:
            colors.append(pop1_color)
    return colors


def tip_display_and_colors(tip_labels, group2_samples, reference_sample, guide):
    pop1 = get_color("Population 1", guide)
    pop1_dark = get_color("Population 1 dark", guide)
    pop2 = get_color("Population 2", guide)

    display_labels = []
    colors = []
    for tip in tip_labels:
        if tip in group2_samples:
            display_labels.append(tip)
            colors.append(pop2)
        else:
            display_labels.append(tip)
            colors.append(pop1)
    return display_labels, colors, pop1, pop2


def add_group_labels(axes, mark, tree, group2_samples, pop1_color, pop2_color):
    ntable = mark.ntable
    xmin = float(min(row[0] for row in ntable))
    label_x = xmin - 0.045

    tip_labels = tree.get_tip_labels()
    group2_indices = [idx for idx, tip in enumerate(tip_labels) if tip in group2_samples]
    group1_indices = [idx for idx, tip in enumerate(tip_labels) if tip not in group2_samples]

    def y_range(indices):
        ys = [float(ntable[idx][1]) for idx in indices]
        return min(ys), max(ys)

    g2_y0, g2_y1 = y_range(group2_indices)
    g1_y0, g1_y1 = y_range(group1_indices)


def add_legend(axes, mark, pop1_color, pop2_color):
    ntable = mark.ntable
    xmin = float(min(row[0] for row in ntable))
    ymax = float(max(row[1] for row in ntable))
    lx = xmin - 0.01
    ly = ymax + 0.3

    entries = [
        (pop1_color, "Population 1"),
        (pop2_color, "Population 2"),
    ]
    for offset, (color, label) in enumerate(entries):
        y = ly - (offset * 0.55)
        axes.text(lx, y, "●", style={"fill": color, "font-size": "13px"})
        axes.text(
            lx + 0.012,
            y,
            label,
            style={"fill": "#333333", "font-size": "12px", "text-anchor": "start"},
        )


def plot_phylogenetic_tree(
    treefile,
    output_prefix,
    color_guide_path,
    width=520,
    height=380,
    layout="r",
    tip_font_size=14,
    reference_sample=REFERENCE_SAMPLE,
):
    guide = load_color_guide(color_guide_path)
    tree = load_tree(treefile)
    tip_labels = tree.get_tip_labels()
    tips_by_node = node_descendant_tips(tree)

    display_tips, tip_colors, pop1_color, pop2_color = tip_display_and_colors(
        tip_labels,
        GROUP2_SAMPLES,
        reference_sample,
        guide,
    )
    edge_colors = population_edge_colors(
        tree,
        tips_by_node,
        GROUP2_SAMPLES,
        pop1_color,
        pop2_color,
    )
    print(f"Loaded tree from: {treefile}")
    print(f"  Tips: {tree.ntips}")
    print(f"  Group 1: {len([t for t in tip_labels if t not in GROUP2_SAMPLES])}")
    print(f"  Group 2: {len([t for t in tip_labels if t in GROUP2_SAMPLES])}")
    print(f"  Layout: {LAYOUT_NAMES[layout]}")

    canvas, axes, mark = tree.draw(
        width=width,
        height=height,
        layout=layout,
        tip_labels=display_tips,
        tip_labels_colors=tip_colors,
        tip_labels_style={
            "font-size": f"{tip_font_size}px",
            "anchor-shift": 22,
        },
        edge_colors=edge_colors,
        edge_style={"stroke-width": 2.2},
        edge_align_style={"stroke": "#E6E6E6", "stroke-width": 0.6},
        node_labels=False,
        node_sizes=0,
        scale_bar=True,
        padding=36,
    )

    canvas.style = {"background-color": "white"}
    scale_label = "Substitutions per 4-fold degenerate site"
    axes.x.label.text = scale_label
    axes.x.label.style = {"font-size": "14px"}
    axes.x.ticks.labels.style = {"font-size": "12px"}

    if layout == "r":
        add_group_labels(axes, mark, tree, GROUP2_SAMPLES, pop1_color, pop2_color)
        add_legend(axes, mark, pop1_color, pop2_color)
        ymin = float(min(row[1] for row in mark.ntable))
        ymax = float(max(row[1] for row in mark.ntable))
        axes.y.domain.min = ymin - 0.8
        axes.y.domain.max = ymax + 1.7

    export_formats(canvas, output_prefix)
    return tree


def export_formats(canvas, output_prefix):
    for renderer, extension in (
        (toyplot.pdf.render, "pdf"),
        (toyplot.svg.render, "svg"),
        (lambda canvas, path: toyplot.png.render(canvas, path, scale=3), "png"),
    ):
        path = f"{output_prefix}.{extension}"
        try:
            renderer(canvas, path)
            print(f"Saved: {path}")
        except Exception as exc:
            print(f"Failed to save {path}: {exc}")


def print_tree_info(tree):
    tip_labels = tree.get_tip_labels()
    group1_samples = [tip for tip in tip_labels if tip not in GROUP2_SAMPLES]

    print("\n" + "=" * 60)
    print("TREE STATISTICS")
    print("=" * 60)
    print(f"Terminal nodes: {tree.ntips}")
    print(f"Internal nodes: {tree.nnodes - tree.ntips}")
    print(f"Total tree length: {tree.treenode.height:.6f}")
    print("\nGroup 1 samples:")
    for sample in sorted(group1_samples):
        print(f"  - {sample}")
    print("\nGroup 2 samples:")
    for sample in sorted(GROUP2_SAMPLES):
        print(f"  - {sample}")
    print("=" * 60 + "\n")


def main():
    args = parse_args()
    tree = load_tree(args.input)

    if args.info_only:
        print_tree_info(tree)
        return

    tree = plot_phylogenetic_tree(
        args.input,
        args.output,
        args.color_guide,
        width=args.width,
        height=args.height,
        layout="r" if args.layout == "rectangular" else "c",
        tip_font_size=args.tip_font_size,
        reference_sample=args.reference_sample,
    )
    print_tree_info(tree)


if __name__ == "__main__":
    main()

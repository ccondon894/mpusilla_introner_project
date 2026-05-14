#!/usr/bin/env python3
"""
Shannon Diversity and Evenness Analysis for Isoforms (Consolidated)
=====================================================================

Dual-mode script:
  default:          Calculate Shannon H and J' per gene from SQANTI3 data
  compare_introner: Compare diversity between introner/non-introner genes

H = -sum(p_i * ln(p_i))  where p_i is proportion of isoform i
J' = H / ln(k)           where k is number of isoforms

Usage (Snakemake):
  Default mode:
    python shannon_diversity_analysis_consolidated.py \
        --sqanti TSV --gtf GTF --output CSV --summary TXT --plot PDF \
        --mt_scaffold STR --mt_start INT --mt_end INT

  Compare mode:
    python shannon_diversity_analysis_consolidated.py \
        --mode compare_introner \
        --diversity CSV --genotype_matrix TSV \
        --output CSV --plot PDF --stats TXT
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import re
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

# Structural categories to include
INCLUDE_CATEGORIES = {'full-splice_match', 'novel_in_catalog', 'novel_not_in_catalog'}


def save_figures_with_png(figures, output_path):
    """Write a multi-page PDF and a single PNG contact sheet of the same figures."""
    output_path = Path(output_path)
    with PdfPages(output_path) as pdf:
        for fig in figures:
            pdf.savefig(fig)

    if output_path.suffix.lower() == ".pdf" and figures:
        rendered = []
        for fig in figures:
            fig.canvas.draw()
            rendered.append(np.asarray(fig.canvas.buffer_rgba()))

        fig_height = max(4, 3.5 * len(rendered))
        summary_fig, axes = plt.subplots(len(rendered), 1, figsize=(12, fig_height), squeeze=False)
        for ax, image in zip(axes.flatten(), rendered):
            ax.imshow(image)
            ax.axis("off")
        summary_fig.tight_layout()
        summary_fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
        plt.close(summary_fig)


def calculate_shannon_index(proportions):
    """
    Calculate Shannon diversity index.
    H = -sum(p_i * ln(p_i))
    """
    p = proportions[proportions > 0]
    H = -np.sum(p * np.log(p))
    return H


def calculate_shannon_evenness(shannon_index, num_species):
    """
    Calculate Shannon evenness (normalized diversity).
    J' = H / ln(k)
    """
    if num_species <= 1:
        return np.nan
    return shannon_index / np.log(num_species)


def extract_mating_type_genes(gtf_path, scaffold, start, end):
    """
    Extract gene IDs in the mating-type region from a GTF file.
    Handles scaffold names with or without pangenome prefixes.
    """
    mating_type_genes = set()

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            attributes = fields[8]

            # Match scaffold name via substring (handles CCMP1545#0#scaffold_2)
            if scaffold not in chrom:
                continue

            if feature_start <= end and feature_end >= start:
                match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', attributes)
                if match:
                    mating_type_genes.add(match.group(1))

    return mating_type_genes


def identify_count_cols(df):
    """Identify replicate count columns in SQANTI3 data."""
    count_cols = [col for col in df.columns
                  if re.match(r'^\d+_[A-Z]\d+$', col)]
    return count_cols


def run_default_mode(args):
    """Calculate Shannon diversity per gene from SQANTI3 data."""
    # Load merged SQANTI3 data
    data = pd.read_csv(args.sqanti, sep='\t')
    n_raw = len(data)

    # Filter by structural category
    data = data[data['structural_category'].isin(INCLUDE_CATEGORIES)]
    n_after_cat = len(data)

    # Identify count columns
    count_cols = identify_count_cols(data)
    if not count_cols:
        # Fallback: try numeric columns that aren't metadata
        meta_cols = {'length', 'exons', 'ORF_length', 'CDS_length'}
        count_cols = [c for c in data.select_dtypes(include=[np.number]).columns
                      if c not in meta_cols]

    # Calculate total count per isoform
    data['total_count'] = data[count_cols].sum(axis=1)

    # Ensure gene_id column exists
    if 'gene_id' not in data.columns:
        data['gene_id'] = data['associated_gene']

    # Filter mating-type region genes
    mt_genes = extract_mating_type_genes(
        args.gtf, args.mt_scaffold, args.mt_start, args.mt_end
    )
    data = data[~data['gene_id'].isin(mt_genes)]
    n_after_mt = len(data)

    # Calculate Shannon metrics per gene
    shannon_results = []

    for gene_id in data['gene_id'].unique():
        gene_data = data[data['gene_id'] == gene_id]
        num_isoforms = len(gene_data)

        if num_isoforms < 2:
            continue

        gene_total = gene_data['total_count'].sum()
        if gene_total == 0:
            continue

        proportions = gene_data['total_count'].values / gene_total
        h = calculate_shannon_index(proportions)
        j = calculate_shannon_evenness(h, num_isoforms)

        dominant_category = gene_data.loc[
            gene_data['total_count'].idxmax(), 'structural_category'
        ]

        shannon_results.append({
            'gene_id': gene_id,
            'num_isoforms': num_isoforms,
            'shannon_index': h,
            'shannon_evenness': j,
            'dominant_category': dominant_category,
            'num_fsm': (gene_data['structural_category'] == 'full-splice_match').sum(),
            'num_novel_in_catalog': (gene_data['structural_category'] == 'novel_in_catalog').sum(),
            'num_novel_not_in_catalog': (gene_data['structural_category'] == 'novel_not_in_catalog').sum()
        })

    shannon_df = pd.DataFrame(shannon_results)
    shannon_df.to_csv(args.output, index=False)

    # Write summary
    with open(args.summary, 'w') as f:
        f.write("Shannon Diversity Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Raw isoforms loaded: {n_raw}\n")
        f.write(f"After structural category filter: {n_after_cat}\n")
        f.write(f"After mating-type filter: {n_after_mt}\n")
        f.write(f"Mating-type genes excluded: {len(mt_genes)}\n")
        f.write(f"Genes analyzed (>=2 isoforms): {len(shannon_df)}\n\n")

        f.write("Shannon Index (H):\n")
        f.write(f"  Mean: {shannon_df['shannon_index'].mean():.4f}\n")
        f.write(f"  Median: {shannon_df['shannon_index'].median():.4f}\n")
        f.write(f"  Std Dev: {shannon_df['shannon_index'].std():.4f}\n\n")

        f.write("Shannon Evenness (J'):\n")
        f.write(f"  Mean: {shannon_df['shannon_evenness'].mean():.4f}\n")
        f.write(f"  Median: {shannon_df['shannon_evenness'].median():.4f}\n")
        f.write(f"  Std Dev: {shannon_df['shannon_evenness'].std():.4f}\n\n")

        high = (shannon_df['shannon_evenness'] > 0.7).sum()
        moderate = ((shannon_df['shannon_evenness'] >= 0.4) &
                    (shannon_df['shannon_evenness'] <= 0.7)).sum()
        low = (shannon_df['shannon_evenness'] < 0.4).sum()
        f.write("Diversity levels:\n")
        f.write(f"  High (J' > 0.7): {high} ({high/len(shannon_df)*100:.1f}%)\n")
        f.write(f"  Moderate (0.4-0.7): {moderate} ({moderate/len(shannon_df)*100:.1f}%)\n")
        f.write(f"  Low (J' < 0.4): {low} ({low/len(shannon_df)*100:.1f}%)\n")

    sns.set_style("whitegrid")
    figures = []

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(shannon_df['shannon_index'], bins=50, color='#3498db',
            alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(shannon_df['shannon_index'].mean(), color='red',
               linestyle='--', linewidth=2,
               label=f"Mean: {shannon_df['shannon_index'].mean():.3f}")
    ax.axvline(shannon_df['shannon_index'].median(), color='green',
               linestyle='--', linewidth=2,
               label=f"Median: {shannon_df['shannon_index'].median():.3f}")
    ax.set_xlabel('Shannon Index (H)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Shannon Diversity Index', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    fig.tight_layout()
    figures.append(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(shannon_df['shannon_evenness'].dropna(), bins=50,
            color='#2ecc71', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(shannon_df['shannon_evenness'].mean(), color='red',
               linestyle='--', linewidth=2,
               label=f"Mean: {shannon_df['shannon_evenness'].mean():.3f}")
    ax.set_xlabel("Shannon Evenness (J')", fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Shannon Evenness', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    ax.legend(fontsize=10)
    fig.tight_layout()
    figures.append(fig)

    fsm_dom = shannon_df[shannon_df['dominant_category'] == 'full-splice_match']
    novel_dom = shannon_df[shannon_df['dominant_category'] != 'full-splice_match']

    if len(fsm_dom) > 0 and len(novel_dom) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, metric, label in [(axes[0], 'shannon_index', 'Shannon Index (H)'),
                                   (axes[1], 'shannon_evenness', "Shannon Evenness (J')")]:
            bp = ax.boxplot([fsm_dom[metric].dropna(), novel_dom[metric].dropna()],
                            labels=['FSM Dominant', 'Novel Dominant'],
                            patch_artist=True)
            for patch, color in zip(bp['boxes'], ['#2ecc71', '#e74c3c']):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            ax.set_ylabel(label, fontsize=12, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        fig.tight_layout()
        figures.append(fig)

    save_figures_with_png(figures, args.plot)
    for fig in figures:
        plt.close(fig)

    print(f"Shannon diversity calculated for {len(shannon_df)} genes")
    print(f"Output: {args.output}")


def run_compare_introner_mode(args):
    """Compare Shannon diversity between introner and non-introner genes."""
    # Load diversity data
    diversity = pd.read_csv(args.diversity)

    # Load genotype matrix
    genotype = pd.read_csv(args.genotype_matrix, sep='\t')

    # Identify introner genes (any gene with presence=1 in any sample)
    introner_genes = set(
        genotype[genotype['presence'] == 1]['gene'].unique()
    )

    # Classify genes
    diversity['has_introner'] = diversity['gene_id'].isin(introner_genes)

    with_introner = diversity[diversity['has_introner']]
    without_introner = diversity[~diversity['has_introner']]

    # Mann-Whitney U tests
    results = {}
    for metric in ['shannon_index', 'shannon_evenness']:
        a = with_introner[metric].dropna()
        b = without_introner[metric].dropna()
        if len(a) > 0 and len(b) > 0:
            u_stat, p_val = stats.mannwhitneyu(a, b, alternative='two-sided')
        else:
            u_stat, p_val = np.nan, np.nan
        results[metric] = {'u_stat': u_stat, 'p_val': p_val}

    # Save comparison CSV
    diversity.to_csv(args.output, index=False)

    # Save stats
    with open(args.stats, 'w') as f:
        f.write("Diversity by Introner Status\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Genes with introners: {len(with_introner)}\n")
        f.write(f"Genes without introners: {len(without_introner)}\n\n")

        for metric in ['shannon_index', 'shannon_evenness']:
            label = "Shannon Index (H)" if metric == 'shannon_index' else "Shannon Evenness (J')"
            f.write(f"{label}:\n")
            f.write(f"  With introners:    mean={with_introner[metric].mean():.4f}, "
                    f"median={with_introner[metric].median():.4f}\n")
            f.write(f"  Without introners: mean={without_introner[metric].mean():.4f}, "
                    f"median={without_introner[metric].median():.4f}\n")
            r = results[metric]
            sig = "SIGNIFICANT" if r['p_val'] < 0.05 else "NOT SIGNIFICANT"
            f.write(f"  Mann-Whitney U: U={r['u_stat']:.1f}, p={r['p_val']:.4e} ({sig})\n\n")

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, metric, label in [
        (axes[0], 'shannon_index', 'Shannon Index (H)'),
        (axes[1], 'shannon_evenness', "Shannon Evenness (J')")
    ]:
        a = with_introner[metric].dropna()
        b = without_introner[metric].dropna()
        bp = ax.boxplot([a, b],
                        labels=['With Introners', 'Without Introners'],
                        patch_artist=True)
        for patch, color in zip(bp['boxes'], ['#e74c3c', '#3498db']):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        r = results[metric]
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.set_title(f'{label}\n(p = {r["p_val"]:.2e})',
                     fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    save_figures_with_png([fig], args.plot)
    plt.close(fig)

    print(f"Comparison complete: {len(with_introner)} introner vs "
          f"{len(without_introner)} non-introner genes")


def main():
    parser = argparse.ArgumentParser(
        description="Shannon diversity analysis for isoforms"
    )
    parser.add_argument('--mode', default='default',
                        choices=['default', 'compare_introner'],
                        help='Analysis mode')

    # Default mode args
    parser.add_argument('--sqanti', help='Merged SQANTI3 TSV')
    parser.add_argument('--gtf', help='GTF file for MT region extraction')
    parser.add_argument('--mt_scaffold', help='Mating-type scaffold name')
    parser.add_argument('--mt_start', type=int, help='Mating-type region start')
    parser.add_argument('--mt_end', type=int, help='Mating-type region end')
    parser.add_argument('--summary', help='Summary text output')

    # Compare mode args
    parser.add_argument('--diversity', help='Shannon diversity CSV (compare mode)')
    parser.add_argument('--genotype_matrix', help='Genotype matrix TSV')
    parser.add_argument('--stats', help='Statistics output (compare mode)')

    # Common args
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--plot', required=True, help='Output PDF plot (also writes sibling PNG)')

    args = parser.parse_args()

    if args.mode == 'default':
        if not all([args.sqanti, args.gtf, args.summary,
                    args.mt_scaffold, args.mt_start is not None,
                    args.mt_end is not None]):
            parser.error("Default mode requires --sqanti, --gtf, --summary, "
                         "--mt_scaffold, --mt_start, --mt_end")
        run_default_mode(args)
    elif args.mode == 'compare_introner':
        if not all([args.diversity, args.genotype_matrix, args.stats]):
            parser.error("Compare mode requires --diversity, "
                         "--genotype_matrix, --stats")
        run_compare_introner_mode(args)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Analyze NMD predictions and introner presence for multi-exon transcript isoforms.

Dual-mode script:
  default:   Analyze NMD × introner association using merged SQANTI3 data
  by_strain: Per-strain NMD analysis with separate SQANTI3 files, GTFs,
             and introner BED files

Usage (Snakemake):
  Default:
    python analyze_nmd_predictions.py \
        --sqanti TSV --gtf GTF --introner_bed BED \
        --output TXT --contingency CSV --plot PDF \
        --categories STR

  By-strain:
    python analyze_nmd_predictions.py \
        --mode by_strain \
        --sqanti_ccmp1545 TXT --sqanti_rcc1614 TXT --sqanti_rcc1749 TXT \
        --gtf_ccmp1545 GTF --gtf_rcc1614 GTF --gtf_rcc1749 GTF \
        --bed_ccmp1545 BED --bed_rcc1614 BED --bed_rcc1749 BED \
        --mt_gtf CCMP1545.gtf \
        --output TXT --plot PDF
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
from collections import defaultdict
from scipy.stats import fisher_exact
from pathlib import Path


# ============================================================
# CORE ANALYSIS FUNCTIONS (preserved from original)
# ============================================================

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

        summary_fig, axes = plt.subplots(
            len(rendered),
            1,
            figsize=(10, max(4, 3.5 * len(rendered))),
            squeeze=False,
        )
        for ax, image in zip(axes.flatten(), rendered):
            ax.imshow(image)
            ax.axis("off")
        summary_fig.tight_layout()
        summary_fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
        plt.close(summary_fig)

def parse_gtf_file(gtf_path):
    """Parse GTF file and extract gene exon coordinates."""
    gene_exons = defaultdict(list)

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue

            chrom = fields[0]
            start_0based = int(fields[3]) - 1
            end_0based = int(fields[4])

            gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                gene_exons[gene_id].append((chrom, start_0based, end_0based))

    return dict(gene_exons)


def parse_bed_file(bed_path):
    """Parse BED file and extract introner coordinates by chromosome."""
    introners_by_chrom = defaultdict(list)

    with open(bed_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            introners_by_chrom[chrom].append((start, end))

    return dict(introners_by_chrom)


def clean_gene_name(associated_gene):
    """Clean SQANTI3 associated_gene to match GTF gene_id format."""
    if pd.isna(associated_gene) or associated_gene == '':
        return []

    candidates = [associated_gene]

    if associated_gene.startswith('novelGene_'):
        cleaned = associated_gene[10:]
        candidates.append(cleaned)
        if cleaned.endswith('_AS'):
            candidates.append(cleaned[:-3])

    if associated_gene.endswith('_AS'):
        candidates.append(associated_gene[:-3])

    gene_patterns = re.findall(r'([A-Za-z0-9_]+\.\d+\.\d+\.\d+\.\d+)', associated_gene)
    candidates.extend(gene_patterns)

    seen = set()
    return [c for c in candidates if not (c in seen or seen.add(c))]


def has_introner_overlap(gene_exons, introners_by_chrom):
    """Check if any introners overlap with gene exons."""
    for chrom, exon_start, exon_end in gene_exons:
        if chrom not in introners_by_chrom:
            continue
        for introner_start, introner_end in introners_by_chrom[chrom]:
            if introner_start < exon_end and introner_end > exon_start:
                return True
    return False


def calculate_fisher_test(a, b, c, d):
    """Perform Fisher's exact test on 2x2 table."""
    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table, alternative='two-sided')

    if a > 0 and b > 0 and c > 0 and d > 0:
        log_or = np.log(odds_ratio)
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        ci_lower = np.exp(log_or - 1.96 * se)
        ci_upper = np.exp(log_or + 1.96 * se)
    else:
        ci_lower = ci_upper = None

    return {
        'odds_ratio': odds_ratio, 'p_value': p_value,
        'ci_lower': ci_lower, 'ci_upper': ci_upper
    }


def extract_mating_type_genes(gtf_path, scaffold='scaffold_2',
                               start=49808, end=1730591):
    """Extract gene IDs in mating-type region."""
    mt_genes = set()
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            if scaffold not in chrom:
                continue
            feat_start = int(fields[3])
            feat_end = int(fields[4])
            if feat_start <= end and feat_end >= start:
                match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', fields[8])
                if match:
                    mt_genes.add(match.group(1))
    return mt_genes


def format_output(results):
    """Format analysis results as text report."""
    lines = []
    lines.append("NMD and Introner Analysis for Multi-Exon Transcript Isoforms")
    lines.append("=" * 70)
    lines.append("")

    for result in results:
        lines.append(f"Sample: {result.get('label', 'All')} "
                     f"(Genome: {result.get('genome', 'merged')})")
        lines.append("-" * 70)

        lines.append("Data Filtering:")
        lines.append(f"  Initial isoforms: {result['n_initial']:,}")
        lines.append(f"  After structural category filter: "
                     f"{result['n_after_category']:,} "
                     f"(removed {result['n_initial'] - result['n_after_category']:,})")
        lines.append(f"  After mating-type region filter: "
                     f"{result['n_after_mt']:,} "
                     f"(removed {result['n_after_category'] - result['n_after_mt']:,})")
        lines.append(f"  Multi-exon transcripts: {result['total_multi_exon']:,} "
                     f"(removed {result['n_after_mt'] - result['total_multi_exon']:,})")
        if 'categories' in result:
            lines.append(f"  Categories included: {', '.join(result['categories'])}")
        lines.append("")

        lines.append("NMD Prediction Analysis:")
        total = result['total_multi_exon']
        nmd_pos = result['count_nmd_positive']
        nmd_neg = result['count_nmd_negative']
        with_pred = result['count_with_prediction']
        pct_pred = (with_pred / total * 100) if total > 0 else 0
        pct_pos = (nmd_pos / total * 100) if total > 0 else 0
        pct_pos_pred = (nmd_pos / with_pred * 100) if with_pred > 0 else 0
        pct_neg = (nmd_neg / total * 100) if total > 0 else 0
        pct_neg_pred = (nmd_neg / with_pred * 100) if with_pred > 0 else 0
        lines.append(f"  Multi-exon with NMD predictions: {with_pred:,} ({pct_pred:.1f}%)")
        lines.append(f"  Multi-exon predicted as NMD+: {nmd_pos:,} "
                     f"({pct_pos:.1f}% of all multi-exon, "
                     f"{pct_pos_pred:.1f}% of those with predictions)")
        lines.append(f"  Multi-exon predicted as NMD-: {nmd_neg:,} "
                     f"({pct_neg:.1f}% of all multi-exon, "
                     f"{pct_neg_pred:.1f}% of those with predictions)")
        lines.append("")

        lines.append("Introner Analysis:")
        pct_with = (result['count_with_introners'] / total * 100) if total > 0 else 0
        pct_without = (result['count_without_introners'] / total * 100) if total > 0 else 0
        lines.append(f"  Multi-exon with introners: {result['count_with_introners']:,} ({pct_with:.1f}%)")
        lines.append(f"  Multi-exon without introners: {result['count_without_introners']:,} ({pct_without:.1f}%)")
        lines.append(f"  Could not map to GTF: {result['not_found_count']:,}")
        lines.append("")

        lines.append("Cross-tabulation (NMD x Introner):")
        lines.append(f"{'':25} {'With Introners':>20} {'Without Introners':>20} {'Total':>10}")
        nmd_true_total = result['nmd_true_with'] + result['nmd_true_without']
        nmd_false_total = result['nmd_false_with'] + result['nmd_false_without']
        no_pred_total = result['no_pred_with'] + result['no_pred_without']
        lines.append(f"{'NMD+ (True)':25} {result['nmd_true_with']:>18,} "
                     f"{result['nmd_true_without']:>18,} {nmd_true_total:>10,}")
        lines.append(f"{'NMD- (False)':25} {result['nmd_false_with']:>18,} "
                     f"{result['nmd_false_without']:>18,} {nmd_false_total:>10,}")
        lines.append(f"{'No prediction':25} {result['no_pred_with']:>18,} "
                     f"{result['no_pred_without']:>18,} {no_pred_total:>10,}")
        lines.append(f"{'Total':25} {result['count_with_introners']:>18,} "
                     f"{result['count_without_introners']:>18,} {total:>10,}")
        lines.append("")

        fisher = result['fisher']
        lines.append("Statistical Test (Fisher's Exact Test - NMD+ vs NMD-):")
        lines.append(f"  Odds Ratio: {fisher['odds_ratio']:.3f}")
        if fisher['ci_lower'] is not None:
            lines.append(f"  95% CI: [{fisher['ci_lower']:.3f}, {fisher['ci_upper']:.3f}]")
        lines.append(f"  p-value: {fisher['p_value']:.4e}")

        if fisher['odds_ratio'] > 1 and fisher['p_value'] < 0.05:
            lines.append(f"  -> NMD+ transcripts are {fisher['odds_ratio']:.2f}x "
                         "more likely to contain introners (significant)")
        elif fisher['odds_ratio'] < 1 and fisher['p_value'] < 0.05:
            lines.append(f"  -> NMD+ transcripts are {1/fisher['odds_ratio']:.2f}x "
                         "less likely to contain introners (significant)")
        else:
            lines.append("  -> No significant association")
        lines.append("")

    if len(results) > 1:
        lines.append("Summary")
        lines.append("-" * 70)
        lines.append(f"Total samples analyzed: {len(results)}")
        lines.append("")
        avg_total = sum(r['total_multi_exon'] for r in results) / len(results)
        avg_with_pred = sum(
            (r['count_with_prediction'] / r['total_multi_exon'] * 100)
            if r['total_multi_exon'] > 0 else 0
            for r in results
        ) / len(results)
        avg_nmd_pos_all = sum(
            (r['count_nmd_positive'] / r['total_multi_exon'] * 100)
            if r['total_multi_exon'] > 0 else 0
            for r in results
        ) / len(results)
        avg_nmd_pos_pred = sum(
            (r['count_nmd_positive'] / r['count_with_prediction'] * 100)
            if r['count_with_prediction'] > 0 else 0
            for r in results
        ) / len(results)
        avg_with_intr = sum(
            (r['count_with_introners'] / r['total_multi_exon'] * 100)
            if r['total_multi_exon'] > 0 else 0
            for r in results
        ) / len(results)
        lines.append(f"Average multi-exon transcripts per sample: {avg_total:,.0f}")
        lines.append(f"Average % with NMD predictions: {avg_with_pred:.1f}%")
        lines.append(f"Average % predicted as NMD+ (of all multi-exon): {avg_nmd_pos_all:.1f}%")
        lines.append(f"Average % predicted as NMD+ (of those with predictions): {avg_nmd_pos_pred:.1f}%")
        lines.append(f"Average % with introners: {avg_with_intr:.1f}%")
        lines.append("")

    lines.append("Notes:")
    lines.append("-" * 70)
    lines.append("- Data filtered for high-confidence isoforms only:")
    lines.append("  * Structural categories: full-splice_match, novel_in_catalog, novel_not_in_catalog")
    lines.append("- Mating-type region genes excluded before multi-exon filtering.")
    lines.append("- Multi-exon transcripts without NMD predictions are typically non-coding transcripts.")
    lines.append("- Introners are detected by checking whether strain-specific introner BED loci overlap gene exon coordinates.")

    return "\n".join(lines)


def analyze_data(df, gene_exons, introners_by_chrom, categories,
                 mating_type_genes=None):
    """Core NMD x introner analysis on a filtered DataFrame."""
    # Filter by structural category
    n_initial = len(df)
    df = df[df['structural_category'].isin(categories)].copy()
    n_after_category = len(df)

    # Filter mating-type region genes after structural category filtering to
    # preserve the original analysis order.
    mating_type_genes = mating_type_genes or set()
    gene_col_for_mt = 'associated_gene' if 'associated_gene' in df.columns else 'gene_id'
    if mating_type_genes and gene_col_for_mt in df.columns:
        df = df[~df[gene_col_for_mt].isin(mating_type_genes)].copy()
    n_after_mt = len(df)

    # Filter to multi-exon
    multi_exon = df[df['exons'] > 1].copy()
    total_multi_exon = len(multi_exon)

    # NMD counts
    count_with_prediction = multi_exon['predicted_NMD'].isin([True, False]).sum()
    count_nmd_positive = (multi_exon['predicted_NMD'] == True).sum()
    count_nmd_negative = (multi_exon['predicted_NMD'] == False).sum()

    # Introner overlap
    multi_exon['has_introner'] = False
    multi_exon['gene_found'] = False
    gene_col = 'associated_gene' if 'associated_gene' in multi_exon.columns else 'gene_id'
    not_found_count = 0

    for idx, row in multi_exon.iterrows():
        candidates = clean_gene_name(row.get(gene_col, ''))
        gene_found = False
        for candidate in candidates:
            if candidate in gene_exons:
                gene_found = True
                multi_exon.at[idx, 'gene_found'] = True
                if has_introner_overlap(gene_exons[candidate], introners_by_chrom):
                    multi_exon.at[idx, 'has_introner'] = True
                break
        if not gene_found:
            not_found_count += 1

    count_with = multi_exon['has_introner'].sum()
    count_without = total_multi_exon - count_with

    # Cross-tabulation
    nmd_true_with = len(multi_exon[(multi_exon['predicted_NMD'] == True) &
                                    (multi_exon['has_introner'] == True)])
    nmd_true_without = len(multi_exon[(multi_exon['predicted_NMD'] == True) &
                                       (multi_exon['has_introner'] == False)])
    nmd_false_with = len(multi_exon[(multi_exon['predicted_NMD'] == False) &
                                     (multi_exon['has_introner'] == True)])
    nmd_false_without = len(multi_exon[(multi_exon['predicted_NMD'] == False) &
                                        (multi_exon['has_introner'] == False)])
    no_pred_with = len(multi_exon[(multi_exon['predicted_NMD'].isna()) &
                                  (multi_exon['has_introner'] == True)])
    no_pred_without = len(multi_exon[(multi_exon['predicted_NMD'].isna()) &
                                     (multi_exon['has_introner'] == False)])

    fisher = calculate_fisher_test(nmd_true_with, nmd_true_without,
                                    nmd_false_with, nmd_false_without)

    return {
        'n_initial': n_initial,
        'n_after_category': n_after_category,
        'n_after_mt': n_after_mt,
        'total_multi_exon': total_multi_exon,
        'count_with_prediction': count_with_prediction,
        'count_nmd_positive': count_nmd_positive,
        'count_nmd_negative': count_nmd_negative,
        'count_with_introners': count_with,
        'count_without_introners': count_without,
        'nmd_true_with': nmd_true_with,
        'nmd_true_without': nmd_true_without,
        'nmd_false_with': nmd_false_with,
        'nmd_false_without': nmd_false_without,
        'no_pred_with': no_pred_with,
        'no_pred_without': no_pred_without,
        'not_found_count': not_found_count,
        'fisher': fisher,
        'multi_exon_df': multi_exon
    }


# ============================================================
# MODE IMPLEMENTATIONS
# ============================================================

def run_default_mode(args):
    """Default mode: analyze merged SQANTI3 data with single GTF + BED."""
    categories = set(args.categories.split(','))

    df = pd.read_csv(args.sqanti, sep='\t')
    gene_exons = parse_gtf_file(args.gtf)
    introners_by_chrom = parse_bed_file(args.introner_bed)
    mating_type_genes = extract_mating_type_genes(args.mt_gtf) if args.mt_gtf else set()

    result = analyze_data(df, gene_exons, introners_by_chrom, categories,
                          mating_type_genes=mating_type_genes)
    result['label'] = 'All strains'
    result['genome'] = 'merged'
    result['categories'] = sorted(categories)

    # Write text report
    report = format_output([result])
    with open(args.output, 'w') as f:
        f.write(report)

    # Write contingency CSV
    contingency = pd.DataFrame({
        'category': ['NMD+', 'NMD-'],
        'with_introners': [result['nmd_true_with'], result['nmd_false_with']],
        'without_introners': [result['nmd_true_without'], result['nmd_false_without']]
    })
    contingency.to_csv(args.contingency, index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(2)
    width = 0.35
    ax.bar(x - width/2,
           [result['nmd_true_with'], result['nmd_false_with']],
           width, label='With Introners', color='#e74c3c', alpha=0.8)
    ax.bar(x + width/2,
           [result['nmd_true_without'], result['nmd_false_without']],
           width, label='Without Introners', color='#3498db', alpha=0.8)
    ax.set_xlabel('NMD Prediction', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('NMD Status by Introner Presence', fontsize=14,
                 fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['NMD+', 'NMD-'])
    ax.legend()

    p = result['fisher']['p_value']
    or_val = result['fisher']['odds_ratio']
    ax.text(0.95, 0.95, f"OR={or_val:.2f}\np={p:.2e}",
            transform=ax.transAxes, ha='right', va='top',
            fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat'))

    fig.tight_layout()
    save_figures_with_png([fig], args.plot)
    plt.close(fig)

    print(f"NMD analysis complete. Report: {args.output}")


def run_by_strain_mode(args):
    """By-strain mode: analyze each strain with its own GTF."""
    categories = {'full-splice_match', 'novel_in_catalog', 'novel_not_in_catalog'}

    strain_inputs = {
        'CCMP1545': {
            'label': '834',
            'sqanti': args.sqanti_ccmp1545 or args.sqanti,
            'gtf': args.gtf_ccmp1545,
            'bed': args.bed_ccmp1545,
        },
        'RCC1614': {
            'label': '1614',
            'sqanti': args.sqanti_rcc1614 or args.sqanti,
            'gtf': args.gtf_rcc1614,
            'bed': args.bed_rcc1614,
        },
        'RCC1749': {
            'label': '1749',
            'sqanti': args.sqanti_rcc1749 or args.sqanti,
            'gtf': args.gtf_rcc1749,
            'bed': args.bed_rcc1749,
        },
    }
    mating_type_genes = extract_mating_type_genes(args.mt_gtf or args.gtf_ccmp1545)

    results = []
    for strain, paths in strain_inputs.items():
        strain_df = pd.read_csv(paths['sqanti'], sep='\t')
        if 'strain' in strain_df.columns:
            strain_df = strain_df[strain_df['strain'] == strain].copy()

        gene_exons = parse_gtf_file(paths['gtf'])
        introners_by_chrom = parse_bed_file(paths['bed'])

        result = analyze_data(
            strain_df,
            gene_exons,
            introners_by_chrom,
            categories,
            mating_type_genes=mating_type_genes,
        )
        result['label'] = paths['label']
        result['genome'] = strain
        result['categories'] = sorted(categories)
        results.append(result)

    # Write combined report
    report = format_output(results)
    with open(args.output, 'w') as f:
        f.write(report)

    if args.contingency:
        rows = []
        for result in results:
            rows.extend([
                {
                    'sample': result['label'],
                    'genome': result['genome'],
                    'category': 'NMD+',
                    'with_introners': result['nmd_true_with'],
                    'without_introners': result['nmd_true_without'],
                },
                {
                    'sample': result['label'],
                    'genome': result['genome'],
                    'category': 'NMD-',
                    'with_introners': result['nmd_false_with'],
                    'without_introners': result['nmd_false_without'],
                },
                {
                    'sample': result['label'],
                    'genome': result['genome'],
                    'category': 'No prediction',
                    'with_introners': result['no_pred_with'],
                    'without_introners': result['no_pred_without'],
                },
            ])
        pd.DataFrame(rows).to_csv(args.contingency, index=False)

    figures = []

    fig, axes = plt.subplots(1, len(results), figsize=(5 * len(results), 5),
                             sharey=False)
    if len(results) == 1:
        axes = [axes]
    for ax, result in zip(axes, results):
        x = np.arange(2)
        width = 0.35
        with_intr = [result['nmd_true_with'], result['nmd_false_with']]
        without_intr = [result['nmd_true_without'], result['nmd_false_without']]
        ax.bar(x - width/2, with_intr, width, label='With introners',
               color='#7DBE9C', edgecolor='black')
        ax.bar(x + width/2, without_intr, width, label='Without introners',
               color='#9DA8C8', edgecolor='black')
        ax.set_title(f"{result['label']} ({result['genome']})")
        ax.set_xticks(x)
        ax.set_xticklabels(['NMD+', 'NMD-'])
        ax.set_ylabel('Transcript count')
        fisher = result['fisher']
        ax.text(0.97, 0.97,
                f"OR={fisher['odds_ratio']:.2f}\np={fisher['p_value']:.2e}",
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85,
                          edgecolor='0.7'))
    axes[0].legend(frameon=False)
    fig.tight_layout()
    figures.append(fig)

    fig2, ax = plt.subplots(figsize=(8, 5))
    strains = [r['label'] for r in results]
    odds = [r['fisher']['odds_ratio'] for r in results]
    colors = ['#7DBE9C' if r['fisher']['p_value'] < 0.05 else '#B8B8B8'
              for r in results]
    ax.bar(strains, odds, color=colors, edgecolor='black')
    ax.axhline(1.0, color='black', linewidth=1, linestyle='--')
    ax.set_ylabel('Odds ratio')
    ax.set_xlabel('Sample')
    ax.set_title('NMD+ association with introner-containing genes')
    for i, result in enumerate(results):
        ax.text(i, odds[i], f"p={result['fisher']['p_value']:.1e}",
                ha='center', va='bottom', fontsize=9)
    fig2.tight_layout()
    figures.append(fig2)

    save_figures_with_png(figures, args.plot)
    for fig in figures:
        plt.close(fig)

    print(f"By-strain NMD analysis complete. Report: {args.output}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze NMD predictions and introner presence"
    )
    parser.add_argument('--mode', default='default',
                        choices=['default', 'by_strain'],
                        help='Analysis mode')
    parser.add_argument('--sqanti',
                        help='Merged SQANTI3 TSV')
    parser.add_argument('--output', required=True,
                        help='Output text report')
    parser.add_argument('--plot', required=True,
                        help='Output PDF plot (also writes sibling PNG)')

    # Default mode args
    parser.add_argument('--gtf', help='GTF file (default mode)')
    parser.add_argument('--introner_bed', help='BED file of introner loci')
    parser.add_argument('--contingency', help='Output contingency CSV')
    parser.add_argument('--categories',
                        default='full-splice_match,novel_in_catalog,novel_not_in_catalog',
                        help='Comma-separated structural categories')
    parser.add_argument('--mt_gtf',
                        help='GTF used to identify mating-type-region genes')

    # By-strain mode args
    parser.add_argument('--sqanti_ccmp1545', help='CCMP1545/834 SQANTI3 classification file')
    parser.add_argument('--sqanti_rcc1614', help='RCC1614/1614 SQANTI3 classification file')
    parser.add_argument('--sqanti_rcc1749', help='RCC1749/1749 SQANTI3 classification file')
    parser.add_argument('--gtf_ccmp1545', help='CCMP1545 GTF')
    parser.add_argument('--gtf_rcc1614', help='RCC1614 GTF')
    parser.add_argument('--gtf_rcc1749', help='RCC1749 GTF')
    parser.add_argument('--bed_ccmp1545', help='CCMP1545 introner BED')
    parser.add_argument('--bed_rcc1614', help='RCC1614 introner BED')
    parser.add_argument('--bed_rcc1749', help='RCC1749 introner BED')

    args = parser.parse_args()

    if args.mode == 'default':
        if not all([args.sqanti, args.gtf, args.introner_bed, args.contingency]):
            parser.error("Default mode requires --sqanti, --gtf, --introner_bed, --contingency")
        run_default_mode(args)
    elif args.mode == 'by_strain':
        required = [
            args.gtf_ccmp1545, args.gtf_rcc1614, args.gtf_rcc1749,
            args.bed_ccmp1545, args.bed_rcc1614, args.bed_rcc1749,
        ]
        has_split_sqanti = all([
            args.sqanti_ccmp1545, args.sqanti_rcc1614, args.sqanti_rcc1749,
        ])
        if not all(required) or not (has_split_sqanti or args.sqanti):
            parser.error(
                "by_strain mode requires per-strain GTFs, per-strain BEDs, "
                "and either --sqanti or all three --sqanti_* files"
            )
        run_by_strain_mode(args)


if __name__ == '__main__':
    main()

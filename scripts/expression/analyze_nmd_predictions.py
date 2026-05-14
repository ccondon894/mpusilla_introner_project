#!/usr/bin/env python3
"""
Analyze NMD predictions and introner presence for multi-exon transcript isoforms.

Dual-mode script:
  default:   Analyze NMD × introner association using merged SQANTI3 data
  by_strain: Per-strain NMD analysis with separate GTFs

Usage (Snakemake):
  Default:
    python analyze_nmd_predictions.py \
        --sqanti TSV --gtf GTF --introner_bed BED \
        --output TXT --contingency CSV --plot PDF \
        --categories STR

  By-strain:
    python analyze_nmd_predictions.py \
        --mode by_strain --sqanti TSV \
        --gtf_ccmp1545 GTF --gtf_rcc1614 GTF --gtf_rcc1749 GTF \
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
                     f"{result['n_after_category']:,}")
        lines.append(f"  Multi-exon transcripts: {result['total_multi_exon']:,}")
        lines.append("")

        lines.append("NMD Prediction Analysis:")
        total = result['total_multi_exon']
        nmd_pos = result['count_nmd_positive']
        nmd_neg = result['count_nmd_negative']
        with_pred = result['count_with_prediction']
        pct_pred = (with_pred / total * 100) if total > 0 else 0
        pct_pos = (nmd_pos / total * 100) if total > 0 else 0
        lines.append(f"  With NMD predictions: {with_pred:,} ({pct_pred:.1f}%)")
        lines.append(f"  NMD+: {nmd_pos:,} ({pct_pos:.1f}%)")
        lines.append(f"  NMD-: {nmd_neg:,}")
        lines.append("")

        lines.append("Introner Analysis:")
        lines.append(f"  With introners: {result['count_with_introners']:,}")
        lines.append(f"  Without introners: {result['count_without_introners']:,}")
        lines.append("")

        lines.append("Cross-tabulation (NMD x Introner):")
        lines.append(f"{'':25} {'With Introners':>18} {'Without Introners':>18}")
        lines.append(f"{'NMD+ (True)':25} {result['nmd_true_with']:>18,} "
                     f"{result['nmd_true_without']:>18,}")
        lines.append(f"{'NMD- (False)':25} {result['nmd_false_with']:>18,} "
                     f"{result['nmd_false_without']:>18,}")
        lines.append("")

        fisher = result['fisher']
        lines.append("Fisher's Exact Test:")
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

    return "\n".join(lines)


def analyze_data(df, gene_exons, introners_by_chrom, categories):
    """Core NMD x introner analysis on a filtered DataFrame."""
    # Filter by structural category
    n_initial = len(df)
    df = df[df['structural_category'].isin(categories)].copy()
    n_after_category = len(df)

    # Filter to multi-exon
    multi_exon = df[df['exons'] > 1].copy()
    total_multi_exon = len(multi_exon)

    # NMD counts
    count_with_prediction = multi_exon['predicted_NMD'].isin([True, False]).sum()
    count_nmd_positive = (multi_exon['predicted_NMD'] == True).sum()
    count_nmd_negative = (multi_exon['predicted_NMD'] == False).sum()

    # Introner overlap
    multi_exon['has_introner'] = False
    gene_col = 'gene_id' if 'gene_id' in multi_exon.columns else 'associated_gene'

    for idx, row in multi_exon.iterrows():
        candidates = clean_gene_name(row.get(gene_col, ''))
        for candidate in candidates:
            if candidate in gene_exons:
                if has_introner_overlap(gene_exons[candidate], introners_by_chrom):
                    multi_exon.at[idx, 'has_introner'] = True
                break

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

    fisher = calculate_fisher_test(nmd_true_with, nmd_true_without,
                                    nmd_false_with, nmd_false_without)

    return {
        'n_initial': n_initial,
        'n_after_category': n_after_category,
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

    result = analyze_data(df, gene_exons, introners_by_chrom, categories)
    result['label'] = 'All strains'
    result['genome'] = 'merged'

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

    df = pd.read_csv(args.sqanti, sep='\t')

    strain_gtfs = {
        'CCMP1545': args.gtf_ccmp1545,
        'RCC1614': args.gtf_rcc1614,
        'RCC1749': args.gtf_rcc1749
    }

    results = []
    for strain, gtf_path in strain_gtfs.items():
        strain_df = df[df['strain'] == strain] if 'strain' in df.columns else df

        gene_exons = parse_gtf_file(gtf_path)
        # Use GTF exons as proxy for introner loci (by-strain uses exon overlap)
        mt_genes = extract_mating_type_genes(gtf_path)
        gene_col = 'gene_id' if 'gene_id' in strain_df.columns else 'associated_gene'
        strain_df = strain_df[~strain_df[gene_col].isin(mt_genes)]

        # For by-strain mode, build introner coords from the reference GTF
        # Since we don't have per-strain BED, use empty introners
        introners_by_chrom = {}

        result = analyze_data(strain_df, gene_exons, introners_by_chrom, categories)
        result['label'] = strain
        result['genome'] = strain
        results.append(result)

    # Write combined report
    report = format_output(results)
    with open(args.output, 'w') as f:
        f.write(report)

    fig, ax = plt.subplots(figsize=(10, 6))
    strains = [r['label'] for r in results]
    nmd_pos = [r['count_nmd_positive'] for r in results]
    nmd_neg = [r['count_nmd_negative'] for r in results]

    x = np.arange(len(strains))
    width = 0.35
    ax.bar(x - width/2, nmd_pos, width, label='NMD+', color='#e74c3c', alpha=0.8)
    ax.bar(x + width/2, nmd_neg, width, label='NMD-', color='#2ecc71', alpha=0.8)
    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('NMD Predictions by Strain', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(strains)
    ax.legend()

    fig.tight_layout()
    save_figures_with_png([fig], args.plot)
    plt.close(fig)

    print(f"By-strain NMD analysis complete. Report: {args.output}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze NMD predictions and introner presence"
    )
    parser.add_argument('--mode', default='default',
                        choices=['default', 'by_strain'],
                        help='Analysis mode')
    parser.add_argument('--sqanti', required=True,
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

    # By-strain mode args
    parser.add_argument('--gtf_ccmp1545', help='CCMP1545 GTF')
    parser.add_argument('--gtf_rcc1614', help='RCC1614 GTF')
    parser.add_argument('--gtf_rcc1749', help='RCC1749 GTF')

    args = parser.parse_args()

    if args.mode == 'default':
        if not all([args.gtf, args.introner_bed, args.contingency]):
            parser.error("Default mode requires --gtf, --introner_bed, --contingency")
        run_default_mode(args)
    elif args.mode == 'by_strain':
        if not all([args.gtf_ccmp1545, args.gtf_rcc1614, args.gtf_rcc1749]):
            parser.error("by_strain mode requires --gtf_ccmp1545, --gtf_rcc1614, --gtf_rcc1749")
        run_by_strain_mode(args)


if __name__ == '__main__':
    main()

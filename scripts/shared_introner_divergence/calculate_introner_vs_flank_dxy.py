#!/usr/bin/env python3
"""
Calculate Dxy on introner body vs flanking regions for shared loci and
test whether introner divergence matches flanking divergence (ancestral)
or differs significantly (potential convergent insertion).

Produces per-category TSV files, a summary report, and comparison plots
stratified by family concordance (same vs different introner family between groups).
"""

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from scipy.stats import wilcoxon, mannwhitneyu, binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Compare introner body Dxy to flanking Dxy for shared loci')
    parser.add_argument('--alignment-dir', required=True)
    parser.add_argument('--classification-json', required=True)
    parser.add_argument('--genotype-matrix', required=True,
                       help='Genotype matrix TSV for family annotations')
    parser.add_argument('--group1-samples', nargs='+', required=True)
    parser.add_argument('--group2-samples', nargs='+', required=True)
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for TSV files')
    parser.add_argument('--figures-dir', required=True,
                       help='Output directory for plots')
    return parser.parse_args()


def is_valid_base(base):
    return base.upper() not in ['-', 'N']


def calculate_pairwise_counts(seq1, seq2):
    """Calculate pairwise differences and valid sites between two sequences."""
    differences = 0
    sites = 0
    for a, b in zip(seq1, seq2):
        if is_valid_base(a) and is_valid_base(b):
            sites += 1
            if a.upper() != b.upper():
                differences += 1
    return differences, sites


def calculate_dxy_ros(seqs_group1, seqs_group2):
    """Calculate between-group Dxy using ratio of sums approach."""
    if not seqs_group1 or not seqs_group2:
        return None

    total_diff = 0
    total_sites = 0

    for s1 in seqs_group1:
        for s2 in seqs_group2:
            diff, sites = calculate_pairwise_counts(s1, s2)
            total_diff += diff
            total_sites += sites

    return total_diff / total_sites if total_sites > 0 else None


def calculate_pi_ros(sequences):
    """Calculate within-group nucleotide diversity (pi) using ratio of sums."""
    if len(sequences) < 2:
        return None

    total_diff = 0
    total_sites = 0

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            diff, sites = calculate_pairwise_counts(sequences[i], sequences[j])
            total_diff += diff
            total_sites += sites

    return total_diff / total_sites if total_sites > 0 else None


def load_alignment_sequences(alignment_dir, ortholog_id, region, group1_samples, group2_samples):
    """Load aligned sequences for a locus and split by group."""
    aligned_file = os.path.join(alignment_dir, f"{ortholog_id}.{region}.mafft.fa")

    if not os.path.exists(aligned_file):
        return None, None

    g1_seqs = []
    g2_seqs = []

    for record in SeqIO.parse(aligned_file, 'fasta'):
        # MAFFT --adjustdirection may prepend _R_ to reversed sequences
        sample_name = record.id.replace('_R_', '')
        seq = str(record.seq)

        if sample_name in group1_samples:
            g1_seqs.append(seq)
        elif sample_name in group2_samples:
            g2_seqs.append(seq)

    return g1_seqs, g2_seqs


def annotate_families(classification, genotype_matrix_path, group1_samples, group2_samples):
    """Annotate each shared locus with family information per group."""
    df = pd.read_csv(genotype_matrix_path, sep='\t')
    group1_set = set(group1_samples)
    group2_set = set(group2_samples)
    shared_ids = set(classification['loci'].keys())

    present = df[(df['ortholog_id'].isin(shared_ids)) & (df['presence'] == 1)]

    family_info = {}
    for oid in shared_ids:
        rows = present[present['ortholog_id'] == oid]
        g1_fams = set(rows[rows['sample'].isin(group1_set)]['family'].unique()) - {-1}
        g2_fams = set(rows[rows['sample'].isin(group2_set)]['family'].unique()) - {-1}

        overlap = g1_fams & g2_fams
        same_family = bool(overlap)

        family_info[oid] = {
            'group1_families': ','.join(map(str, sorted(g1_fams))) if g1_fams else 'NA',
            'group2_families': ','.join(map(str, sorted(g2_fams))) if g2_fams else 'NA',
            'family_concordance': 'same_family' if same_family else 'different_family',
        }

    return family_info


def calculate_metrics_for_loci(alignment_dir, classification, group1_set, group2_set):
    """Calculate Dxy metrics for all shared loci."""
    results = []

    for oid, locus_info in sorted(classification['loci'].items()):
        g1_body, g2_body = load_alignment_sequences(
            alignment_dir, oid, 'introner_body', group1_set, group2_set)
        g1_left, g2_left = load_alignment_sequences(
            alignment_dir, oid, 'left_flank', group1_set, group2_set)
        g1_right, g2_right = load_alignment_sequences(
            alignment_dir, oid, 'right_flank', group1_set, group2_set)

        dxy_introner = calculate_dxy_ros(g1_body, g2_body) if g1_body and g2_body else None
        dxy_left = calculate_dxy_ros(g1_left, g2_left) if g1_left and g2_left else None
        dxy_right = calculate_dxy_ros(g1_right, g2_right) if g1_right and g2_right else None

        flank_values = [v for v in [dxy_left, dxy_right] if v is not None]
        dxy_flank_mean = np.mean(flank_values) if flank_values else None

        pi_g1_introner = calculate_pi_ros(g1_body) if g1_body and len(g1_body) >= 2 else None
        pi_g2_introner = calculate_pi_ros(g2_body) if g2_body and len(g2_body) >= 2 else None

        results.append({
            'ortholog_id': oid,
            'category': locus_info['category'],
            'cross_group_status': locus_info.get('cross_group_status', 'NA'),
            'ancestry_class': locus_info.get('ancestry_class', 'unclassified'),
            'within_group_status': locus_info.get('within_group_status', ''),
            'group1_present_count': locus_info['group1_present_count'],
            'group2_present_count': locus_info['group2_present_count'],
            'dxy_introner': dxy_introner,
            'dxy_left_flank': dxy_left,
            'dxy_right_flank': dxy_right,
            'dxy_flank_mean': dxy_flank_mean,
            'pi_group1_introner': pi_g1_introner,
            'pi_group2_introner': pi_g2_introner,
            'n_group1_seqs': len(g1_body) if g1_body else 0,
            'n_group2_seqs': len(g2_body) if g2_body else 0,
        })

    return pd.DataFrame(results)


def generate_plots(df, category_name, figures_dir):
    """Generate comparison plots for a category, stratified by family concordance."""
    plot_df = df.dropna(subset=['dxy_introner', 'dxy_flank_mean'])

    if len(plot_df) < 3:
        print(f"  Skipping plots for {category_name}: only {len(plot_df)} loci with complete data")
        return

    introner_dxy = plot_df['dxy_introner'].values
    flank_dxy = plot_df['dxy_flank_mean'].values

    # Wilcoxon signed-rank test (paired, all loci)
    try:
        stat, p_value = wilcoxon(introner_dxy, flank_dxy)
    except ValueError:
        stat, p_value = None, None

    # --- Figure 1: Scatter plot colored by family concordance ---
    fig, ax = plt.subplots(figsize=(9, 8))

    same_mask = plot_df['family_concordance'] == 'same_family'
    diff_mask = plot_df['family_concordance'] == 'different_family'

    n_same = same_mask.sum()
    n_diff = diff_mask.sum()

    if same_mask.any():
        ax.scatter(flank_dxy[same_mask], introner_dxy[same_mask],
                  alpha=0.6, s=45, edgecolors='black', linewidth=0.5,
                  color='#2ca02c', label=f'Same family (n={n_same})', zorder=3)
    if diff_mask.any():
        ax.scatter(flank_dxy[diff_mask], introner_dxy[diff_mask],
                  alpha=0.6, s=45, edgecolors='black', linewidth=0.5,
                  color='#d62728', label=f'Different family (n={n_diff})', zorder=3)

    max_val = max(max(flank_dxy), max(introner_dxy)) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1, label='y = x')

    ax.set_xlabel('Flanking Region Dxy', fontsize=13)
    ax.set_ylabel('Introner Body Dxy', fontsize=13)
    ax.set_title(f'Introner vs Flanking Divergence ({category_name})\n'
                f'n = {len(plot_df)} loci', fontsize=14)
    ax.set_xlim(-0.005, max_val)
    ax.set_ylim(-0.005, max_val)
    ax.set_aspect('equal')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(True, alpha=0.3)

    stats_text = (f'Median introner Dxy: {np.median(introner_dxy):.4f}\n'
                 f'Median flank Dxy: {np.median(flank_dxy):.4f}')
    if p_value is not None:
        stats_text += f'\nWilcoxon p = {p_value:.2e}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    scatter_path = os.path.join(figures_dir, f"dxy_scatter_{category_name}.png")
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()

    # --- Figure 2: Violin plot, 4 groups ---
    fig, ax = plt.subplots(figsize=(8, 8))

    groups = []
    labels = []
    colors = []

    for concordance, color, label_prefix in [
        ('same_family', '#2ca02c', 'Same family'),
        ('different_family', '#d62728', 'Diff. family')
    ]:
        mask = plot_df['family_concordance'] == concordance
        if mask.sum() < 2:
            continue
        subset = plot_df[mask]
        groups.append(subset['dxy_introner'].values)
        labels.append(f'{label_prefix}\nintroner')
        colors.append(color)
        groups.append(subset['dxy_flank_mean'].values)
        labels.append(f'{label_prefix}\nflank')
        colors.append(color)

    if groups:
        positions = list(range(1, len(groups) + 1))
        parts = ax.violinplot(groups, positions=positions, showmedians=True)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.6)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel('Dxy (between-group divergence)', fontsize=13)
        ax.set_title(f'Divergence by Family Concordance ({category_name})',
                    fontsize=14)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    violin_path = os.path.join(figures_dir, f"dxy_violin_{category_name}.png")
    plt.savefig(violin_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Plots saved: {scatter_path}, {violin_path}")

    # Print summary
    n_introner_higher = (introner_dxy > flank_dxy).sum()
    n_flank_higher = (flank_dxy > introner_dxy).sum()
    print(f"  Introner Dxy > Flank Dxy: {n_introner_higher}/{len(plot_df)} loci")
    print(f"  Flank Dxy > Introner Dxy: {n_flank_higher}/{len(plot_df)} loci")
    if p_value is not None:
        print(f"  Wilcoxon signed-rank test: p = {p_value:.2e}")


def calculate_family_concordance_chance(genotype_matrix_path, group1_samples, group2_samples):
    """
    Calculate the probability that a random introner from Group 1 and a random
    introner from Group 2 belong to the same family, based on each group's
    overall family frequency distribution.

    P(match by chance) = sum_f P(G1=f) * P(G2=f)
    """
    df = pd.read_csv(genotype_matrix_path, sep='\t')
    group1_set = set(group1_samples)
    group2_set = set(group2_samples)

    present = df[(df['presence'] == 1) & (df['family'] != -1)]

    # Count unique loci per family per group (one entry per ortholog_id)
    g1 = present[present['sample'].isin(group1_set)]
    g2 = present[present['sample'].isin(group2_set)]

    g1_fam_freq = g1.groupby('ortholog_id')['family'].first().value_counts(normalize=True)
    g2_fam_freq = g2.groupby('ortholog_id')['family'].first().value_counts(normalize=True)

    all_families = set(g1_fam_freq.index) | set(g2_fam_freq.index)
    p_match = sum(g1_fam_freq.get(f, 0) * g2_fam_freq.get(f, 0) for f in all_families)

    return p_match, g1_fam_freq, g2_fam_freq


def write_summary_report(all_metrics, output_path, genotype_matrix_path,
                         group1_samples, group2_samples):
    """Write a human-readable summary report of shared loci characteristics."""
    valid = all_metrics.dropna(subset=['dxy_introner', 'dxy_flank_mean'])

    with open(output_path, 'w') as f:
        f.write("=" * 72 + "\n")
        f.write("SHARED INTRONER DIVERGENCE ANALYSIS - SUMMARY REPORT\n")
        f.write("=" * 72 + "\n\n")

        # --- Section 1: Locus counts ---
        f.write("1. SHARED LOCUS COUNTS\n")
        f.write("-" * 40 + "\n")

        n_total = len(all_metrics)
        n_fixed = (all_metrics['category'] == 'fixed_shared').sum()
        n_poly = (all_metrics['category'] == 'polymorphic_shared').sum()
        n_same = (all_metrics['family_concordance'] == 'same_family').sum()
        n_diff = (all_metrics['family_concordance'] == 'different_family').sum()

        f.write(f"Total shared loci:                {n_total}\n")
        f.write(f"  Fixed in all samples:           {n_fixed}\n")
        f.write(f"  Polymorphic (not all samples):  {n_poly}\n")
        f.write(f"\n")
        f.write(f"Same introner family in both groups:      {n_same} ({100*n_same/n_total:.1f}%)\n")
        f.write(f"Different introner family between groups:  {n_diff} ({100*n_diff/n_total:.1f}%)\n")
        f.write(f"\n")

        # Cross-tabulation
        f.write("Cross-tabulation (fixation x family concordance):\n")
        for cat in ['fixed_shared', 'polymorphic_shared']:
            cat_df = all_metrics[all_metrics['category'] == cat]
            n_cat_same = (cat_df['family_concordance'] == 'same_family').sum()
            n_cat_diff = (cat_df['family_concordance'] == 'different_family').sum()
            f.write(f"  {cat}: {n_cat_same} same family, {n_cat_diff} different family\n")
        f.write(f"\n")

        # --- Section 2: Family pair breakdown ---
        f.write("2. FAMILY PAIR BREAKDOWN (different-family loci)\n")
        f.write("-" * 40 + "\n")
        diff_loci = all_metrics[all_metrics['family_concordance'] == 'different_family']
        pair_counts = diff_loci.groupby(['group1_families', 'group2_families']).size()
        pair_counts = pair_counts.sort_values(ascending=False)
        for (g1f, g2f), count in pair_counts.items():
            f.write(f"  Group1={g1f} vs Group2={g2f}: {count} loci\n")
        f.write(f"\n")

        # --- Section 3: Dxy comparison ---
        f.write("3. DXY COMPARISON (introner body vs flanking regions)\n")
        f.write("-" * 40 + "\n")

        for label, subset in [
            ("All shared loci", valid),
            ("Same family", valid[valid['family_concordance'] == 'same_family']),
            ("Different family", valid[valid['family_concordance'] == 'different_family']),
            ("Fixed shared", valid[valid['category'] == 'fixed_shared']),
            ("Polymorphic shared", valid[valid['category'] == 'polymorphic_shared']),
        ]:
            if len(subset) < 2:
                continue

            intr = subset['dxy_introner']
            flank = subset['dxy_flank_mean']

            try:
                _, p = wilcoxon(intr.values, flank.values)
                p_str = f"{p:.2e}"
            except ValueError:
                p_str = "NA"

            n_intr_higher = (intr.values > flank.values).sum()

            f.write(f"\n  {label} (n={len(subset)}):\n")
            f.write(f"    Introner Dxy:  median={intr.median():.4f}, mean={intr.mean():.4f}\n")
            f.write(f"    Flank Dxy:     median={flank.median():.4f}, mean={flank.mean():.4f}\n")
            f.write(f"    Introner > Flank: {n_intr_higher}/{len(subset)} loci\n")
            f.write(f"    Wilcoxon signed-rank p-value: {p_str}\n")

        f.write(f"\n")

        # --- Section 4: Family concordance vs chance ---
        f.write("4. FAMILY CONCORDANCE VS CHANCE EXPECTATION\n")
        f.write("-" * 40 + "\n")

        p_match, g1_freq, g2_freq = calculate_family_concordance_chance(
            genotype_matrix_path, group1_samples, group2_samples)

        f.write("Family frequency distributions (per unique locus):\n")
        f.write(f"  {'Family':>8s}  {'Group1':>8s}  {'Group2':>8s}\n")
        all_fams = sorted(set(g1_freq.index) | set(g2_freq.index))
        for fam in all_fams:
            f1 = g1_freq.get(fam, 0)
            f2 = g2_freq.get(fam, 0)
            f.write(f"  {fam:>8d}  {f1:>8.3f}  {f2:>8.3f}\n")
        f.write(f"\n")

        expected_same = p_match * n_total
        f.write(f"P(same family by chance) = {p_match:.4f} ({p_match*100:.1f}%)\n")
        f.write(f"Expected same-family matches by chance: {expected_same:.1f}/{n_total}\n")
        f.write(f"Observed same-family matches:           {n_same}/{n_total} ({100*n_same/n_total:.1f}%)\n")
        f.write(f"\n")

        binom_result = binomtest(n_same, n_total, p_match, alternative='greater')
        f.write(f"Binomial test (one-sided, observed > expected):\n")
        f.write(f"  p-value = {binom_result.pvalue:.4e}\n")
        f.write(f"\n")

        excess = n_same - expected_same
        f.write(f"Excess same-family loci beyond chance: ~{excess:.0f}\n")
        f.write(f"\n")

        # --- Section 5: Interpretation ---
        f.write("5. INTERPRETATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"Different-family loci ({n_diff}/{n_total}) are strong candidates\n")
        f.write(f"for convergent insertion: the same genomic position was\n")
        f.write(f"independently targeted by different introner families in\n")
        f.write(f"the two lineages.\n\n")

        f.write(f"Of the {n_same} same-family loci, ~{expected_same:.0f} are expected\n")
        f.write(f"by chance given each group's family frequency distribution.\n")
        f.write(f"The significant excess of ~{excess:.0f} same-family loci (binomial\n")
        f.write(f"p = {binom_result.pvalue:.2e}) suggests that a subset are\n")
        f.write(f"genuinely ancestral (inherited from the common ancestor).\n\n")

        f.write(f"Rough decomposition of {n_total} shared loci:\n")
        f.write(f"  ~{n_diff} different-family: clearly convergent\n")
        f.write(f"  ~{expected_same:.0f} same-family by chance: likely convergent\n")
        f.write(f"  ~{excess:.0f} same-family in excess: candidate ancestral\n")

        same_valid = valid[valid['family_concordance'] == 'same_family']
        if len(same_valid) > 0 and same_valid['dxy_flank_mean'].median() > 0:
            same_ratio = same_valid['dxy_introner'].median() / same_valid['dxy_flank_mean'].median()
            f.write(f"\nSame-family loci show introner Dxy {same_ratio:.1f}x higher than\n")
            f.write(f"flanking Dxy. This likely reflects relaxed selective constraint\n")
            f.write(f"on introner body sequences compared to (likely exonic)\n")
            f.write(f"flanking regions.\n")

        f.write("\n" + "=" * 72 + "\n")

    print(f"Summary report saved to: {output_path}")


def main():
    args = parse_arguments()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.figures_dir, exist_ok=True)

    print("Loading classification...")
    with open(args.classification_json) as f:
        classification = json.load(f)

    group1_set = set(args.group1_samples)
    group2_set = set(args.group2_samples)

    # Annotate families
    print("Annotating introner families per locus...")
    family_info = annotate_families(
        classification, args.genotype_matrix, args.group1_samples, args.group2_samples)

    print(f"Calculating Dxy metrics for {len(classification['loci'])} shared loci...")
    all_metrics = calculate_metrics_for_loci(
        args.alignment_dir, classification, group1_set, group2_set)

    # Add family annotations
    all_metrics['group1_families'] = all_metrics['ortholog_id'].map(
        lambda x: family_info[x]['group1_families'])
    all_metrics['group2_families'] = all_metrics['ortholog_id'].map(
        lambda x: family_info[x]['group2_families'])
    all_metrics['family_concordance'] = all_metrics['ortholog_id'].map(
        lambda x: family_info[x]['family_concordance'])

    # Write summary report
    summary_path = os.path.join(args.output_dir, "shared_introner_divergence_summary.txt")
    write_summary_report(all_metrics, summary_path, args.genotype_matrix,
                         args.group1_samples, args.group2_samples)

    # Process each category
    categories = {
        'all_shared': all_metrics,
        'fixed_shared': all_metrics[all_metrics['category'] == 'fixed_shared'],
        'polymorphic_shared': all_metrics[all_metrics['category'] == 'polymorphic_shared'],
    }

    for cat_name, cat_df in categories.items():
        print(f"\n{'='*60}")
        print(f"Category: {cat_name} ({len(cat_df)} loci)")
        print(f"{'='*60}")

        if len(cat_df) == 0:
            print("  No loci in this category, skipping.")
            continue

        # Save TSV
        tsv_path = os.path.join(args.output_dir, f"introner_vs_flank_dxy_{cat_name}.tsv")
        cat_df.to_csv(tsv_path, sep='\t', index=False)
        print(f"  Metrics saved to: {tsv_path}")

        # Summary stats
        valid = cat_df.dropna(subset=['dxy_introner', 'dxy_flank_mean'])
        print(f"  Loci with complete Dxy data: {len(valid)}")
        if len(valid) > 0:
            print(f"  Mean introner Dxy: {valid['dxy_introner'].mean():.4f}")
            print(f"  Mean flank Dxy: {valid['dxy_flank_mean'].mean():.4f}")

        # Generate plots
        generate_plots(cat_df, cat_name, args.figures_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()

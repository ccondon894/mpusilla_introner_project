#!/usr/bin/env python3
"""
Generate summary statistics for the final genotype matrix.

Reports ortholog group counts, derived per-group callability/pattern
breakdowns, legacy within- and cross-group status counts, gene occupancy,
family frequency composition, and shared locus statistics.
"""

import argparse
import sys
from collections import Counter
import pandas as pd
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate summary statistics for the genotype matrix')
    parser.add_argument('--genotype-matrix', required=True)
    parser.add_argument('--group1-samples', nargs='+', required=True)
    parser.add_argument('--group2-samples', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    return parser.parse_args()


def main():
    args = parse_arguments()

    df = pd.read_csv(args.genotype_matrix, sep='\t')
    group1_set = set(args.group1_samples)
    group2_set = set(args.group2_samples)
    all_loci = set(df['ortholog_id'].unique())

    present = df[df['presence'] == 1]
    g1_df = df[df['sample'].isin(group1_set)]
    g2_df = df[df['sample'].isin(group2_set)]

    # Loci with at least one present call per group
    g1_loci = set(present[present['sample'].isin(group1_set)]['ortholog_id'])
    g2_loci = set(present[present['sample'].isin(group2_set)]['ortholog_id'])
    shared_loci = g1_loci & g2_loci
    neither_loci = all_loci - g1_loci - g2_loci

    # Fixed shared: present in ALL samples of both groups
    fixed_shared = set()
    for oid in shared_loci:
        rows = df[df['ortholog_id'] == oid]
        g1_rows = rows[rows['sample'].isin(group1_set)]
        g2_rows = rows[rows['sample'].isin(group2_set)]
        if (g1_rows['presence'] == 1).all() and (g2_rows['presence'] == 1).all():
            fixed_shared.add(oid)

    # Classification breakdowns derived from the new status columns
    # One classification per ortholog group (collapse duplicate rows)
    group_status = df.drop_duplicates('ortholog_id').set_index('ortholog_id')
    total_groups = len(all_loci)
    within_status_counts = Counter(group_status['within_group_status'].fillna(''))
    cross_status_counts = Counter(group_status['cross_group_status'].fillna('NA'))
    # Treat empty string as NA for cross_group_status (annotate step may
    # replace literal 'NA' with empty)
    if '' in cross_status_counts:
        cross_status_counts['NA'] = cross_status_counts.pop('') + cross_status_counts.get('NA', 0)

    # For within-group-only groups (cross_group_status == NA), split by
    # which clade(s) have presence=1 members
    na_groups = set(group_status[
        group_status['cross_group_status'].fillna('NA').isin(['NA', ''])
    ].index)
    na_g1_only = len({oid for oid in na_groups
                      if oid in g1_loci and oid not in g2_loci})
    na_g2_only = len({oid for oid in na_groups
                      if oid in g2_loci and oid not in g1_loci})
    na_both_presence = len({oid for oid in na_groups
                            if oid in g1_loci and oid in g2_loci})
    na_neither = len({oid for oid in na_groups
                      if oid not in g1_loci and oid not in g2_loci})

    # Gene occupancy (per group, using present introners)
    def gene_occupancy(loci, present_df, group_samples):
        """% of loci that have a non-empty gene annotation in at least one sample."""
        group_present = present_df[present_df['sample'].isin(group_samples)]
        loci_genes = group_present[group_present['ortholog_id'].isin(loci)]
        has_gene = loci_genes.groupby('ortholog_id')['gene'].apply(
            lambda x: any(pd.notna(v) and v != '' for v in x))
        return has_gene.sum(), len(has_gene)

    g1_in_gene, g1_gene_total = gene_occupancy(g1_loci, present, group1_set)
    g2_in_gene, g2_gene_total = gene_occupancy(g2_loci, present, group2_set)

    # Family frequency composition (per group, one entry per locus using mode)
    def family_frequencies(loci, present_df, group_samples):
        """Family composition based on majority vote across samples per locus."""
        group_present = present_df[(present_df['sample'].isin(group_samples)) &
                                    (present_df['family'] != -1)]
        loci_fam = group_present[group_present['ortholog_id'].isin(loci)]
        # Majority family per locus
        locus_family = loci_fam.groupby('ortholog_id')['family'].agg(
            lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else x.iloc[0])
        return Counter(locus_family)

    g1_fam_counts = family_frequencies(g1_loci, present, group1_set)
    g2_fam_counts = family_frequencies(g2_loci, present, group2_set)

    has_derived = {
        'group1_pattern', 'group2_pattern',
        'group1_callability', 'group2_callability',
        'within_group_orthology_confidence',
        'cross_group_origin', 'cross_group_confidence',
    }.issubset(group_status.columns)

    # Write output
    with open(args.output, 'w') as f:
        f.write("=" * 64 + "\n")
        f.write("GENOTYPE MATRIX SUMMARY STATISTICS\n")
        f.write("=" * 64 + "\n\n")

        # Section 1: Ortholog group counts
        f.write("1. ORTHOLOG GROUP COUNTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total ortholog groups:   {len(all_loci):,}\n")
        f.write(f"Group 1 loci:            {len(g1_loci):,}\n")
        f.write(f"Group 2 loci:            {len(g2_loci):,}\n")
        f.write(f"Present in neither:      {len(neither_loci):,}\n")
        f.write(f"\n")

        # Section 2: Shared / fixed
        f.write("2. SHARED AND FIXED LOCI\n")
        f.write("-" * 40 + "\n")
        f.write(f"Shared (present in both groups):  {len(shared_loci):,} "
                f"({100*len(shared_loci)/len(all_loci):.1f}%)\n")
        f.write(f"Fixed shared (all samples):       {len(fixed_shared):,} "
                f"({100*len(fixed_shared)/len(all_loci):.1f}%)\n")
        f.write(f"\n")

        # Section 3: Within-group classification status
        f.write("3. DERIVED GROUP PATTERNS AND CONFIDENCE\n")
        f.write("-" * 40 + "\n")
        if has_derived:
            for label, col in [('Group 1 pattern', 'group1_pattern'),
                               ('Group 2 pattern', 'group2_pattern'),
                               ('Group 1 callability', 'group1_callability'),
                               ('Group 2 callability', 'group2_callability'),
                               ('Within-group orthology confidence',
                                'within_group_orthology_confidence'),
                               ('Cross-group origin', 'cross_group_origin'),
                               ('Cross-group confidence', 'cross_group_confidence')]:
                f.write(f"{label}:\n")
                counts = Counter(group_status[col].fillna('NA'))
                for status, count in sorted(counts.items(),
                                            key=lambda kv: (-kv[1], str(kv[0]))):
                    pct = 100 * count / total_groups if total_groups else 0
                    f.write(f"  {status:<26s} {count:>7,d}  ({pct:5.1f}%)\n")
                f.write("\n")
        else:
            f.write("Derived classification columns not present in matrix.\n\n")

        # Section 4: Legacy within-group classification status
        f.write("4. LEGACY WITHIN-GROUP STATUS\n")
        f.write("-" * 40 + "\n")
        within_order = ['consistent', 'singleton', 'low_identity',
                        'discordant', 'uncertain']
        for status in within_order:
            count = within_status_counts.get(status, 0)
            if count == 0 and status in ('discordant', 'uncertain'):
                continue  # these are filtered out upstream; skip zeros
            pct = 100 * count / total_groups if total_groups else 0
            f.write(f"  {status:<22s}  {count:>7,d}  ({pct:5.1f}%)\n")
        # Catch any unexpected values
        for status, count in within_status_counts.items():
            if status not in within_order:
                pct = 100 * count / total_groups if total_groups else 0
                f.write(f"  {status:<22s}  {count:>7,d}  ({pct:5.1f}%)\n")
        f.write(f"  {'TOTAL':<22s}  {total_groups:>7,d}\n")
        f.write(f"\n")

        # Section 5: Cross-group classification (ancestral vs independent)
        f.write("5. LEGACY CROSS-GROUP STATUS (ancestral vs independent)\n")
        f.write("-" * 40 + "\n")
        na_count = cross_status_counts.get('NA', 0)
        cross_group_total = total_groups - na_count
        f.write(f"Cross-group classification "
                f"(groups with both G1 and G2 members, n={cross_group_total:,}):\n")
        cross_order = ['ancestral', 'likely_ancestral', 'ancestral_low_identity',
                       'independent', 'likely_independent',
                       'uncertain']
        for status in cross_order:
            count = cross_status_counts.get(status, 0)
            if count == 0:
                continue
            pct = 100 * count / cross_group_total if cross_group_total else 0
            f.write(f"  {status:<24s}  {count:>5,d}  ({pct:5.1f}%)\n")
        # Catch any unexpected non-NA values
        for status, count in cross_status_counts.items():
            if status not in cross_order and status != 'NA':
                pct = 100 * count / cross_group_total if cross_group_total else 0
                f.write(f"  {status:<24s}  {count:>5,d}  ({pct:5.1f}%)\n")
        f.write(f"\n")
        f.write(f"Within-group-only groups (cross_group_status = NA): "
                f"{na_count:,}\n")
        f.write(f"  Present only in G1:   {na_g1_only:,}\n")
        f.write(f"  Present only in G2:   {na_g2_only:,}\n")
        if na_both_presence:
            f.write(f"  Present in both*:     {na_both_presence:,}\n")
            f.write(f"    * classification is NA despite presence in both\n"
                    f"      clades; likely lacked fingerprints for comparison\n")
        if na_neither:
            f.write(f"  Present in neither:   {na_neither:,}\n")
        f.write(f"\n")

        # Section 6: Gene occupancy
        f.write("6. GENE OCCUPANCY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Note: gene annotations from validation step; may be unreliable.\n")
        f.write(f"Group 1: {g1_in_gene:,}/{g1_gene_total:,} "
                f"({100*g1_in_gene/g1_gene_total:.1f}%) within annotated genes\n")
        f.write(f"Group 2: {g2_in_gene:,}/{g2_gene_total:,} "
                f"({100*g2_in_gene/g2_gene_total:.1f}%) within annotated genes\n")
        f.write(f"\n")

        # Section 7: Family frequency composition
        f.write("7. FAMILY FREQUENCY COMPOSITION\n")
        f.write("-" * 40 + "\n")
        all_families = sorted(set(g1_fam_counts.keys()) | set(g2_fam_counts.keys()))
        g1_fam_total = sum(g1_fam_counts.values())
        g2_fam_total = sum(g2_fam_counts.values())

        f.write(f"{'Family':>8s}  {'Group1':>8s} {'(%)':>7s}  {'Group2':>8s} {'(%)':>7s}\n")
        for fam in all_families:
            g1_c = g1_fam_counts.get(fam, 0)
            g2_c = g2_fam_counts.get(fam, 0)
            g1_pct = 100 * g1_c / g1_fam_total if g1_fam_total else 0
            g2_pct = 100 * g2_c / g2_fam_total if g2_fam_total else 0
            f.write(f"{int(fam):>8d}  {g1_c:>8,} {g1_pct:>6.1f}%  {g2_c:>8,} {g2_pct:>6.1f}%\n")
        f.write(f"{'Total':>8s}  {g1_fam_total:>8,}          {g2_fam_total:>8,}\n")

        f.write("\n" + "=" * 64 + "\n")

    print(f"Summary written to: {args.output}")


if __name__ == "__main__":
    main()

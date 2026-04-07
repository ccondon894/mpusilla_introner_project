#!/usr/bin/env python3
"""
Generate summary statistics for the final genotype matrix.

Reports ortholog group counts, cross-group presence/absence patterns,
gene occupancy, family frequency composition, and shared locus statistics.
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

    # Cross-group status
    def cross_group_status(source_loci, target_df, target_samples):
        """For each source locus, classify its status in the target group."""
        n_present = 0
        n_absent = 0
        n_not_callable = 0
        for oid in source_loci:
            target_rows = target_df[target_df['ortholog_id'] == oid]
            presences = set(target_rows['presence'])
            if 1 in presences:
                n_present += 1  # shared
            elif 2 in presences:
                n_absent += 1   # absent in at least one target sample
            else:
                n_not_callable += 1  # all target samples are 3
        return n_present, n_absent, n_not_callable

    g1_in_g2_present, g1_in_g2_absent, g1_in_g2_nocall = cross_group_status(g1_loci, g2_df, group2_set)
    g2_in_g1_present, g2_in_g1_absent, g2_in_g1_nocall = cross_group_status(g2_loci, g1_df, group1_set)

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

        # Section 3: Cross-group status
        f.write("3. CROSS-GROUP STATUS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Group 1 loci ({len(g1_loci):,}) in Group 2:\n")
        f.write(f"  Present (shared):    {g1_in_g2_present:,} ({100*g1_in_g2_present/len(g1_loci):.1f}%)\n")
        f.write(f"  Absent:              {g1_in_g2_absent:,} ({100*g1_in_g2_absent/len(g1_loci):.1f}%)\n")
        f.write(f"  Not callable:        {g1_in_g2_nocall:,} ({100*g1_in_g2_nocall/len(g1_loci):.1f}%)\n")
        f.write(f"\n")
        f.write(f"Group 2 loci ({len(g2_loci):,}) in Group 1:\n")
        f.write(f"  Present (shared):    {g2_in_g1_present:,} ({100*g2_in_g1_present/len(g2_loci):.1f}%)\n")
        f.write(f"  Absent:              {g2_in_g1_absent:,} ({100*g2_in_g1_absent/len(g2_loci):.1f}%)\n")
        f.write(f"  Not callable:        {g2_in_g1_nocall:,} ({100*g2_in_g1_nocall/len(g2_loci):.1f}%)\n")
        f.write(f"\n")

        # Section 4: Gene occupancy
        f.write("4. GENE OCCUPANCY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Note: gene annotations from validation step; may be unreliable.\n")
        f.write(f"Group 1: {g1_in_gene:,}/{g1_gene_total:,} "
                f"({100*g1_in_gene/g1_gene_total:.1f}%) within annotated genes\n")
        f.write(f"Group 2: {g2_in_gene:,}/{g2_gene_total:,} "
                f"({100*g2_in_gene/g2_gene_total:.1f}%) within annotated genes\n")
        f.write(f"\n")

        # Section 5: Family frequency composition
        f.write("5. FAMILY FREQUENCY COMPOSITION\n")
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

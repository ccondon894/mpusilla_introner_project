#!/usr/bin/env python3
"""
Test whether the observed rate of same-family introner matches among shared
loci (present in both Group 1 and Group 2) exceeds what would be expected
by chance, given each group's introner family frequency distribution.

P(same family by chance) = sum_f P(G1=f) * P(G2=f)

A significant excess of same-family matches suggests that a subset of shared
loci are genuinely ancestral (inherited from the common ancestor), rather
than all being convergent insertions that happen to share a family.
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import binomtest


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Test family concordance of shared introner loci against chance')
    parser.add_argument('--genotype-matrix', required=True,
                       help='Genotype matrix TSV')
    parser.add_argument('--group1-samples', nargs='+', required=True)
    parser.add_argument('--group2-samples', nargs='+', required=True)
    return parser.parse_args()


def main():
    args = parse_arguments()

    df = pd.read_csv(args.genotype_matrix, sep='\t')
    group1_set = set(args.group1_samples)
    group2_set = set(args.group2_samples)

    # ---------------------------------------------------------------
    # 1. Identify shared loci
    # ---------------------------------------------------------------
    present = df[df['presence'] == 1]
    g1_loci = set(present[present['sample'].isin(group1_set)]['ortholog_id'])
    g2_loci = set(present[present['sample'].isin(group2_set)]['ortholog_id'])
    shared_loci = g1_loci & g2_loci

    print(f"Group 1 loci (presence=1): {len(g1_loci)}")
    print(f"Group 2 loci (presence=1): {len(g2_loci)}")
    print(f"Shared loci:               {len(shared_loci)}")

    # ---------------------------------------------------------------
    # 2. Family frequency distributions (genome-wide, per unique locus)
    # ---------------------------------------------------------------
    present_no_neg = present[present['family'] != -1]

    g1_all = present_no_neg[present_no_neg['sample'].isin(group1_set)]
    g2_all = present_no_neg[present_no_neg['sample'].isin(group2_set)]

    # One entry per ortholog_id (take the first sample's family assignment)
    g1_fam_freq = g1_all.groupby('ortholog_id')['family'].first().value_counts(normalize=True)
    g2_fam_freq = g2_all.groupby('ortholog_id')['family'].first().value_counts(normalize=True)

    print(f"\nFamily frequency distributions (per unique locus):")
    print(f"  {'Family':>8s}  {'Group1':>8s}  {'Group2':>8s}")
    all_families = sorted(set(g1_fam_freq.index) | set(g2_fam_freq.index))
    for fam in all_families:
        f1 = g1_fam_freq.get(fam, 0)
        f2 = g2_fam_freq.get(fam, 0)
        print(f"  {fam:>8d}  {f1:>8.3f}  {f2:>8.3f}")

    # ---------------------------------------------------------------
    # 3. Expected P(same family by chance)
    # ---------------------------------------------------------------
    p_match = sum(g1_fam_freq.get(f, 0) * g2_fam_freq.get(f, 0) for f in all_families)
    print(f"\nP(same family by chance) = {p_match:.4f} ({p_match*100:.1f}%)")

    # ---------------------------------------------------------------
    # 4. Observed family concordance among shared loci
    # ---------------------------------------------------------------
    shared_df = present_no_neg[present_no_neg['ortholog_id'].isin(shared_loci)]

    n_same = 0
    n_diff = 0
    for oid in shared_loci:
        rows = shared_df[shared_df['ortholog_id'] == oid]
        g1_fams = set(rows[rows['sample'].isin(group1_set)]['family'].unique())
        g2_fams = set(rows[rows['sample'].isin(group2_set)]['family'].unique())
        if g1_fams & g2_fams:
            n_same += 1
        else:
            n_diff += 1

    n_total = n_same + n_diff
    print(f"\nObserved among {n_total} shared loci:")
    print(f"  Same family:      {n_same} ({100*n_same/n_total:.1f}%)")
    print(f"  Different family: {n_diff} ({100*n_diff/n_total:.1f}%)")

    # ---------------------------------------------------------------
    # 5. Binomial test
    # ---------------------------------------------------------------
    expected_same = p_match * n_total

    result_greater = binomtest(n_same, n_total, p_match, alternative='greater')
    result_two = binomtest(n_same, n_total, p_match, alternative='two-sided')

    print(f"\nExpected same-family by chance: {expected_same:.1f}/{n_total}")
    print(f"Observed same-family:          {n_same}/{n_total}")
    print(f"Excess beyond chance:          ~{n_same - expected_same:.0f}")

    print(f"\nBinomial test:")
    print(f"  One-sided (greater) p-value: {result_greater.pvalue:.4e}")
    print(f"  Two-sided p-value:           {result_two.pvalue:.4e}")

    # ---------------------------------------------------------------
    # 6. Decomposition
    # ---------------------------------------------------------------
    excess = n_same - expected_same
    print(f"\nRough decomposition of {n_total} shared loci:")
    print(f"  ~{n_diff} different-family loci: clearly convergent insertions")
    print(f"  ~{expected_same:.0f} same-family loci expected by chance: likely convergent")
    print(f"  ~{excess:.0f} same-family loci in excess of chance: candidate ancestral")


if __name__ == "__main__":
    main()

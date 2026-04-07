#!/usr/bin/env python3
"""
Count ortholog groups that meet criteria for RCC1749-polarized analysis.

Requirements for each ortholog group:
1. RCC1749 must have valid data (presence in [1,2], NOT 3)
2. ALL Group1 samples must have valid data (presence in [1,2], NOT 3)
3. Group1 must be polymorphic (at least one sample with 1 AND at least one with 2)
4. RCC3052 can have any value (1, 2, or 3/missing)

Breakdown by RCC1749 state:
- RCC1749 = 1 (present)
- RCC1749 = 2 (absent)
"""

import pandas as pd
import argparse

# Define Group1 samples (all samples except RCC1749 and RCC3052)
GROUP1_SAMPLES = {
    'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482',
    'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
}

def check_ortholog(orth_group, group1_samples=GROUP1_SAMPLES):
    """
    Check if an ortholog group meets all criteria.

    Args:
        orth_group: DataFrame subset for one ortholog_id
        group1_samples: Set of Group1 sample names

    Returns:
        tuple: (passes_criteria, rcc1749_state) where rcc1749_state is 1 or 2
    """

    # Get RCC1749 data
    rcc1749_rows = orth_group[orth_group['sample'] == 'RCC1749']
    if len(rcc1749_rows) == 0:
        return False, None  # No RCC1749 data

    rcc1749_presence = int(rcc1749_rows.iloc[0]['presence'])

    # RCC1749 must have valid data
    if rcc1749_presence not in [1, 2]:
        return False, None

    # Check Group1 samples
    group1_data = orth_group[orth_group['sample'].isin(group1_samples)]

    # All Group1 samples must be present
    if len(group1_data) != len(group1_samples):
        return False, None

    # All Group1 samples must have valid data
    presences = group1_data['presence'].values
    if not all(p in [1, 2] for p in presences):
        return False, None

    # Check for polymorphism in Group1
    has_present = 1 in presences
    has_absent = 2 in presences

    if not (has_present and has_absent):
        return False, None

    # All criteria met
    return True, rcc1749_presence


def main():
    parser = argparse.ArgumentParser(
        description="Count RCC1749-polarized ortholog groups"
    )
    parser.add_argument(
        '--matrix',
        required=True,
        help='Genotype matrix TSV file'
    )
    parser.add_argument(
        '--output',
        default='rcc1749_polarized_orthologs_summary.txt',
        help='Output summary file'
    )
    args = parser.parse_args()

    print(f"Reading genotype matrix from {args.matrix}...")
    df = pd.read_csv(args.matrix, sep='\t')

    total_orthologs = df['ortholog_id'].nunique()
    print(f"Total orthologs in matrix: {total_orthologs}")

    # Count useable orthologs
    useable_total = 0
    useable_rcc1749_present = 0
    useable_rcc1749_absent = 0

    print("\nProcessing orthologs...")
    for orth_id, group in df.groupby('ortholog_id'):
        passes, rcc1749_state = check_ortholog(group)

        if passes:
            useable_total += 1
            if rcc1749_state == 1:
                useable_rcc1749_present += 1
            elif rcc1749_state == 2:
                useable_rcc1749_absent += 1

    # Print summary
    print("\n" + "="*70)
    print("RESULTS: RCC1749-Polarized Ortholog Groups")
    print("="*70)
    print(f"\nTotal useable ortholog groups: {useable_total}")
    print(f"  - RCC1749 = present (1): {useable_rcc1749_present}")
    print(f"  - RCC1749 = absent (2):  {useable_rcc1749_absent}")
    print(f"\nFiltering criteria:")
    print(f"  ✓ RCC1749 has valid data (1 or 2)")
    print(f"  ✓ All 11 Group1 samples have valid data (1 or 2)")
    print(f"  ✓ Group1 is polymorphic (≥1 present AND ≥1 absent)")
    print(f"  ✓ RCC3052 allowed to be missing (3)")
    print("="*70)

    # Save to file
    with open(args.output, 'w') as f:
        f.write("RCC1749-Polarized Ortholog Groups Summary\n")
        f.write("="*70 + "\n\n")
        f.write(f"Total useable ortholog groups: {useable_total}\n")
        f.write(f"  - RCC1749 = present (1): {useable_rcc1749_present}\n")
        f.write(f"  - RCC1749 = absent (2):  {useable_rcc1749_absent}\n\n")
        f.write("Filtering criteria:\n")
        f.write("  ✓ RCC1749 has valid data (1 or 2)\n")
        f.write("  ✓ All 11 Group1 samples have valid data (1 or 2)\n")
        f.write("  ✓ Group1 is polymorphic (≥1 present AND ≥1 absent)\n")
        f.write("  ✓ RCC3052 allowed to be missing (3)\n")

    print(f"\nResults saved to {args.output}")


if __name__ == '__main__':
    main()

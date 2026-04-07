#!/usr/bin/env python3
"""
Classify sharing status for introner ortholog groups.

For each ortholog group in the genotype matrix, compares codon-level insertion
fingerprints among presence=1 members to distinguish:
  - ancestral:    same (codon, offset) + same family → shared before divergence
  - independent:  different (codon, offset) OR different family → separate insertions
  - ambiguous:    codon positions differ by ≤CODON_TOLERANCE, same family → unclear
  - uncertain:    insufficient fingerprint data to classify

Reads:
  - Oriented genotype matrix (TSV)
  - Per-sample insertion fingerprint files (TSV)

Outputs:
  - Verified genotype matrix with added 'sharing_status' column
  - Summary statistics (printed to stderr/stdout)
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


# Codon tolerance for cross-group comparisons
# Small offsets (1-3 codons) may arise from protein annotation differences
# between divergent genomes, not from different insertion sites
CODON_TOLERANCE = 3

GROUP2_SAMPLES = {'RCC1749', 'RCC3052'}


def load_fingerprints(fingerprint_files):
    """Load fingerprints from multiple per-sample TSV files.

    Returns dict keyed by (sample, contig, start, end).
    """
    fps = {}
    for fp_file in fingerprint_files:
        with open(fp_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                key = (row['sample'], row['contig'],
                       int(row['start']), int(row['end']))
                fps[key] = row
    return fps


def classify_ortholog_group(members, codon_tolerance=CODON_TOLERANCE):
    """Classify the sharing status of an ortholog group.

    Args:
        members: list of dicts, each with keys:
            sample, family, codon_number, codon_offset, confidence, group
        codon_tolerance: max codon difference to consider "same position"

    Returns:
        sharing_status string for the group
    """
    # Only compare presence=1 members with high-confidence fingerprints
    typed = [m for m in members if m['confidence'] == 'high']

    if len(typed) < 2:
        return 'uncertain'

    # Check if there are members from both groups
    groups_present = {m['group'] for m in typed}
    has_both_groups = 'G1' in groups_present and 'G2' in groups_present

    if not has_both_groups:
        # All typed members are from the same group — within-group comparison
        codon_nums = [int(m['codon_number']) for m in typed]
        spread = max(codon_nums) - min(codon_nums)

        if spread <= codon_tolerance:
            return 'consistent'
        else:
            return 'within_group_discordant'

    # Cross-group comparison
    g1_members = [m for m in typed if m['group'] == 'G1']
    g2_members = [m for m in typed if m['group'] == 'G2']

    if not g1_members or not g2_members:
        return 'uncertain'

    # Get representative codon positions for each group
    # Use the most common (codon, offset) within each group
    def get_consensus(group_members):
        """Get the most common (codon, offset, family) from a group."""
        from collections import Counter
        codon_counts = Counter(
            (int(m['codon_number']), int(m['codon_offset'])) for m in group_members)
        family_counts = Counter(m['family'] for m in group_members)
        return codon_counts.most_common(1)[0][0], family_counts.most_common(1)[0][0]

    g1_codon, g1_family = get_consensus(g1_members)
    g2_codon, g2_family = get_consensus(g2_members)

    g1_codon_num, g1_offset = g1_codon
    g2_codon_num, g2_offset = g2_codon

    codon_diff = abs(g1_codon_num - g2_codon_num)
    same_offset = (g1_offset == g2_offset)
    same_family = (g1_family == g2_family)

    # Classification logic
    if codon_diff == 0 and same_offset and same_family:
        return 'ancestral'

    if codon_diff == 0 and same_offset and not same_family:
        # Same position, different family — could be family misannotation
        # or convergent insertion at a hotspot
        return 'same_site_diff_family'

    if codon_diff <= codon_tolerance and same_offset and same_family:
        # Small positional difference likely from annotation imprecision
        return 'ambiguous'

    if codon_diff <= codon_tolerance and same_offset and not same_family:
        return 'ambiguous_diff_family'

    if not same_family:
        return 'independent'

    if codon_diff > codon_tolerance:
        return 'independent'

    if not same_offset and codon_diff <= codon_tolerance:
        # Different phase but close — likely independent
        return 'independent'

    return 'ambiguous'


def main():
    parser = argparse.ArgumentParser(
        description='Classify sharing status of introner ortholog groups')
    parser.add_argument('--matrix', required=True,
                        help='Oriented genotype matrix TSV')
    parser.add_argument('--fingerprints', required=True, nargs='+',
                        help='Per-sample fingerprint TSV files')
    parser.add_argument('--output', required=True,
                        help='Output verified genotype matrix TSV')
    parser.add_argument('--codon-tolerance', type=int, default=CODON_TOLERANCE,
                        help=f'Max codon difference for "same position" '
                             f'(default: {CODON_TOLERANCE})')
    parser.add_argument('--summary', default=None,
                        help='Optional summary output TSV')
    args = parser.parse_args()

    # Load fingerprints
    print(f"Loading fingerprints from {len(args.fingerprints)} files...")
    fps = load_fingerprints(args.fingerprints)
    print(f"  Loaded {len(fps)} fingerprint records")

    # Read genotype matrix and group by ortholog_id
    print("Reading genotype matrix...")
    matrix_rows = []
    ortholog_groups = defaultdict(list)

    with open(args.matrix) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            idx = len(matrix_rows)
            matrix_rows.append(row)
            ortholog_groups[row['ortholog_id']].append(idx)

    print(f"  {len(matrix_rows)} rows in {len(ortholog_groups)} ortholog groups")

    # Classify each ortholog group
    status_counts = defaultdict(int)
    group_summaries = []

    for oid, row_indices in ortholog_groups.items():
        # Build member list for classification
        members = []
        for idx in row_indices:
            row = matrix_rows[idx]
            if row['presence'] != '1':
                continue

            sample = row['sample']
            group = 'G2' if sample in GROUP2_SAMPLES else 'G1'

            # Look up fingerprint
            key = (sample, row['contig'], int(row['start']), int(row['end']))
            fp = fps.get(key, {})

            members.append({
                'sample': sample,
                'group': group,
                'family': row['family'],
                'codon_number': fp.get('codon_number', ''),
                'codon_offset': fp.get('codon_offset', ''),
                'confidence': fp.get('confidence', 'no_fingerprint'),
                'gene_id': fp.get('gene_id', ''),
            })

        # Classify
        status = classify_ortholog_group(members, args.codon_tolerance)
        status_counts[status] += 1

        # Store summary for this group
        if members:
            typed = [m for m in members if m['confidence'] == 'high']
            group_summaries.append({
                'ortholog_id': oid,
                'sharing_status': status,
                'n_present': len(members),
                'n_fingerprinted': len(typed),
                'families': ';'.join(sorted(set(m['family'] for m in members))),
                'codons': ';'.join(
                    f"{m['sample']}:{m['codon_number']}.{m['codon_offset']}"
                    for m in typed),
            })

        # Assign status to all rows in this group
        for idx in row_indices:
            matrix_rows[idx]['sharing_status'] = status

    # Print summary
    print(f"\nSharing status classification (codon_tolerance={args.codon_tolerance}):")
    for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")
    print(f"  {'TOTAL':30s}  {len(ortholog_groups):5d}")

    # Write output matrix
    out_fieldnames = fieldnames + ['sharing_status']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(matrix_rows)

    print(f"\nWrote {len(matrix_rows)} rows to {args.output}")

    # Write optional summary
    if args.summary:
        summary_fields = ['ortholog_id', 'sharing_status', 'n_present',
                          'n_fingerprinted', 'families', 'codons']
        with open(args.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(group_summaries)
        print(f"Wrote {len(group_summaries)} group summaries to {args.summary}")


if __name__ == '__main__':
    main()

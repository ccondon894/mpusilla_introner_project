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


def get_locus_key(member):
    """Build a comparison key for an introner based on its location type.

    For CDS introners:    ('cds', codon_number, codon_offset, exon_number)
    For intron introners: ('intron', intron_number)
    Returns None if the member has no usable location data.
    """
    loc = member.get('location_type', '')
    if loc == 'cds':
        try:
            return ('cds', int(member['codon_number']),
                    int(member['codon_offset']), int(member['exon_number']))
        except (ValueError, TypeError):
            return None
    if loc == 'intron':
        try:
            return ('intron', int(member['intron_number']))
        except (ValueError, TypeError):
            return None
    return None


def cluster_keys(keys, codon_tolerance):
    """Cluster a list of locus keys using union-find with chained merging.

    Two keys are in the same cluster if they're compatible (via
    keys_within_tolerance), allowing transitive merging. This handles
    both:
      - Continuous noisy distributions (e.g. codons 422-433 within tolerance 3
        with chained merging form one cluster), AND
      - Clean discrete clusters (e.g. codons 282 and 313 with gap > tolerance
        form two separate clusters).

    Returns list of clusters, where each cluster is a list of indices into keys.
    """
    n = len(keys)
    if n == 0:
        return []
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i in range(n):
        for j in range(i + 1, n):
            if keys_within_tolerance(keys[i], keys[j], codon_tolerance):
                union(i, j)

    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)
    return list(clusters.values())


def keys_within_tolerance(k1, k2, codon_tolerance):
    """Check whether two locus keys point to the same biological locus.

    - Both 'intron': same intron_number
    - Both 'cds': codon_number within tolerance (offset is not required to
      match, since small offset differences arise from annotation noise
      in the splice site detection)
    - Mixed 'cds' and 'intron': compatible if cds is in exon E and intron
      is intron E (immediately after exon E) or intron E-1 (immediately
      before exon E). Handles miniprot annotation differences where some
      samples annotated through the introner (treating it as exonic) and
      others split the gene around it.
    """
    if k1 is None or k2 is None:
        return False

    # Same location type
    if k1[0] == k2[0]:
        if k1[0] == 'intron':
            return k1[1] == k2[1]
        if k1[0] == 'cds':
            return abs(k1[1] - k2[1]) <= codon_tolerance

    # Mixed cds + intron — check exon/intron number compatibility
    cds_key = k1 if k1[0] == 'cds' else k2 if k2[0] == 'cds' else None
    intron_key = k1 if k1[0] == 'intron' else k2 if k2[0] == 'intron' else None

    if cds_key is not None and intron_key is not None:
        exon_num = cds_key[3]
        intron_num = intron_key[1]
        return intron_num == exon_num or intron_num == exon_num - 1

    return False


def classify_ortholog_group(members, codon_tolerance=CODON_TOLERANCE):
    """Classify the sharing status of an ortholog group.

    Uses location-aware comparison: introners in CDS exons are compared by
    codon position (with tolerance), introners in introns are compared by
    intron number (which is invariant to miniprot exon boundary differences).
    """
    # Only compare presence=1 members with high-confidence fingerprints
    typed = [m for m in members if m['confidence'] == 'high']

    if len(typed) < 2:
        return 'uncertain'

    # Check if there are members from both groups
    groups_present = {m['group'] for m in typed}
    has_both_groups = 'G1' in groups_present and 'G2' in groups_present

    if not has_both_groups:
        # Within-group comparison: cluster members and check if they all
        # collapse into a single cluster. Chained union-find clustering
        # correctly handles continuous distributions (annotation noise)
        # by transitively merging members within tolerance, while still
        # separating clear discrete clusters with gaps > tolerance.
        keys = [get_locus_key(m) for m in typed]
        keys = [k for k in keys if k is not None]
        if len(keys) < 2:
            return 'uncertain'

        clusters = cluster_keys(keys, codon_tolerance)
        if len(clusters) == 1:
            return 'consistent'
        return 'within_group_discordant'

    # Cross-group comparison
    g1_members = [m for m in typed if m['group'] == 'G1']
    g2_members = [m for m in typed if m['group'] == 'G2']

    if not g1_members or not g2_members:
        return 'uncertain'

    # Get representative locus key and family for each group
    def get_consensus(group_members):
        from collections import Counter
        keys = [get_locus_key(m) for m in group_members]
        keys = [k for k in keys if k is not None]
        if not keys:
            return None, None
        key_counts = Counter(keys)
        family_counts = Counter(m['family'] for m in group_members)
        return key_counts.most_common(1)[0][0], family_counts.most_common(1)[0][0]

    g1_key, g1_family = get_consensus(g1_members)
    g2_key, g2_family = get_consensus(g2_members)

    if g1_key is None or g2_key is None:
        return 'uncertain'

    same_family = (g1_family == g2_family)

    # Check if the two keys are compatible (same locus)
    compatible = keys_within_tolerance(g1_key, g2_key, codon_tolerance)

    if not compatible:
        return 'independent'

    # Compatible: determine how strong the match is
    # Exact match: same location type AND same position
    same_type = g1_key[0] == g2_key[0]
    exact_match = (g1_key == g2_key)

    if exact_match and same_family:
        return 'ancestral'

    if exact_match and not same_family:
        return 'same_site_diff_family'

    # Compatible but not exact (e.g., cds vs intron compatibility, or
    # cds with codon difference within tolerance)
    if same_family:
        return 'ambiguous'

    return 'ambiguous_diff_family'


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
                'location_type': fp.get('location_type', ''),
                'codon_number': fp.get('codon_number', ''),
                'codon_offset': fp.get('codon_offset', ''),
                'exon_number': fp.get('exon_number', ''),
                'intron_number': fp.get('intron_number', ''),
                'confidence': fp.get('confidence', 'no_fingerprint'),
                'gene_id': fp.get('gene_id', ''),
            })

        # Classify
        status = classify_ortholog_group(members, args.codon_tolerance)
        status_counts[status] += 1

        # Store summary for this group
        if members:
            typed = [m for m in members if m['confidence'] == 'high']

            def format_locus(m):
                if m['location_type'] == 'cds':
                    return f"{m['sample']}:cds.{m['codon_number']}.{m['codon_offset']}"
                if m['location_type'] == 'intron':
                    return f"{m['sample']}:intron.{m['intron_number']}"
                return f"{m['sample']}:?"

            group_summaries.append({
                'ortholog_id': oid,
                'sharing_status': status,
                'n_present': len(members),
                'n_fingerprinted': len(typed),
                'families': ';'.join(sorted(set(m['family'] for m in members))),
                'loci': ';'.join(format_locus(m) for m in typed),
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
                          'n_fingerprinted', 'families', 'loci']
        with open(args.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(group_summaries)
        print(f"Wrote {len(group_summaries)} group summaries to {args.summary}")


if __name__ == '__main__':
    main()

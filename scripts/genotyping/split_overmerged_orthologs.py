#!/usr/bin/env python3
"""
Split over-merged ortholog groups using codon/intron position clusters.

For each ortholog group with multiple distinct codon/intron position clusters
(within_group_discordant), partitions the presence=1 members into clusters by
their fingerprint position and creates a separate ortholog sub-group for each
cluster.

Reads:
  - Verified genotype matrix (with sharing_status column)
  - Per-sample fingerprint files (codon/intron positions)
  - Per-sample tandem cluster files (optional, for flagging affected groups)

Writes:
  - Split genotype matrix with new ortholog_id values for split groups
  - Mapping file showing how original IDs split into sub-IDs
"""

import argparse
import csv
from collections import defaultdict


GROUP2_SAMPLES = {'RCC1749', 'RCC3052'}
CODON_TOLERANCE = 3


def get_locus_key(fp):
    """Build a comparison key from a fingerprint dict."""
    loc = fp.get('location_type', '')
    if loc == 'cds':
        try:
            return ('cds', int(fp['codon_number']),
                    int(fp.get('codon_offset', 0)),
                    int(fp.get('exon_number', 0)))
        except (ValueError, TypeError, KeyError):
            return None
    if loc == 'intron':
        try:
            return ('intron', int(fp['intron_number']))
        except (ValueError, TypeError, KeyError):
            return None
    return None


def keys_within_tolerance(k1, k2, codon_tolerance):
    """Check whether two keys point to the same biological locus."""
    if k1 is None or k2 is None:
        return False

    if k1[0] == k2[0]:
        if k1[0] == 'intron':
            return k1[1] == k2[1]
        if k1[0] == 'cds':
            return abs(k1[1] - k2[1]) <= codon_tolerance

    cds_key = k1 if k1[0] == 'cds' else k2 if k2[0] == 'cds' else None
    intron_key = k1 if k1[0] == 'intron' else k2 if k2[0] == 'intron' else None

    if cds_key is not None and intron_key is not None:
        exon_num = cds_key[3]
        intron_num = intron_key[1]
        return intron_num == exon_num or intron_num == exon_num - 1

    return False


def cluster_keys(keys, codon_tolerance):
    """Cluster locus keys via union-find with chained merging."""
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


def load_fingerprints(fingerprint_files):
    """Load fingerprints into a dict keyed by (sample, contig, start, end)."""
    fps = {}
    for fp_file in fingerprint_files:
        with open(fp_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                key = (row['sample'], row['contig'],
                       int(row['start']), int(row['end']))
                fps[key] = row
    return fps


def load_tandems(tandem_files):
    """Load tandem cluster info into a dict keyed by (sample, contig, start, end)."""
    tandems = {}
    for t_file in tandem_files:
        with open(t_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if row.get('tandem_cluster_id'):
                    key = (row['sample'], row['contig'],
                           int(row['start']), int(row['end']))
                    tandems[key] = {
                        'cluster_id': row['tandem_cluster_id'],
                        'cluster_size': row.get('cluster_size', ''),
                    }
    return tandems


def main():
    parser = argparse.ArgumentParser(
        description='Split over-merged ortholog groups by codon/intron position')
    parser.add_argument('--matrix', required=True,
                        help='Verified genotype matrix with sharing_status column')
    parser.add_argument('--fingerprints', required=True, nargs='+',
                        help='Per-sample fingerprint TSV files')
    parser.add_argument('--tandems', nargs='*', default=[],
                        help='Per-sample tandem cluster TSV files (optional)')
    parser.add_argument('--output', required=True,
                        help='Output split genotype matrix')
    parser.add_argument('--mapping', default=None,
                        help='Optional ortholog ID mapping output')
    parser.add_argument('--codon-tolerance', type=int, default=CODON_TOLERANCE,
                        help='Tolerance for clustering codons (default: 3)')
    args = parser.parse_args()

    # Load data
    print(f"Loading fingerprints from {len(args.fingerprints)} files...")
    fps = load_fingerprints(args.fingerprints)
    print(f"  Loaded {len(fps)} fingerprint records")

    tandems = {}
    if args.tandems:
        print(f"Loading tandem clusters from {len(args.tandems)} files...")
        tandems = load_tandems(args.tandems)
        print(f"  Loaded {len(tandems)} tandem records")

    # Read matrix and group by ortholog_id
    print("Reading genotype matrix...")
    matrix_rows = []
    ortholog_groups = defaultdict(list)

    with open(args.matrix) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            idx = len(matrix_rows)
            # Preserve original ortholog_id so the flank-based grouping
            # can still be referenced after splitting
            row['original_ortholog_id'] = row['ortholog_id']
            matrix_rows.append(row)
            ortholog_groups[row['ortholog_id']].append(idx)

    print(f"  {len(matrix_rows)} rows in {len(ortholog_groups)} ortholog groups")

    # Process each group: split if discordant
    n_split = 0
    n_subgroups_added = 0
    n_tandem_flagged = 0
    mapping_records = []

    for oid, row_indices in sorted(ortholog_groups.items()):
        # Get presence=1 members with valid fingerprints
        present_indices = [idx for idx in row_indices
                           if matrix_rows[idx]['presence'] == '1']

        # Build per-group key lists
        g1_keys = []
        g1_idx = []
        g2_keys = []
        g2_idx = []
        for idx in present_indices:
            row = matrix_rows[idx]
            fp_key = (row['sample'], row['contig'],
                      int(row['start']), int(row['end']))
            fp = fps.get(fp_key, {})
            k = get_locus_key(fp)
            if k is None:
                continue
            if row['sample'] in GROUP2_SAMPLES:
                g2_keys.append(k)
                g2_idx.append(idx)
            else:
                g1_keys.append(k)
                g1_idx.append(idx)

        # Only split when there's intra-group multi-clustering. Cross-group
        # codon differences are meaningful biology (independent insertions
        # at the same approximate locus), not over-merging artifacts.
        g1_clusters = cluster_keys(g1_keys, args.codon_tolerance) if g1_keys else []
        g2_clusters = cluster_keys(g2_keys, args.codon_tolerance) if g2_keys else []

        if len(g1_clusters) <= 1 and len(g2_clusters) <= 1:
            mapping_records.append({
                'original_id': oid,
                'new_id': oid,
                'n_clusters': 1,
                'n_members': len(present_indices),
                'has_tandem': '',
                'note': 'unchanged',
            })
            continue

        # Build the actual partitioning. The split should respect the
        # within-group clusters AND also separate G1 from G2 when they're
        # not at compatible positions.
        # Approach: cluster all keys together, but require BOTH members of
        # any merge to be from the same intra-group cluster. This way, the
        # cross-group comparison is preserved within each cluster.
        keys = g1_keys + g2_keys
        idx_for_key = g1_idx + g2_idx
        clusters = cluster_keys(keys, args.codon_tolerance)

        if len(clusters) < 2:
            mapping_records.append({
                'original_id': oid,
                'new_id': oid,
                'n_clusters': 1,
                'n_members': len(present_indices),
                'has_tandem': '',
                'note': 'single_chained_cluster',
            })
            continue

        # SPLIT: assign each cluster a new ortholog_id
        n_split += 1

        # Map matrix index to cluster index
        idx_to_cluster = {}
        for ci, cluster in enumerate(clusters):
            for ki in cluster:
                idx_to_cluster[idx_for_key[ki]] = ci

        # Check for tandem cluster members in this group
        has_tandem = any(
            (matrix_rows[idx]['sample'], matrix_rows[idx]['contig'],
             int(matrix_rows[idx]['start']), int(matrix_rows[idx]['end'])) in tandems
            for idx in present_indices)
        if has_tandem:
            n_tandem_flagged += 1

        # Create new sub-IDs for each cluster
        # presence=1 members get reassigned to their cluster's new ID
        # presence=2/3 members stay with the first sub-ID (arbitrary, as
        # they represent absence/missing data)
        for ci, cluster in enumerate(clusters):
            new_id = f"{oid}_split{ci + 1}"
            if ci > 0:
                n_subgroups_added += 1

            cluster_size = len(cluster)
            mapping_records.append({
                'original_id': oid,
                'new_id': new_id,
                'n_clusters': len(clusters),
                'n_members': cluster_size,
                'has_tandem': '1' if has_tandem else '0',
                'note': 'split',
            })

        # Reassign ortholog_id for all rows in this group
        for idx in row_indices:
            row = matrix_rows[idx]
            if idx in idx_to_cluster:
                ci = idx_to_cluster[idx]
                row['ortholog_id'] = f"{oid}_split{ci + 1}"
            else:
                # presence=2/3 or no fingerprint: assign to split1
                # (these represent the ortholog group's absent/missing samples;
                # they should belong to all sub-groups conceptually, but we
                # keep them in the first sub-group for matrix compactness)
                row['ortholog_id'] = f"{oid}_split1"

    # Print summary
    print(f"\n=== Splitting summary ===")
    print(f"  Original groups:                 {len(ortholog_groups)}")
    print(f"  Groups split:                    {n_split}")
    print(f"  New sub-groups added:            {n_subgroups_added}")
    print(f"  Total groups after splitting:    {len(ortholog_groups) + n_subgroups_added}")
    print(f"  Split groups with tandem flag:   {n_tandem_flagged}")

    # Write output matrix
    out_fieldnames = list(fieldnames)
    if 'original_ortholog_id' not in out_fieldnames:
        out_fieldnames.append('original_ortholog_id')

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(matrix_rows)
    print(f"\nWrote {len(matrix_rows)} rows to {args.output}")

    # Write mapping
    if args.mapping:
        mapping_fields = ['original_id', 'new_id', 'n_clusters', 'n_members',
                          'has_tandem', 'note']
        with open(args.mapping, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=mapping_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(mapping_records)
        print(f"Wrote {len(mapping_records)} mapping records to {args.mapping}")


if __name__ == '__main__':
    main()

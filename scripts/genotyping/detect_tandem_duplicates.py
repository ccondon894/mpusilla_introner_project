#!/usr/bin/env python3
"""
Detect tandem duplications in introner candidate loci.

For each sample, identifies pairs of introners that are within a configurable
genomic window of each other and have very similar body sequences. Groups
them into tandem clusters via union-find.

Tandem clusters arise when an introner-containing region of a gene has been
duplicated, resulting in multiple introners with nearly identical sequences
(and identical flanking regions). The flank-based ortholog detection cannot
distinguish copies in a tandem cluster, leading to over-merging or
incorrect cross-pairing of introners across samples.

Output: TSV with one row per BED entry, annotated with:
  tandem_cluster_id  (integer or empty for singletons)
  cluster_size       (number of members in the cluster, or empty)
"""

import argparse
import csv
from collections import defaultdict

import pysam
from Bio.Align import PairwiseAligner


# Default parameters
DEFAULT_WINDOW = 10000     # max distance between tandem copies (bp)
DEFAULT_IDENTITY = 0.95    # minimum body sequence identity
DEFAULT_FLANK = 100        # flanking length to trim from BED coordinates
MIN_BODY_LENGTH = 30


def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def compute_identity(s1, s2, aligner):
    """Compute pairwise identity between two sequences."""
    if not s1 or not s2:
        return 0.0
    if abs(len(s1) - len(s2)) > max(len(s1), len(s2)) * 0.5:
        # Length difference too large to be a tandem copy
        return 0.0
    aln = aligner.align(s1.upper(), s2.upper())[0]
    counts = aln.counts()
    total = counts.identities + counts.mismatches + counts.gaps
    return counts.identities / total if total else 0.0


def extract_body(genome, contig, start, end, flank):
    """Extract introner body (excluding flanks) from indexed genome."""
    body_start = start + flank
    body_end = end - flank
    if body_end <= body_start:
        return None
    try:
        seq = genome.fetch(contig, body_start, body_end).upper()
    except (ValueError, KeyError):
        return None
    if len(seq) < MIN_BODY_LENGTH:
        return None
    return seq


def find_tandem_clusters(introners, genome, window_bp, min_identity,
                          aligner, flank):
    """Find tandem duplicate clusters within a sample's introners.

    Args:
        introners: list of dicts with contig, start, end
        genome: pysam.FastaFile
        window_bp: maximum genomic distance between tandem copies
        min_identity: minimum body identity to call as tandem
        aligner: PairwiseAligner instance
        flank: flanking length to exclude from body extraction

    Returns:
        dict {introner_index -> cluster_id} (only for cluster members)
    """
    # Extract body sequences
    bodies = {}
    for i, intr in enumerate(introners):
        body = extract_body(genome, intr['contig'], intr['start'],
                            intr['end'], flank)
        if body is not None:
            bodies[i] = body

    # Group by contig
    by_contig = defaultdict(list)
    for i in bodies:
        by_contig[introners[i]['contig']].append(i)

    # Sort each contig's introners by start position
    for contig in by_contig:
        by_contig[contig].sort(key=lambda i: introners[i]['start'])

    # Union-Find for clustering
    parent = {i: i for i in bodies}

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        # Path compression
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    # Compare nearby pairs within each contig
    n_compared = 0
    for contig, indices in by_contig.items():
        for i in range(len(indices)):
            idx1 = indices[i]
            start1 = introners[idx1]['start']
            end1 = introners[idx1]['end']

            for j in range(i + 1, len(indices)):
                idx2 = indices[j]
                start2 = introners[idx2]['start']
                # Check window: distance from end of idx1 to start of idx2
                if start2 - end1 > window_bp:
                    break

                identity = compute_identity(bodies[idx1], bodies[idx2], aligner)
                n_compared += 1
                if identity >= min_identity:
                    union(idx1, idx2)

    # Collect cluster members
    cluster_members = defaultdict(list)
    for i in bodies:
        cluster_members[find(i)].append(i)

    # Assign cluster IDs only to clusters with >= 2 members
    cluster_ids = {}
    next_id = 0
    for root, members in cluster_members.items():
        if len(members) >= 2:
            for m in members:
                cluster_ids[m] = next_id
            next_id += 1

    return cluster_ids, n_compared


def main():
    parser = argparse.ArgumentParser(
        description='Detect tandem duplications in introner BED files')
    parser.add_argument('--bed', required=True,
                        help='Candidate loci BED file')
    parser.add_argument('--genome', required=True,
                        help='Indexed genome FASTA')
    parser.add_argument('--sample', required=True,
                        help='Sample name')
    parser.add_argument('--output', required=True,
                        help='Output tandem cluster TSV')
    parser.add_argument('--window-bp', type=int, default=DEFAULT_WINDOW,
                        help=f'Max distance between tandem copies '
                             f'(default: {DEFAULT_WINDOW})')
    parser.add_argument('--min-identity', type=float, default=DEFAULT_IDENTITY,
                        help=f'Minimum body identity for tandem '
                             f'(default: {DEFAULT_IDENTITY})')
    parser.add_argument('--flanking-length', type=int, default=DEFAULT_FLANK,
                        help=f'Flanking length to trim (default: {DEFAULT_FLANK})')
    args = parser.parse_args()

    # Read BED
    introners = []
    with open(args.bed) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            introners.append({
                'contig': row['#chrom'],
                'start': int(row['start']),
                'end': int(row['end']),
                'introner_id': row['name'],
                'family': row.get('family', ''),
                'gene_info': row.get('gene_info', ''),
            })

    print(f"Read {len(introners)} introners from {args.bed}")

    # Open genome and aligner
    genome = pysam.FastaFile(args.genome)
    aligner = make_aligner()

    # Find clusters
    cluster_ids, n_compared = find_tandem_clusters(
        introners, genome, args.window_bp, args.min_identity,
        aligner, args.flanking_length)

    genome.close()

    # Compute cluster sizes
    cluster_sizes = defaultdict(int)
    for cid in cluster_ids.values():
        cluster_sizes[cid] += 1

    n_in_clusters = len(cluster_ids)
    n_clusters = len(cluster_sizes)
    print(f"Compared {n_compared} pairs")
    print(f"Found {n_clusters} tandem clusters with {n_in_clusters} member introners")
    if cluster_sizes:
        sizes = sorted(cluster_sizes.values())
        print(f"  Cluster sizes: min={min(sizes)}, "
              f"median={sizes[len(sizes)//2]}, max={max(sizes)}")

    # Write output
    fieldnames = ['sample', 'introner_id', 'contig', 'start', 'end',
                  'family', 'gene_info', 'tandem_cluster_id', 'cluster_size']

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for i, intr in enumerate(introners):
            cid = cluster_ids.get(i, '')
            writer.writerow({
                'sample': args.sample,
                'introner_id': intr['introner_id'],
                'contig': intr['contig'],
                'start': intr['start'],
                'end': intr['end'],
                'family': intr['family'],
                'gene_info': intr['gene_info'],
                'tandem_cluster_id': cid,
                'cluster_size': cluster_sizes[cid] if cid != '' else '',
            })

    print(f"Wrote {len(introners)} records to {args.output}")


if __name__ == '__main__':
    main()

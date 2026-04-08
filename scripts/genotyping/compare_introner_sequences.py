#!/usr/bin/env python3
"""
Compare introner body sequences within ortholog groups to refine sharing status.

For each ortholog group with ≥2 presence=1 members, extracts introner body
sequences from indexed genome FASTAs, computes pairwise sequence identity,
and combines with the codon-level sharing_status to produce a refined
classification.

All codon-level results are preserved. This script adds:
  - cross_group_identity:    median G1-vs-G2 pairwise identity
  - within_group_identity:   median within-group pairwise identity
  - refined_sharing_status:  final verdict combining both evidence types
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from statistics import median

import pysam
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


GROUP2_SAMPLES = {'RCC1749', 'RCC3052'}
G1_PREFERRED = 'CCMP1545'
G2_PREFERRED = 'RCC1749'

MIN_BODY_LENGTH = 30
MAX_N_FRACTION = 0.20


def load_genomes(genome_dir, samples):
    """Open pysam FastaFile handles for all samples."""
    genomes = {}
    for sample in samples:
        fa_path = os.path.join(genome_dir, f"{sample}.vg_paths.fa")
        if os.path.exists(fa_path):
            genomes[sample] = pysam.FastaFile(fa_path)
        else:
            print(f"  Warning: genome not found for {sample}: {fa_path}")
    return genomes


def extract_body(genomes, sample, contig, start, end, flanking):
    """Extract introner body sequence (excluding flanking regions)."""
    if sample not in genomes:
        return None

    body_start = start + flanking
    body_end = end - flanking
    if body_end <= body_start:
        return None

    try:
        seq = genomes[sample].fetch(contig, body_start, body_end).upper()
    except (ValueError, KeyError):
        return None

    if len(seq) < MIN_BODY_LENGTH:
        return None

    n_count = seq.count('N')
    if n_count / len(seq) > MAX_N_FRACTION:
        return None

    return seq


def make_aligner():
    """Create a PairwiseAligner with standard introner comparison parameters."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def compute_identity_from_alignment(alignment):
    """Compute sequence identity from a PairwiseAligner alignment.

    Identity = matching non-gap positions / total non-gap positions.
    """
    aligned = alignment.format().split('\n')
    # PairwiseAligner.format() gives: seq1\n annotation\n seq2\n
    # But the format may vary; use the aligned sequences directly
    seq1_aligned = aligned[0]
    seq2_aligned = aligned[2]

    matches = 0
    non_gap = 0
    for a, b in zip(seq1_aligned, seq2_aligned):
        if a != '-' or b != '-':
            non_gap += 1
            if a == b and a != '-':
                matches += 1

    if non_gap == 0:
        return 0.0
    return matches / non_gap


def pairwise_identity(seq1, seq2, aligner):
    """Compute pairwise identity, testing both forward and reverse complement.

    Returns (identity, orientation) where identity is a float [0, 1].
    """
    seq1_upper = seq1.upper()
    seq2_upper = seq2.upper()

    # Forward alignment
    alignments_fwd = aligner.align(seq1_upper, seq2_upper)
    fwd_identity = compute_identity_from_alignment(alignments_fwd[0])

    # Reverse complement alignment
    seq1_rc = str(Seq(seq1_upper).reverse_complement())
    alignments_rc = aligner.align(seq1_rc, seq2_upper)
    rc_identity = compute_identity_from_alignment(alignments_rc[0])

    if fwd_identity >= rc_identity:
        return fwd_identity, 'forward'
    return rc_identity, 'reverse'


def pick_representative(members, preferred):
    """Pick a representative member, preferring a specific sample."""
    for m in members:
        if m['sample'] == preferred:
            return m
    return members[0]


def compute_group_identities(members, genomes, aligner, flanking):
    """Compute within-group and cross-group identity for an ortholog group.

    Args:
        members: list of dicts with sample, contig, start, end, group keys
        genomes: dict of pysam.FastaFile handles
        aligner: PairwiseAligner instance
        flanking: flanking length to trim

    Returns:
        (within_group_identity, cross_group_identity, n_extracted)
        where identities are floats or None if not computable.
    """
    # Extract sequences for all members
    seqs = []
    for m in members:
        seq = extract_body(genomes, m['sample'], m['contig'],
                           m['start'], m['end'], flanking)
        if seq:
            seqs.append({'sample': m['sample'], 'group': m['group'], 'seq': seq})

    if len(seqs) < 2:
        return None, None, len(seqs)

    # Split by group
    g1_seqs = [s for s in seqs if s['group'] == 'G1']
    g2_seqs = [s for s in seqs if s['group'] == 'G2']

    # Cross-group identity
    cross_identity = None
    if g1_seqs and g2_seqs:
        cross_identities = []
        g1_rep = pick_representative(g1_seqs, G1_PREFERRED)
        g2_rep = pick_representative(g2_seqs, G2_PREFERRED)

        # Always compute representative pair
        ident, _ = pairwise_identity(g1_rep['seq'], g2_rep['seq'], aligner)
        cross_identities.append(ident)

        # For small groups, compute all G1-vs-G2 pairs for a robust median
        if len(g1_seqs) <= 5 and len(g2_seqs) <= 2:
            for s1 in g1_seqs:
                for s2 in g2_seqs:
                    if s1['sample'] == g1_rep['sample'] and s2['sample'] == g2_rep['sample']:
                        continue  # already computed
                    ident, _ = pairwise_identity(s1['seq'], s2['seq'], aligner)
                    cross_identities.append(ident)

        cross_identity = median(cross_identities)

    # Within-group identity (representative pair from largest group)
    within_identity = None
    largest_group = g1_seqs if len(g1_seqs) >= len(g2_seqs) else g2_seqs
    if len(largest_group) >= 2:
        rep = largest_group[0]
        other = largest_group[1]
        ident, _ = pairwise_identity(rep['seq'], other['seq'], aligner)
        within_identity = ident

    return within_identity, cross_identity, len(seqs)


def refine_status(sharing_status, within_identity, cross_identity,
                  within_threshold, cross_threshold):
    """Combine codon-level sharing_status with sequence identity evidence.

    Uses two distinct thresholds:
      - within_threshold: applied to within-group-only ortholog groups, where
        high identity is expected (~92% median in practice)
      - cross_threshold: applied to cross-group ortholog groups, where
        substantial divergence between groups is expected

    Returns refined_sharing_status string.
    """
    # Within-group-only ortholog groups: use within-group identity check
    if sharing_status == 'consistent':
        if within_identity is None:
            return 'consistent'
        if within_identity >= within_threshold:
            return 'consistent_confirmed'
        return 'consistent_low_identity'

    if sharing_status == 'within_group_discordant':
        return 'within_group_discordant'

    # Cross-group ortholog groups: use cross-group identity check
    if cross_identity is None:
        return sharing_status

    high_cross = cross_identity >= cross_threshold

    mapping = {
        ('ancestral', True): 'ancestral_confirmed',
        ('ancestral', False): 'convergent_same_site',
        ('same_site_diff_family', True): 'ancestral_diff_family',
        ('same_site_diff_family', False): 'convergent_same_site',
        ('ambiguous', True): 'ancestral_resolved',
        ('ambiguous', False): 'independent_resolved',
        ('ambiguous_diff_family', True): 'ambiguous',
        ('ambiguous_diff_family', False): 'independent_resolved',
        ('independent', True): 'independent',
        ('independent', False): 'independent',
        ('uncertain', True): 'likely_ancestral',
        ('uncertain', False): 'likely_independent',
    }

    return mapping.get((sharing_status, high_cross), sharing_status)


def main():
    parser = argparse.ArgumentParser(
        description='Compare introner sequences within ortholog groups')
    parser.add_argument('--matrix', required=True,
                        help='Verified genotype matrix with sharing_status column')
    parser.add_argument('--genome-dir', required=True,
                        help='Directory with indexed genome FASTAs')
    parser.add_argument('--output', required=True,
                        help='Output matrix with refined sharing status')
    parser.add_argument('--summary', default=None,
                        help='Optional per-group summary TSV')
    parser.add_argument('--flanking-length', type=int, default=100,
                        help='Flanking length to trim from coordinates (default: 100)')
    parser.add_argument('--within-group-identity-threshold', type=float, default=0.80,
                        help='Identity threshold for within-group comparisons (default: 0.80)')
    parser.add_argument('--cross-group-identity-threshold', type=float, default=0.60,
                        help='Identity threshold for cross-group comparisons (default: 0.60)')
    parser.add_argument('--coverage-threshold', type=float, default=0.80,
                        help='Minimum alignment coverage (default: 0.80)')
    args = parser.parse_args()

    # Read matrix and group by ortholog_id
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

    # Collect all samples and load genomes
    all_samples = sorted({row['sample'] for row in matrix_rows})
    print(f"Loading genomes for {len(all_samples)} samples...")
    genomes = load_genomes(args.genome_dir, all_samples)
    print(f"  Loaded {len(genomes)} genomes")

    # Set up aligner
    aligner = make_aligner()

    # Process each ortholog group
    status_counts = defaultdict(int)
    group_summaries = []
    groups_processed = 0

    for oid, row_indices in ortholog_groups.items():
        # Collect presence=1 members
        members = []
        sharing_status = matrix_rows[row_indices[0]].get('sharing_status', '')

        for idx in row_indices:
            row = matrix_rows[idx]
            if row['presence'] != '1':
                continue
            sample = row['sample']
            members.append({
                'sample': sample,
                'contig': row['contig'],
                'start': int(row['start']),
                'end': int(row['end']),
                'group': 'G2' if sample in GROUP2_SAMPLES else 'G1',
                'family': row.get('family', ''),
            })

        # Compute identities
        within_ident = None
        cross_ident = None
        n_extracted = 0

        if len(members) >= 2:
            within_ident, cross_ident, n_extracted = compute_group_identities(
                members, genomes, aligner, args.flanking_length)

        # Refine status
        refined = refine_status(
            sharing_status, within_ident, cross_ident,
            args.within_group_identity_threshold,
            args.cross_group_identity_threshold)
        status_counts[refined] += 1

        # Format identity values for output
        cross_str = f"{cross_ident:.4f}" if cross_ident is not None else ''
        within_str = f"{within_ident:.4f}" if within_ident is not None else ''

        # Apply to all rows in the group
        for idx in row_indices:
            matrix_rows[idx]['cross_group_identity'] = cross_str
            matrix_rows[idx]['within_group_identity'] = within_str
            matrix_rows[idx]['refined_sharing_status'] = refined

        # Summary record
        families = sorted(set(m['family'] for m in members if m['family']))
        g1_rep = [m['sample'] for m in members if m['group'] == 'G1'][:1]
        g2_rep = [m['sample'] for m in members if m['group'] == 'G2'][:1]

        group_summaries.append({
            'ortholog_id': oid,
            'sharing_status': sharing_status,
            'refined_sharing_status': refined,
            'n_present': len(members),
            'n_extracted': n_extracted,
            'cross_group_identity': cross_str,
            'within_group_identity': within_str,
            'families': ';'.join(families),
            'g1_representative': g1_rep[0] if g1_rep else '',
            'g2_representative': g2_rep[0] if g2_rep else '',
        })

        groups_processed += 1
        if groups_processed % 1000 == 0:
            print(f"  Processed {groups_processed}/{len(ortholog_groups)} groups...")

    # Close genomes
    for fa in genomes.values():
        fa.close()

    # Print summary
    print(f"\nRefined sharing status (within_threshold={args.within_group_identity_threshold}, "
          f"cross_threshold={args.cross_group_identity_threshold}):")
    for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")
    print(f"  {'TOTAL':30s}  {len(ortholog_groups):5d}")

    # Write output matrix
    out_fieldnames = fieldnames + ['cross_group_identity', 'within_group_identity',
                                    'refined_sharing_status']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(matrix_rows)

    print(f"\nWrote {len(matrix_rows)} rows to {args.output}")

    # Write optional summary
    if args.summary:
        summary_fields = ['ortholog_id', 'sharing_status', 'refined_sharing_status',
                          'n_present', 'n_extracted', 'cross_group_identity',
                          'within_group_identity', 'families',
                          'g1_representative', 'g2_representative']
        with open(args.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(group_summaries)
        print(f"Wrote {len(group_summaries)} group summaries to {args.summary}")


if __name__ == '__main__':
    main()

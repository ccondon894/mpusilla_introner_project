#!/usr/bin/env python3
"""
Compare introner body sequences within ortholog groups to refine classification.

For each ortholog group with >=2 presence=1 members, extracts introner body
sequences from indexed genome FASTAs, computes pairwise sequence identity,
and refines the codon-level within_group_status and cross_group_status.

Key refinements:
  - Resolves 'uncertain' within-group groups (mostly intergenic) to
    'consistent' when body identity >= threshold AND same family
  - Resolves 'uncertain' cross-group cases to 'likely_ancestral' or
    'likely_independent' based on body identity
  - Flags 'consistent' groups with low within-group identity
  - Flags 'ancestral' groups with low cross-group identity

Overwrites within_group_status and cross_group_status with refined values.
Adds within_group_identity and cross_group_identity as audit columns.
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
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

    Identity = identities / (identities + mismatches + gaps).
    """
    counts = alignment.counts()
    total = counts.identities + counts.mismatches + counts.gaps
    if total == 0:
        return 0.0
    return counts.identities / total


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

    Returns (within_group_identity, cross_group_identity, n_extracted)
    where identities are floats or None if not computable.
    """
    # Extract sequences for all members
    seqs = []
    for m in members:
        seq = extract_body(genomes, m['sample'], m['contig'],
                           m['start'], m['end'], flanking)
        if seq:
            seqs.append({'sample': m['sample'], 'group': m['group'],
                         'seq': seq, 'family': m.get('family', '')})

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

        ident, _ = pairwise_identity(g1_rep['seq'], g2_rep['seq'], aligner)
        cross_identities.append(ident)

        # For small groups, compute all G1-vs-G2 pairs for a robust median
        if len(g1_seqs) <= 5 and len(g2_seqs) <= 2:
            for s1 in g1_seqs:
                for s2 in g2_seqs:
                    if (s1['sample'] == g1_rep['sample'] and
                            s2['sample'] == g2_rep['sample']):
                        continue
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


def refine_within_group(within_status, within_identity, within_threshold,
                        families_consistent):
    """Refine within_group_status using sequence identity.

    Args:
        within_status: codon-level within_group_status
        within_identity: float or None
        within_threshold: identity threshold (e.g. 0.80)
        families_consistent: bool, True if all members share the same family
    """
    if within_status == 'singleton':
        return 'singleton'

    if within_status == 'discordant':
        return 'discordant'

    if within_status == 'consistent':
        if within_identity is None:
            return 'consistent'
        if within_identity >= within_threshold:
            return 'consistent'
        return 'low_identity'

    if within_status == 'uncertain':
        # Resolve using sequence identity + family check
        if within_identity is None:
            return 'uncertain'
        if within_identity >= within_threshold and families_consistent:
            return 'consistent'
        if within_identity < within_threshold:
            return 'low_identity'
        # High identity but different families — keep uncertain
        return 'uncertain'

    return within_status


def refine_cross_group(cross_status, cross_identity, cross_threshold):
    """Refine cross_group_status using sequence identity.

    Args:
        cross_status: codon-level cross_group_status
        cross_identity: float or None
        cross_threshold: identity threshold (e.g. 0.60)
    """
    if cross_status == 'NA':
        return 'NA'

    if cross_status == 'ancestral':
        if cross_identity is not None and cross_identity < cross_threshold:
            return 'ancestral_low_identity'
        return 'ancestral'

    if cross_status == 'independent':
        return 'independent'

    if cross_status == 'uncertain':
        if cross_identity is None:
            return 'uncertain'
        if cross_identity >= cross_threshold:
            return 'likely_ancestral'
        return 'likely_independent'

    return cross_status


def check_family_consistency(members):
    """Check if all presence=1 members within each clade share the same family.

    Returns True if families are consistent within each clade (ignoring
    empty/missing families). Returns True for single-member clades.
    """
    by_clade = defaultdict(list)
    for m in members:
        fam = m.get('family', '')
        if fam:
            by_clade[m['group']].append(fam)

    for clade, families in by_clade.items():
        if len(families) < 2:
            continue
        if len(set(families)) > 1:
            return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description='Compare introner sequences within ortholog groups')
    parser.add_argument('--matrix', required=True,
                        help='Verified genotype matrix with within_group_status '
                             'and cross_group_status columns')
    parser.add_argument('--genome-dir', required=True,
                        help='Directory with indexed genome FASTAs')
    parser.add_argument('--output', required=True,
                        help='Output matrix with refined status columns')
    parser.add_argument('--summary', default=None,
                        help='Optional per-group summary TSV')
    parser.add_argument('--flanking-length', type=int, default=100,
                        help='Flanking length to trim from coordinates (default: 100)')
    parser.add_argument('--within-group-identity-threshold', type=float, default=0.80,
                        help='Identity threshold for within-group comparisons (default: 0.80)')
    parser.add_argument('--cross-group-identity-threshold', type=float, default=0.60,
                        help='Identity threshold for cross-group comparisons (default: 0.60)')
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
    within_counts = defaultdict(int)
    cross_counts = defaultdict(int)
    within_transitions = defaultdict(int)
    cross_transitions = defaultdict(int)
    group_summaries = []
    groups_processed = 0

    for oid, row_indices in ortholog_groups.items():
        # Read codon-level statuses
        first_row = matrix_rows[row_indices[0]]
        within_status = first_row.get('within_group_status', '')
        cross_status = first_row.get('cross_group_status', '')

        # Collect presence=1 members
        members = []
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

        # Check family consistency for within-group resolution
        families_consistent = check_family_consistency(members)

        # Refine both statuses
        refined_within = refine_within_group(
            within_status, within_ident,
            args.within_group_identity_threshold, families_consistent)
        refined_cross = refine_cross_group(
            cross_status, cross_ident,
            args.cross_group_identity_threshold)

        within_counts[refined_within] += 1
        cross_counts[refined_cross] += 1

        # Track transitions
        if within_status != refined_within:
            within_transitions[f"{within_status} -> {refined_within}"] += 1
        if cross_status != refined_cross:
            cross_transitions[f"{cross_status} -> {refined_cross}"] += 1

        # Format identity values for output
        cross_str = f"{cross_ident:.4f}" if cross_ident is not None else ''
        within_str = f"{within_ident:.4f}" if within_ident is not None else ''

        # Apply refined statuses to all rows in the group
        for idx in row_indices:
            matrix_rows[idx]['within_group_status'] = refined_within
            matrix_rows[idx]['cross_group_status'] = refined_cross
            matrix_rows[idx]['within_group_identity'] = within_str
            matrix_rows[idx]['cross_group_identity'] = cross_str

        # Summary record
        families = sorted(set(m['family'] for m in members if m['family']))
        g1_rep = [m['sample'] for m in members if m['group'] == 'G1'][:1]
        g2_rep = [m['sample'] for m in members if m['group'] == 'G2'][:1]

        group_summaries.append({
            'ortholog_id': oid,
            'within_group_status': refined_within,
            'cross_group_status': refined_cross,
            'n_present': len(members),
            'n_extracted': n_extracted,
            'cross_group_identity': cross_str,
            'within_group_identity': within_str,
            'families': ';'.join(families),
            'families_consistent': str(families_consistent),
        })

        groups_processed += 1
        if groups_processed % 1000 == 0:
            print(f"  Processed {groups_processed}/{len(ortholog_groups)} groups...")

    # Close genomes
    for fa in genomes.values():
        fa.close()

    # Print summary
    print(f"\nWithin-group status (refined, threshold={args.within_group_identity_threshold}):")
    for status, count in sorted(within_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")

    print(f"\nCross-group status (refined, threshold={args.cross_group_identity_threshold}):")
    for status, count in sorted(cross_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")

    print(f"  {'TOTAL':30s}  {len(ortholog_groups):5d}")

    if within_transitions:
        print(f"\nWithin-group transitions:")
        for t, count in sorted(within_transitions.items(), key=lambda x: -x[1]):
            print(f"  {t:45s}  {count:5d}")

    if cross_transitions:
        print(f"\nCross-group transitions:")
        for t, count in sorted(cross_transitions.items(), key=lambda x: -x[1]):
            print(f"  {t:45s}  {count:5d}")

    # Write output matrix
    # Remove old columns if re-running, add new identity columns
    out_fieldnames = [f for f in fieldnames
                      if f not in ('cross_group_identity', 'within_group_identity',
                                   'refined_sharing_status', 'sharing_status')]
    # Ensure within_group_status and cross_group_status are present
    if 'within_group_status' not in out_fieldnames:
        out_fieldnames.append('within_group_status')
    if 'cross_group_status' not in out_fieldnames:
        out_fieldnames.append('cross_group_status')
    out_fieldnames += ['within_group_identity', 'cross_group_identity']

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(matrix_rows)

    print(f"\nWrote {len(matrix_rows)} rows to {args.output}")

    # Write optional summary
    if args.summary:
        summary_fields = ['ortholog_id', 'within_group_status',
                          'cross_group_status', 'n_present', 'n_extracted',
                          'cross_group_identity', 'within_group_identity',
                          'families', 'families_consistent']
        with open(args.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(group_summaries)
        print(f"Wrote {len(group_summaries)} group summaries to {args.summary}")


if __name__ == '__main__':
    main()

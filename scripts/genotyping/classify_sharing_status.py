#!/usr/bin/env python3
"""
Classify sharing status for introner ortholog groups.

For each ortholog group in the genotype matrix, compares codon-level insertion
fingerprints among presence=1 members to produce two independent classifications:

  within_group_status:  consistency of insertion site within each clade
    - consistent:   all fingerprinted members within each clade cluster together
    - discordant:   at least one clade has members at multiple distinct sites
    - uncertain:    insufficient fingerprint data (<2 typed members total)
    - singleton:    only 1 presence=1 member in the entire group

  cross_group_status:   ancestry relationship between G1 and G2
    - ancestral:    same insertion site (exact or aa-context match) + same family
    - independent:  different insertion site, different families, or only
                    legacy codon match without conserved protein context
    - uncertain:    insufficient fingerprint data for cross-group comparison
    - NA:           no cross-group members (within-group-only group)

Reads:
  - Oriented genotype matrix (TSV)
  - Per-sample insertion fingerprint files (TSV)

Outputs:
  - Verified genotype matrix with within_group_status and cross_group_status columns
  - Summary statistics (printed to stdout)
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict


# Codon tolerance for cross-group comparisons
# Small offsets (1-3 codons) may arise from protein annotation differences
# between divergent genomes, not from different insertion sites
CODON_TOLERANCE = 3

# The flanking aa context is extracted as 10 amino acids (5 on each side
# of the insertion). For comparison, we slide a 6-mer window across both
# contexts looking for a match with ≤1 mismatch. The wider extracted
# context + sliding window lets us tolerate up to ±2 codon shifts in
# the computed insertion position (which arise from splice site detection
# variability across samples) while being specific enough to distinguish
# genuinely different insertion sites.
AA_CONTEXT_WINDOW_LEN = 6
AA_CONTEXT_MISMATCH_TOLERANCE = 1

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
    """Build a hybrid comparison key for an introner.

    Returns a ('hybrid', aa_ctx_or_none, legacy_key_or_none) tuple containing
    both the amino acid flanking context AND the legacy codon/intron key.
    Amino-acid context is authoritative for mixed CDS/intron comparisons when
    both contexts are available. The legacy exon/intron-number fallback is
    retained when either context is unavailable.

    Rationale: the aa context is the primary signal (annotation-invariant),
    but it breaks when miniprot produces inconsistent CDS lengths across
    samples (frameshifts cause totally different downstream amino acids).
    The legacy codon/intron fallback catches these cases.

    Returns None if the member has no usable location data at all.
    """
    aa_ctx = member.get('flanking_aa_context', '')
    if not aa_ctx:
        aa_ctx = None

    legacy_key = None
    loc = member.get('location_type', '')
    if loc == 'cds':
        try:
            legacy_key = ('cds', int(member['codon_number']),
                          int(member['codon_offset']),
                          int(member['exon_number']))
        except (ValueError, TypeError):
            pass
    elif loc == 'intron':
        try:
            legacy_key = ('intron', int(member['intron_number']))
        except (ValueError, TypeError):
            pass

    if aa_ctx is None and legacy_key is None:
        return None
    return ('hybrid', aa_ctx, legacy_key)


def aa_contexts_match(ctx1, ctx2,
                       window_len=AA_CONTEXT_WINDOW_LEN,
                       max_mismatches=AA_CONTEXT_MISMATCH_TOLERANCE):
    """Check if two aa context strings match via sliding-window comparison.

    Slides a window of `window_len` amino acids across both contexts and
    returns True if any pair of windows matches with ≤max_mismatches
    position-wise mismatches. This tolerates small shifts in the computed
    insertion position (since splice site detection is sometimes off by
    ±1-2 codons) while remaining specific: a 6-mer exact match has
    ~1 in 20^6 chance by random, so the sliding window adds at most
    a factor of 25 more comparisons per pair.

    Padding characters ('-') are not treated specially — they count as
    mismatches if they don't align, which is the desired behavior at
    gene edges.
    """
    if not ctx1 or not ctx2:
        return False
    if len(ctx1) < window_len or len(ctx2) < window_len:
        # Fall back to position-wise comparison on the shorter length
        n = min(len(ctx1), len(ctx2))
        mm = sum(1 for i in range(n) if ctx1[i] != ctx2[i])
        return mm <= max_mismatches

    for i in range(len(ctx1) - window_len + 1):
        w1 = ctx1[i:i + window_len]
        # Skip windows that are mostly padding
        if w1.count('-') > max_mismatches:
            continue
        for j in range(len(ctx2) - window_len + 1):
            w2 = ctx2[j:j + window_len]
            if w2.count('-') > max_mismatches:
                continue
            mm = sum(1 for a, b in zip(w1, w2) if a != b)
            if mm <= max_mismatches:
                return True
    return False


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


def _legacy_keys_match(lk1, lk2, codon_tolerance):
    """Check two legacy codon/intron keys for compatibility."""
    if lk1 is None or lk2 is None:
        return False

    if lk1[0] == lk2[0]:
        if lk1[0] == 'intron':
            return lk1[1] == lk2[1]
        if lk1[0] == 'cds':
            return abs(lk1[1] - lk2[1]) <= codon_tolerance

    # Mixed cds + intron
    cds_key = lk1 if lk1[0] == 'cds' else lk2 if lk2[0] == 'cds' else None
    intron_key = lk1 if lk1[0] == 'intron' else lk2 if lk2[0] == 'intron' else None
    if cds_key is not None and intron_key is not None:
        exon_num = cds_key[3]
        intron_num = intron_key[1]
        return intron_num == exon_num or intron_num == exon_num - 1

    return False


def keys_within_tolerance(k1, k2, codon_tolerance,
                            aa_mismatch_tolerance=AA_CONTEXT_MISMATCH_TOLERANCE):
    """Check whether two hybrid locus keys point to the same biological locus.

    Each key is a ('hybrid', aa_ctx, legacy_key) tuple. Two keys match if:
      - Both have aa contexts and they match via sliding-window comparison, OR
      - Both have legacy codon/intron keys and they're compatible, provided a
        mixed CDS/intron comparison does not have contradictory AA contexts.

    The hybrid logic makes the comparison robust to miniprot annotation
    inconsistencies: when CDS lengths differ between samples (causing
    frameshifts that break aa context comparison), the legacy fallback
    can still match. When the CDS is consistent and the only difference
    is codon numbering (annotation noise), the aa context matches.
    """
    if k1 is None or k2 is None:
        return False

    # Handle legacy keys directly (shouldn't happen with new code, but safe)
    if k1[0] != 'hybrid' or k2[0] != 'hybrid':
        return _legacy_keys_match(k1, k2, codon_tolerance)

    _, aa1, legacy1 = k1
    _, aa2, legacy2 = k2

    # For mixed CDS/intron annotations, explicit AA disagreement is evidence
    # against equivalence and may not be overridden by the permissive exon N
    # <-> intron N/N-1 fallback. Preserve that fallback when AA evidence is
    # genuinely unavailable.
    if aa1 and aa2:
        if aa_contexts_match(aa1, aa2, max_mismatches=aa_mismatch_tolerance):
            return True
        if legacy1 is not None and legacy2 is not None:
            if {legacy1[0], legacy2[0]} == {'cds', 'intron'}:
                return False

    # Fall back to legacy key match
    return _legacy_keys_match(legacy1, legacy2, codon_tolerance)


def _assess_within_group(typed_members, codon_tolerance):
    """Assess within-group consistency of fingerprinted members.

    Checks each clade independently. If any clade has >=2 typed members
    that cluster to multiple distinct sites, returns 'discordant'.
    Otherwise returns 'consistent'.

    Returns 'uncertain' if fewer than 2 typed members total.
    """
    if len(typed_members) < 2:
        return 'uncertain'

    # Check each clade independently
    by_clade = defaultdict(list)
    for m in typed_members:
        by_clade[m['group']].append(m)

    for clade, clade_members in by_clade.items():
        if len(clade_members) < 2:
            continue  # trivially consistent with 0-1 members
        keys = [get_locus_key(m) for m in clade_members]
        keys = [k for k in keys if k is not None]
        if len(keys) < 2:
            continue
        clusters = cluster_keys(keys, codon_tolerance)
        if len(clusters) > 1:
            return 'discordant'

    return 'consistent'


def _assess_cross_group(typed_members, codon_tolerance):
    """Assess cross-group relationship between G1 and G2.

    Returns a (status, reason) tuple:
      status: one of 'ancestral', 'independent', 'uncertain'
      reason: fine-grained label for which classification branch fired
              ('ancestral_exact', 'ancestral_aa_match',
               'different_position', 'close_diff_aa_ctx',
               'compatible_diff_family', 'exact_diff_family',
               'uncertain_no_keys', 'uncertain_single_clade')

    Different-family cases are classified as independent: family assignment
    is strong evidence that the introners are separate insertion events
    regardless of positional similarity. The reason column distinguishes
    'different_position' and 'close_diff_aa_ctx' from the genuinely-
    same-site cases; downstream splitting can use this to separate
    mispaired clade-specific insertions from real convergent same-site
    insertions.
    """
    g1_members = [m for m in typed_members if m['group'] == 'G1']
    g2_members = [m for m in typed_members if m['group'] == 'G2']

    if not g1_members or not g2_members:
        return 'uncertain', 'uncertain_single_clade'

    # Get representative locus key and family for each group
    def get_consensus(group_members):
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
        return 'uncertain', 'uncertain_no_keys'

    same_family = (g1_family == g2_family)

    # Check if the two keys are compatible (same locus)
    compatible = keys_within_tolerance(g1_key, g2_key, codon_tolerance)

    if not compatible:
        # Neither aa context nor legacy keys match → truly different sites
        return 'independent', 'different_position'

    # Compatible: determine how strong the match is.
    # Exact match: aa contexts are identical, OR (if no aa context) the
    # legacy keys are exactly equal.
    exact_match = False
    if g1_key[0] == 'hybrid' and g2_key[0] == 'hybrid':
        _, aa1, legacy1 = g1_key
        _, aa2, legacy2 = g2_key
        if aa1 and aa2 and aa1 == aa2:
            exact_match = True
        elif aa1 is None and aa2 is None and legacy1 == legacy2:
            exact_match = True
    else:
        exact_match = (g1_key == g2_key)

    if exact_match and same_family:
        return 'ancestral', 'ancestral_exact'

    if exact_match and not same_family:
        return 'independent', 'exact_diff_family'

    if not same_family:
        # Compatible positions but different families (not exact match)
        return 'independent', 'compatible_diff_family'

    # Compatible but not exact, same family. Distinguish by HOW the keys
    # matched: if the amino acid contexts around the insertion site are
    # conserved (sliding window match), the insertion site is genuinely
    # the same and this is ancestral. If only the legacy codon/intron
    # numbers matched, the protein context is completely different and
    # the position match is likely annotation noise → independent.
    if g1_key[0] == 'hybrid' and g2_key[0] == 'hybrid':
        _, aa1, _ = g1_key
        _, aa2, _ = g2_key
        if aa1 and aa2 and aa_contexts_match(aa1, aa2):
            return 'ancestral', 'ancestral_aa_match'

    return 'independent', 'close_diff_aa_ctx'


def classify_ortholog_group(members, codon_tolerance=CODON_TOLERANCE):
    """Classify the sharing status of an ortholog group.

    Returns (within_group_status, cross_group_status, cross_group_reason) tuple.

    within_group_status: 'consistent', 'discordant', 'uncertain', or 'singleton'
    cross_group_status:  'ancestral', 'independent', 'uncertain', or 'NA'
    cross_group_reason:  fine-grained reason for the cross_group_status.
                         Empty string when cross_group_status is NA.
    """
    n_present = len(members)

    if n_present == 0:
        return ('uncertain', 'NA', '')
    if n_present == 1:
        return ('singleton', 'NA', '')

    typed = [m for m in members if m['confidence'] == 'high']

    # Determine group composition from ALL presence=1 members
    groups_present = {m['group'] for m in members}
    has_both = 'G1' in groups_present and 'G2' in groups_present

    # Within-group consistency (checked per-clade)
    within_status = _assess_within_group(typed, codon_tolerance)

    # Cross-group status
    if has_both:
        cross_status, cross_reason = _assess_cross_group(typed, codon_tolerance)
    else:
        cross_status, cross_reason = 'NA', ''

    return (within_status, cross_status, cross_reason)


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
    within_counts = defaultdict(int)
    cross_counts = defaultdict(int)
    reason_counts = defaultdict(int)
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
                'flanking_aa_context': fp.get('flanking_aa_context', ''),
                'confidence': fp.get('confidence', 'no_fingerprint'),
                'gene_id': fp.get('gene_id', ''),
            })

        # Classify
        within_status, cross_status, cross_reason = classify_ortholog_group(
            members, args.codon_tolerance)
        within_counts[within_status] += 1
        cross_counts[cross_status] += 1
        if cross_reason:
            reason_counts[cross_reason] += 1

        # Store summary for this group
        if members:
            typed = [m for m in members if m['confidence'] == 'high']

            def format_locus(m):
                ctx = m.get('flanking_aa_context', '')
                if ctx:
                    return f"{m['sample']}:{ctx}"
                if m['location_type'] == 'cds':
                    return f"{m['sample']}:cds.{m['codon_number']}.{m['codon_offset']}"
                if m['location_type'] == 'intron':
                    return f"{m['sample']}:intron.{m['intron_number']}"
                return f"{m['sample']}:?"

            group_summaries.append({
                'ortholog_id': oid,
                'within_group_status': within_status,
                'cross_group_status': cross_status,
                'cross_group_reason': cross_reason,
                'n_present': len(members),
                'n_fingerprinted': len(typed),
                'families': ';'.join(sorted(set(m['family'] for m in members))),
                'loci': ';'.join(format_locus(m) for m in typed),
            })

        # Assign status to all rows in this group
        for idx in row_indices:
            matrix_rows[idx]['within_group_status'] = within_status
            matrix_rows[idx]['cross_group_status'] = cross_status
            matrix_rows[idx]['cross_group_reason'] = cross_reason

    # Print summary
    print(f"\nWithin-group status (codon_tolerance={args.codon_tolerance}):")
    for status, count in sorted(within_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")

    print(f"\nCross-group status:")
    for status, count in sorted(cross_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / max(1, len(ortholog_groups))
        print(f"  {status:30s}  {count:5d}  ({pct:.1f}%)")

    if reason_counts:
        print(f"\nCross-group reason breakdown:")
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"  {reason:30s}  {count:5d}")

    print(f"  {'TOTAL':30s}  {len(ortholog_groups):5d}")

    # Write output matrix
    # Remove old sharing_status if present, add new columns
    out_fieldnames = [f for f in fieldnames
                      if f not in ('sharing_status', 'within_group_status',
                                   'cross_group_status', 'cross_group_reason')]
    out_fieldnames += ['within_group_status', 'cross_group_status',
                       'cross_group_reason']

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(matrix_rows)

    print(f"\nWrote {len(matrix_rows)} rows to {args.output}")

    # Write optional summary
    if args.summary:
        summary_fields = ['ortholog_id', 'within_group_status',
                          'cross_group_status', 'cross_group_reason',
                          'n_present', 'n_fingerprinted', 'families', 'loci']
        with open(args.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(group_summaries)
        print(f"Wrote {len(group_summaries)} group summaries to {args.summary}")


if __name__ == '__main__':
    main()

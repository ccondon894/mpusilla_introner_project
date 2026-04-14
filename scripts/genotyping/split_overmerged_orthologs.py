#!/usr/bin/env python3
"""
Split over-merged ortholog groups using codon/intron position clusters.

For each ortholog group with multiple distinct codon/intron position clusters
within Group 1 or Group 2 (within-group over-merge), partitions the present
members into clusters and creates a separate ortholog sub-group for each.

Lost introner recovery: when a sub-group is created, the splitter searches
each sample's BED file for additional introners that match the sub-group's
consensus locus. If found, those "lost" introners are recovered as presence=1
rows in the appropriate sub-group. If not found, the sample is marked as
presence=2 (genuinely absent at this sub-locus).

This ensures each sub-group has a complete sample set with biologically
correct presence/absence calls — not just the introners that happened to
make it through the original flank-based ortholog detection.

Cross-group differences (G1 vs G2 at distinct positions in the same gene
region) are NOT split — they represent meaningful biology (independent
insertions at the same approximate locus).

Reads:
  --matrix:       Verified genotype matrix (with within_group_status column)
  --fingerprints: Per-sample fingerprint TSV files
  --beds:         Per-sample BED files (for coordinate/splice/orientation lookup)
  --tandems:      Per-sample tandem cluster TSV files (optional, for flagging)

Writes:
  --output:       Split genotype matrix
  --mapping:      Original-to-new ortholog ID mapping
"""

import argparse
import csv
import os
from collections import Counter, defaultdict


GROUP2_SAMPLES = {'RCC1749', 'RCC3052'}
CODON_TOLERANCE = 3
AA_CONTEXT_WINDOW_LEN = 6
AA_CONTEXT_MISMATCH_TOLERANCE = 1
FLANK_LENGTH = 100  # used to back-compute BED-style coordinates from body span


def get_locus_key(fp):
    """Build a hybrid comparison key from a fingerprint dict.

    Returns a ('hybrid', aa_ctx_or_none, legacy_key_or_none) tuple.
    Two hybrid keys match if EITHER their aa contexts match OR their
    legacy codon/intron keys are compatible.
    """
    aa_ctx = fp.get('flanking_aa_context', '') or None

    legacy_key = None
    loc = fp.get('location_type', '')
    if loc == 'cds':
        try:
            legacy_key = ('cds', int(fp['codon_number']),
                          int(fp.get('codon_offset', 0)),
                          int(fp.get('exon_number', 0)))
        except (ValueError, TypeError, KeyError):
            pass
    elif loc == 'intron':
        try:
            legacy_key = ('intron', int(fp['intron_number']))
        except (ValueError, TypeError, KeyError):
            pass

    if aa_ctx is None and legacy_key is None:
        return None
    return ('hybrid', aa_ctx, legacy_key)


def aa_contexts_match(ctx1, ctx2,
                       window_len=AA_CONTEXT_WINDOW_LEN,
                       max_mismatches=AA_CONTEXT_MISMATCH_TOLERANCE):
    """Check if two aa context strings match via sliding-window comparison.

    See classify_sharing_status.aa_contexts_match for details.
    """
    if not ctx1 or not ctx2:
        return False
    if len(ctx1) < window_len or len(ctx2) < window_len:
        n = min(len(ctx1), len(ctx2))
        mm = sum(1 for i in range(n) if ctx1[i] != ctx2[i])
        return mm <= max_mismatches

    for i in range(len(ctx1) - window_len + 1):
        w1 = ctx1[i:i + window_len]
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


def _legacy_keys_match(lk1, lk2, codon_tolerance):
    """Check two legacy codon/intron keys for compatibility."""
    if lk1 is None or lk2 is None:
        return False
    if lk1[0] == lk2[0]:
        if lk1[0] == 'intron':
            return lk1[1] == lk2[1]
        if lk1[0] == 'cds':
            return abs(lk1[1] - lk2[1]) <= codon_tolerance
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

    Each key is a ('hybrid', aa_ctx, legacy_key) tuple. Two keys match if
    either the aa contexts match via sliding window OR the legacy codon/
    intron keys are compatible.
    """
    if k1 is None or k2 is None:
        return False

    # Handle legacy direct keys if any sneak through
    if k1[0] != 'hybrid' or k2[0] != 'hybrid':
        return _legacy_keys_match(k1, k2, codon_tolerance)

    _, aa1, legacy1 = k1
    _, aa2, legacy2 = k2

    # Try aa context match first
    if aa1 and aa2:
        if aa_contexts_match(aa1, aa2, max_mismatches=aa_mismatch_tolerance):
            return True

    # Fall back to legacy key match
    return _legacy_keys_match(legacy1, legacy2, codon_tolerance)


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
    """Load fingerprints into multiple indices for efficient lookup."""
    by_coord = {}
    by_sample_gene = defaultdict(list)
    for fp_file in fingerprint_files:
        with open(fp_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                key = (row['sample'], row['contig'],
                       int(row['start']), int(row['end']))
                by_coord[key] = row
                gene = row.get('gene_id', '')
                if gene and row.get('confidence') == 'high':
                    by_sample_gene[(row['sample'], gene)].append(row)
    return by_coord, by_sample_gene


def load_beds(bed_files):
    """Load BED entries indexed by (sample, contig, start, end).

    The sample is parsed from the filename (e.g. CCMP1545.candidate_loci.filtered.bed).
    """
    by_coord = {}
    for bed_file in bed_files:
        sample = os.path.basename(bed_file).split('.')[0]
        with open(bed_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                key = (sample, row['#chrom'],
                       int(row['start']), int(row['end']))
                by_coord[key] = row
    return by_coord


def load_tandems(tandem_files):
    """Load tandem cluster info."""
    tandems = {}
    for t_file in tandem_files:
        with open(t_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if row.get('tandem_cluster_id'):
                    key = (row['sample'], row['contig'],
                           int(row['start']), int(row['end']))
                    tandems[key] = row['tandem_cluster_id']
    return tandems


def load_gtf_cds(gtf_files):
    """Load CDS exons by (sample, gene_id) for absent-sample coordinate lookup.

    Returns dict: (sample, gene_id) -> list of CDS dicts with contig, start,
    end, strand. Sample is parsed from the filename.
    """
    cds_by_sample_gene = defaultdict(list)
    for gtf_file in gtf_files:
        sample = os.path.basename(gtf_file).split('.')[0]
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'CDS':
                    continue
                gene_id = None
                for attr in fields[8].split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                        break
                if gene_id:
                    cds_by_sample_gene[(sample, gene_id)].append({
                        'contig': fields[0],
                        'start': int(fields[3]),  # 1-based GTF
                        'end': int(fields[4]),
                        'strand': fields[6],
                    })
    return cds_by_sample_gene


def locate_consensus_in_gtf(sample, gene_id, locus_key, cds_by_sample_gene):
    """Find the genomic coordinates of a consensus locus in a sample's GTF.

    Args:
        sample: sample name
        gene_id: target gene
        locus_key: ('cds', codon, offset, exon) or ('intron', intron_num)
        cds_by_sample_gene: GTF index

    Returns:
        (contig, start, end, strand) where start/end are 0-based BED coords
        for a 305bp window matching the typical introner BED format
        (200bp body + 100bp flanks on each side, but here just 200bp body
        plus a 100bp window on each side as a placeholder for visualization).
        Returns None if the gene is not in the sample or the locus can't be
        located.
    """
    cds = cds_by_sample_gene.get((sample, gene_id), [])
    if not cds:
        return None

    strand = cds[0]['strand']
    if strand == '+':
        sorted_exons = sorted(cds, key=lambda x: x['start'])
    else:
        sorted_exons = sorted(cds, key=lambda x: x['start'], reverse=True)

    if locus_key[0] == 'intron':
        intron_num = locus_key[1]  # 1-indexed
        # Intron N is between transcript exon N and exon N+1
        if intron_num < 1 or intron_num >= len(sorted_exons):
            return None
        exon_before = sorted_exons[intron_num - 1]
        exon_after = sorted_exons[intron_num]

        if strand == '+':
            # Intron is between exon_before.end and exon_after.start
            intron_start = exon_before['end']  # 1-based, last exonic nt
            intron_end = exon_after['start']   # 1-based, first exonic nt
        else:
            # On - strand, exon_before is at higher genomic coords
            intron_start = exon_after['end']
            intron_end = exon_before['start']

        # Build a BED-style window around the intron region
        # Use a flank-padded window of the intron
        bed_start = max(0, intron_start - 1 - FLANK_LENGTH)
        bed_end = intron_end + FLANK_LENGTH
        return (sorted_exons[0]['contig'], bed_start, bed_end, strand)

    if locus_key[0] == 'cds':
        codon_num = locus_key[1]
        codon_offset = locus_key[2]
        target_cds_pos = codon_num * 3 + codon_offset

        cds_offset = 0
        for exon in sorted_exons:
            exon_len = exon['end'] - exon['start'] + 1
            if cds_offset + exon_len > target_cds_pos:
                pos_in_exon = target_cds_pos - cds_offset
                if strand == '+':
                    insertion_genomic_1 = exon['start'] + pos_in_exon
                else:
                    insertion_genomic_1 = exon['end'] - pos_in_exon
                # 1-based -> 0-based BED
                bed_start = max(0, insertion_genomic_1 - 1 - FLANK_LENGTH)
                bed_end = insertion_genomic_1 + FLANK_LENGTH
                return (exon['contig'], bed_start, bed_end, strand)
            cds_offset += exon_len
        return None

    return None


def find_matching_fingerprint(sample_gene_fps, target_key, codon_tolerance,
                                used_coords):
    """Search a sample's fingerprints for one matching a target locus key.

    Args:
        sample_gene_fps: list of fingerprint records for (sample, gene)
        target_key: locus key tuple to match against
        codon_tolerance: tolerance for the comparison
        used_coords: set of (contig, start, end) tuples already assigned to
                     other sub-groups (so we don't double-assign)

    Returns:
        The matching fingerprint record, or None.
    """
    for fp in sample_gene_fps:
        coord = (fp['contig'], int(fp['start']), int(fp['end']))
        if coord in used_coords:
            continue
        fp_key = get_locus_key(fp)
        if keys_within_tolerance(fp_key, target_key, codon_tolerance):
            return fp
    return None


def get_consensus_key_and_gene(cluster_indices, keys, members):
    """Get the consensus key, gene_id, and a legacy location key for a cluster.

    Returns:
        consensus_key: hybrid key for matching via keys_within_tolerance
        consensus_gene: most common gene_id
        fallback_location_key: the legacy codon/intron part of the consensus
            hybrid key, used for coordinate lookup in absent samples (since
            aa_ctx doesn't directly map to genomic positions).
    """
    cluster_keys = [keys[i] for i in cluster_indices if keys[i] is not None]
    cluster_members = [members[i] for i in cluster_indices]
    cluster_genes = [m.get('gene_id', '') for m in cluster_members]

    if not cluster_keys:
        return None, '', None

    # Hybrid keys have the form ('hybrid', aa_ctx, legacy_key). Pick the most
    # common one for the consensus.
    key_counts = Counter(cluster_keys)
    gene_counts = Counter(g for g in cluster_genes if g)
    consensus_key = key_counts.most_common(1)[0][0]
    consensus_gene = gene_counts.most_common(1)[0][0] if gene_counts else ''

    # Extract the legacy part of the consensus key, or fall back to
    # searching for any member with a legacy key.
    fallback_key = None
    if consensus_key and consensus_key[0] == 'hybrid' and consensus_key[2] is not None:
        fallback_key = consensus_key[2]
    else:
        for k in cluster_keys:
            if k and k[0] == 'hybrid' and k[2] is not None:
                fallback_key = k[2]
                break

    return consensus_key, consensus_gene, fallback_key


def build_recovered_row(template_row, bed_row, fp_row, presence='1'):
    """Build a new matrix row from BED + fingerprint data, using template
    for the column structure and inheriting non-coordinate fields."""
    new_row = dict(template_row)
    new_row['contig'] = bed_row['#chrom']
    new_row['start'] = bed_row['start']
    new_row['end'] = bed_row['end']
    new_row['sequence_id'] = f"{bed_row['#chrom']}:{bed_row['start']}-{bed_row['end']}"
    new_row['family'] = bed_row.get('family', '')
    new_row['gene'] = bed_row.get('gene_info', '')
    new_row['splice_site'] = bed_row.get('splice_info', '')
    new_row['orientation'] = bed_row.get('orientation', 'forward')
    new_row['presence'] = presence
    # left/right_reverse will be re-computed by check_orientation if re-run;
    # default to False here
    new_row['left_reverse'] = 'False'
    new_row['right_reverse'] = 'False'
    return new_row


def build_absent_row(template_row, presence='2', locus_coords=None,
                      consensus_gene=''):
    """Build a new 'absent' or 'missing' row for a sample at a sub-locus.

    For presence=2 (absent): records the genomic coordinates where the
    introner WOULD have been if it were present, looked up via the consensus
    locus position in the sample's GTF. This makes the absent calls
    biologically locatable rather than just marked as missing.

    For presence=3 (missing data): leaves coordinates blank.

    Args:
        template_row: an existing matrix row from this sample to inherit
            sample/identity fields from
        presence: '2' for absent, '3' for missing
        locus_coords: tuple (contig, bed_start, bed_end, strand) from
            locate_consensus_in_gtf, or None
        consensus_gene: the gene_id at the consensus locus
    """
    new_row = dict(template_row)
    new_row['presence'] = presence
    new_row['family'] = ''
    new_row['splice_site'] = ''
    new_row['left_reverse'] = ''
    new_row['right_reverse'] = ''

    if presence == '2' and locus_coords is not None:
        contig, bed_start, bed_end, strand = locus_coords
        new_row['contig'] = contig
        new_row['start'] = bed_start
        new_row['end'] = bed_end
        new_row['sequence_id'] = f"{contig}:{bed_start}-{bed_end}"
        new_row['gene'] = consensus_gene
        new_row['orientation'] = 'forward' if strand == '+' else 'reverse'
    else:
        # presence=3 (missing data) or no GTF lookup possible
        new_row['contig'] = ''
        new_row['start'] = ''
        new_row['end'] = ''
        new_row['sequence_id'] = ''
        new_row['gene'] = ''
        new_row['orientation'] = ''

    return new_row


def split_group_with_recovery(oid, row_indices, matrix_rows, fps_by_coord,
                                fps_by_sample_gene, beds_by_coord, tandems,
                                cds_by_sample_gene, codon_tolerance):
    """Split an ortholog group, recovering lost introners.

    Returns a list of (new_id, [matrix_rows]) tuples. Each tuple represents
    a sub-group with a complete sample set.
    """
    # Step 1: Get presence=1 members and their locus keys
    present_members = []  # list of dicts with idx, sample, fp, key
    for idx in row_indices:
        row = matrix_rows[idx]
        if row['presence'] != '1':
            continue
        coord = (row['sample'], row['contig'],
                 int(row['start']), int(row['end']))
        fp = fps_by_coord.get(coord, {})
        key = get_locus_key(fp)
        present_members.append({
            'idx': idx,
            'sample': row['sample'],
            'fp': fp,
            'key': key,
            'gene_id': fp.get('gene_id', ''),
        })

    # Filter to those with valid fingerprints for clustering
    typed = [m for m in present_members if m['key'] is not None]

    if len(typed) < 2:
        return None  # not splittable

    # Step 2: Check for within-group multi-clustering
    g1_keys = [m['key'] for m in typed if m['sample'] not in GROUP2_SAMPLES]
    g2_keys = [m['key'] for m in typed if m['sample'] in GROUP2_SAMPLES]
    g1_clusters = cluster_keys(g1_keys, codon_tolerance) if g1_keys else []
    g2_clusters = cluster_keys(g2_keys, codon_tolerance) if g2_keys else []

    if len(g1_clusters) <= 1 and len(g2_clusters) <= 1:
        return None  # no within-group over-merge

    # Step 3: Compute the global cluster assignment
    keys = [m['key'] for m in typed]
    clusters = cluster_keys(keys, codon_tolerance)

    if len(clusters) < 2:
        return None  # nothing to split

    # Map from member idx -> cluster index
    member_to_cluster = {}
    for ci, cluster in enumerate(clusters):
        for ki in cluster:
            member_to_cluster[typed[ki]['idx']] = ci

    # Step 4: For each cluster, get consensus locus key, gene, and fallback
    # legacy key (used for coordinate lookup in absent samples, since the
    # aa_ctx key doesn't map directly to genomic positions).
    sub_consensuses = []
    for ci, cluster in enumerate(clusters):
        ck, cgene, fallback = get_consensus_key_and_gene(cluster, keys, typed)
        sub_consensuses.append({
            'cluster_idx': ci,
            'key': ck,
            'gene': cgene,
            'fallback_key': fallback,
        })

    # Step 5: Get all samples that were in the original group
    all_samples = []
    seen_samples = set()
    sample_to_template = {}
    for idx in row_indices:
        s = matrix_rows[idx]['sample']
        if s not in seen_samples:
            seen_samples.add(s)
            all_samples.append(s)
            sample_to_template[s] = matrix_rows[idx]

    # Step 6: Track which BED entries are already assigned to a sub-group
    used_coords_per_sample = defaultdict(set)

    # First pass: assign existing presence=1 rows to their cluster
    sub_group_rows = [[] for _ in clusters]
    samples_with_pres1_in_cluster = [set() for _ in clusters]

    for m in present_members:
        if m['idx'] in member_to_cluster:
            ci = member_to_cluster[m['idx']]
            row = dict(matrix_rows[m['idx']])
            sub_group_rows[ci].append(row)
            samples_with_pres1_in_cluster[ci].add(m['sample'])
            coord = (matrix_rows[m['idx']]['contig'],
                     int(matrix_rows[m['idx']]['start']),
                     int(matrix_rows[m['idx']]['end']))
            used_coords_per_sample[m['sample']].add(coord)
        else:
            # presence=1 but no fingerprint — assign to first cluster as fallback
            row = dict(matrix_rows[m['idx']])
            sub_group_rows[0].append(row)
            samples_with_pres1_in_cluster[0].add(m['sample'])
            coord = (matrix_rows[m['idx']]['contig'],
                     int(matrix_rows[m['idx']]['start']),
                     int(matrix_rows[m['idx']]['end']))
            used_coords_per_sample[m['sample']].add(coord)

    # Step 7: For each cluster, fill in missing samples
    has_tandem_flags = [False] * len(clusters)

    for ci, consensus in enumerate(sub_consensuses):
        ckey = consensus['key']
        cgene = consensus['gene']
        fallback_key = consensus.get('fallback_key')

        # Check if any presence=1 row in this cluster has a tandem flag
        for row in sub_group_rows[ci]:
            coord = (row['sample'], row['contig'],
                     int(row['start']), int(row['end']))
            if coord in tandems:
                has_tandem_flags[ci] = True
                break

        for sample in all_samples:
            if sample in samples_with_pres1_in_cluster[ci]:
                continue  # already has presence=1 row in this cluster

            template = sample_to_template[sample]

            # Try to find a matching fingerprint in the sample's BED
            matching_fp = None
            if cgene:
                sample_gene_fps = fps_by_sample_gene.get((sample, cgene), [])
                matching_fp = find_matching_fingerprint(
                    sample_gene_fps, ckey, codon_tolerance,
                    used_coords_per_sample[sample])

            if matching_fp:
                # Recovered a "lost" introner
                bed_coord = (sample, matching_fp['contig'],
                             int(matching_fp['start']), int(matching_fp['end']))
                bed_row = beds_by_coord.get(bed_coord)
                if bed_row:
                    new_row = build_recovered_row(template, bed_row,
                                                   matching_fp, presence='1')
                    sub_group_rows[ci].append(new_row)
                    used_coords_per_sample[sample].add(
                        (matching_fp['contig'], int(matching_fp['start']),
                         int(matching_fp['end'])))
                    if bed_coord in tandems:
                        has_tandem_flags[ci] = True
                    continue

            # No recovery possible — mark as absent at this sub-locus.
            # Inherit original presence value: if it was 3 (missing), keep 3.
            # For presence=2, look up the consensus locus in this sample's
            # GTF to record where the introner WOULD have been. Uses the
            # fallback legacy key (codon/intron) since locate_consensus_in_gtf
            # needs a codon/intron position, not an aa context string.
            original_presence = template['presence']
            if original_presence == '3':
                new_row = build_absent_row(template, presence='3')
            else:
                locus_coords = None
                if cgene and fallback_key is not None:
                    locus_coords = locate_consensus_in_gtf(
                        sample, cgene, fallback_key, cds_by_sample_gene)
                new_row = build_absent_row(
                    template, presence='2',
                    locus_coords=locus_coords,
                    consensus_gene=cgene)
            sub_group_rows[ci].append(new_row)

    # Step 8: Build sub-group results
    result = []
    for ci, rows in enumerate(sub_group_rows):
        new_id = f"{oid}_split{ci + 1}"
        for row in rows:
            row['ortholog_id'] = new_id
        result.append((new_id, rows, len(clusters), has_tandem_flags[ci]))

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Split over-merged ortholog groups with lost introner recovery')
    parser.add_argument('--matrix', required=True,
                        help='Verified genotype matrix with within_group_status column')
    parser.add_argument('--fingerprints', required=True, nargs='+',
                        help='Per-sample fingerprint TSV files')
    parser.add_argument('--beds', required=True, nargs='+',
                        help='Per-sample BED files (candidate_loci.filtered.bed)')
    parser.add_argument('--gtfs', required=True, nargs='+',
                        help='Per-sample GTF files (for absent-row coordinate lookup)')
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
    fps_by_coord, fps_by_sample_gene = load_fingerprints(args.fingerprints)
    print(f"  Loaded {len(fps_by_coord)} fingerprint records")
    print(f"  Indexed across {len(fps_by_sample_gene)} (sample, gene) pairs")

    print(f"Loading BED files from {len(args.beds)} files...")
    beds_by_coord = load_beds(args.beds)
    print(f"  Loaded {len(beds_by_coord)} BED records")

    print(f"Loading GTF files from {len(args.gtfs)} files...")
    cds_by_sample_gene = load_gtf_cds(args.gtfs)
    print(f"  Loaded CDS for {len(cds_by_sample_gene)} (sample, gene) pairs")

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
            row['original_ortholog_id'] = row['ortholog_id']
            matrix_rows.append(row)
            ortholog_groups[row['ortholog_id']].append(idx)
    print(f"  {len(matrix_rows)} rows in {len(ortholog_groups)} ortholog groups")

    # Process each group
    n_split = 0
    n_subgroups_added = 0
    n_tandem_flagged = 0
    n_recovered = 0
    output_rows = []
    mapping_records = []

    for oid in sorted(ortholog_groups.keys()):
        row_indices = ortholog_groups[oid]
        result = split_group_with_recovery(
            oid, row_indices, matrix_rows, fps_by_coord, fps_by_sample_gene,
            beds_by_coord, tandems, cds_by_sample_gene, args.codon_tolerance)

        if result is None:
            # No split — keep original rows
            for idx in row_indices:
                output_rows.append(matrix_rows[idx])
            n_present = sum(1 for idx in row_indices
                            if matrix_rows[idx]['presence'] == '1')
            mapping_records.append({
                'original_id': oid,
                'new_id': oid,
                'n_clusters': 1,
                'n_members': n_present,
                'n_recovered': 0,
                'has_tandem': '',
                'note': 'unchanged',
            })
            continue

        # Splitting happened
        n_split += 1
        n_clusters = len(result)
        n_subgroups_added += n_clusters - 1

        # Count recovered introners (rows that didn't exist in original matrix)
        original_pres1_count = sum(1 for idx in row_indices
                                   if matrix_rows[idx]['presence'] == '1')
        new_pres1_count = sum(1 for _, rows, _, _ in result
                              for r in rows if r['presence'] == '1')
        recovered_in_group = new_pres1_count - original_pres1_count
        n_recovered += recovered_in_group

        for new_id, rows, _, has_tandem in result:
            output_rows.extend(rows)
            if has_tandem:
                n_tandem_flagged += 1
            mapping_records.append({
                'original_id': oid,
                'new_id': new_id,
                'n_clusters': n_clusters,
                'n_members': sum(1 for r in rows if r['presence'] == '1'),
                'n_recovered': recovered_in_group if new_id == result[0][0] else '',
                'has_tandem': '1' if has_tandem else '0',
                'note': 'split',
            })

    # Print summary
    print(f"\n=== Splitting summary ===")
    print(f"  Original groups:                 {len(ortholog_groups)}")
    print(f"  Groups split:                    {n_split}")
    print(f"  New sub-groups added:            {n_subgroups_added}")
    print(f"  Total groups after splitting:    {len(ortholog_groups) + n_subgroups_added}")
    print(f"  Recovered introners (presence=1): {n_recovered}")
    print(f"  Sub-groups with tandem flag:     {n_tandem_flagged}")
    print(f"  Total output rows:               {len(output_rows)}")
    print(f"  (Original matrix rows: {len(matrix_rows)})")

    # Write output matrix
    out_fieldnames = list(fieldnames)
    if 'original_ortholog_id' not in out_fieldnames:
        out_fieldnames.append('original_ortholog_id')

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames,
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"\nWrote {len(output_rows)} rows to {args.output}")

    # Write mapping
    if args.mapping:
        mapping_fields = ['original_id', 'new_id', 'n_clusters', 'n_members',
                          'n_recovered', 'has_tandem', 'note']
        with open(args.mapping, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=mapping_fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(mapping_records)
        print(f"Wrote {len(mapping_records)} mapping records to {args.mapping}")


if __name__ == '__main__':
    main()

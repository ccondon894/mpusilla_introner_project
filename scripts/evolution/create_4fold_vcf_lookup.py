#!/usr/bin/env python3
"""
VCF-based 4-fold degenerate site diversity calculation for haplotype diversity analysis.

Parses the 4-fold degenerate VCF file to calculate sample-set-specific π values
by computing diversity site-by-site for each unique sample set, then averaging.

Optimizations over original version:
- Interval tree for introner region overlap checking (O(log n) per site)
- Vectorized π calculation using numpy arrays instead of per-site Python loops
"""

import argparse
import pandas as pd
import numpy as np
from bisect import bisect_left, bisect_right
from collections import defaultdict
import gzip
import os


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create 4-fold diversity lookup tables from VCF file')
    parser.add_argument('--vcf', required=True, help='Input VCF file with 4-fold degenerate sites')
    parser.add_argument('--genotype_matrix', required=True, help='Genotype matrix TSV file')
    parser.add_argument('--output_dir', required=True, help='Output directory for lookup tables')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    return parser.parse_args()


def format_sample_set(sample_list):
    """Convert list of samples to sorted set format with curly braces."""
    if not sample_list:
        return ''
    return '{' + ','.join(sorted(sample_list)) + '}'


def build_introner_region_index(genotype_matrix_file, verbose=False):
    """
    Build a sorted interval index of CCMP1545 introner regions for fast overlap queries.

    Returns dict: chrom -> sorted list of (start, end) tuples
    """
    if verbose:
        print(f"Loading CCMP1545 introner regions from {genotype_matrix_file}")

    df = pd.read_csv(genotype_matrix_file, sep='\t')
    ccmp = df[(df['sample'] == 'CCMP1545') & (df['presence'] != 3)]

    if verbose:
        print(f"  Found {len(ccmp)} CCMP1545 introner regions")

    index = defaultdict(list)
    for _, row in ccmp.iterrows():
        index[row['contig']].append((int(row['start']), int(row['end'])))

    for chrom in index:
        index[chrom].sort()

    return index


def site_overlaps_introner(chrom, pos, introner_index):
    """Check if a site overlaps any introner region using bisect for O(log n) lookup."""
    regions = introner_index.get(chrom)
    if not regions:
        return False

    # Find the rightmost region whose start <= pos
    idx = bisect_right(regions, (pos, float('inf'))) - 1

    if idx >= 0 and regions[idx][0] <= pos <= regions[idx][1]:
        return True

    return False


def extract_unique_sample_sets(genotype_matrix_file, verbose=False):
    """Extract all unique sample sets from the genotype matrix."""
    if verbose:
        print(f"Loading genotype matrix from {genotype_matrix_file}")

    df = pd.read_csv(genotype_matrix_file, sep='\t')

    group2_samples = {'RCC1749', 'RCC3052'}
    df_filtered = df[~df['sample'].isin(group2_samples)]

    sample_sets = set()

    # Only include truly polymorphic loci (both present and absent calls)
    for ortholog_id, ortholog_data in df_filtered.groupby('ortholog_id'):
        presences = set(ortholog_data['presence'])
        if 1 not in presences or 2 not in presences:
            continue  # Skip monomorphic or not-callable-only loci

        present_samples = ortholog_data[ortholog_data['presence'] == 1]['sample'].tolist()
        absent_samples = ortholog_data[ortholog_data['presence'] == 2]['sample'].tolist()

        if len(present_samples) >= 2:
            sample_sets.add(format_sample_set(present_samples))
        if len(absent_samples) >= 2:
            sample_sets.add(format_sample_set(absent_samples))

    if verbose:
        print(f"  Found {len(sample_sets)} unique sample sets with >= 2 samples (polymorphic loci only)")

    return sample_sets


def parse_vcf_header(file_handle):
    """Parse VCF header to get sample names. Returns sample names list."""
    for line in file_handle:
        if isinstance(line, bytes):
            line = line.decode()
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            return fields[9:]  # Sample names start at column 9
        if not line.startswith('#'):
            break
    return None


def process_vcf_vectorized(vcf_file, sample_sets, introner_index, verbose=False):
    """
    Process VCF file and compute π for each sample set using vectorized operations.

    Instead of computing π per-site in Python, we:
    1. Parse all genotypes into a numpy matrix (sites x samples)
    2. For each sample set, compute π across all sites at once
    """
    if verbose:
        print(f"Processing VCF file: {vcf_file}")

    # Parse sample set strings into sample lists and indices
    sample_set_info = {}
    for ss_str in sample_sets:
        sample_list = ss_str.strip('{}').split(',')
        sample_set_info[ss_str] = sample_list

    # First pass: parse VCF into arrays
    if vcf_file.endswith('.gz'):
        fh = gzip.open(vcf_file, 'rt')
    else:
        fh = open(vcf_file, 'r')

    # Get sample names from header
    sample_names = None
    for line in fh:
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            sample_names = fields[9:]
            break
        if not line.startswith('#'):
            break

    if sample_names is None:
        raise ValueError("Could not parse VCF header")

    if verbose:
        print(f"  VCF samples: {sample_names}")

    n_samples = len(sample_names)
    sample_to_idx = {s: i for i, s in enumerate(sample_names)}

    # Pre-compute sample index arrays for each sample set
    sample_set_indices = {}
    for ss_str, sample_list in sample_set_info.items():
        indices = [sample_to_idx[s] for s in sample_list if s in sample_to_idx]
        if len(indices) >= 2:
            sample_set_indices[ss_str] = np.array(indices, dtype=np.int32)

    # Parse genotypes into lists (will convert to numpy after)
    genotype_rows = []
    sites_excluded = 0
    sites_processed = 0

    for line in fh:
        if line.startswith('#'):
            continue

        fields = line.split('\t', 10)  # Only split enough to get chrom, pos, and start of genotypes
        chrom = fields[0]
        pos = int(fields[1])

        # Fast overlap check
        if site_overlaps_introner(chrom, pos, introner_index):
            sites_excluded += 1
            continue

        # Parse genotypes for all samples
        gt_fields = line.strip().split('\t')[9:]
        row = np.full(n_samples, -1, dtype=np.int8)  # -1 = missing

        for i, gt_field in enumerate(gt_fields):
            gt = gt_field.split(':')[0]
            if gt in ('0', '0/0', '0|0'):
                row[i] = 0
            elif gt in ('1', '1/1', '1|1'):
                row[i] = 1

        genotype_rows.append(row)
        sites_processed += 1

        if verbose and sites_processed % 500000 == 0:
            print(f"  Parsed {sites_processed} sites...")

    fh.close()

    if verbose:
        print(f"  Parsed {sites_processed} sites, excluded {sites_excluded} overlapping introner regions")

    if not genotype_rows:
        return {}

    # Convert to numpy matrix: (n_sites x n_samples)
    genotype_matrix = np.array(genotype_rows, dtype=np.int8)
    n_sites = genotype_matrix.shape[0]

    if verbose:
        print(f"  Genotype matrix: {n_sites} sites x {n_samples} samples")
        print(f"  Computing π for {len(sample_set_indices)} sample sets...")

    # Vectorized π calculation for each sample set
    results = {}
    for ss_str, indices in sample_set_indices.items():
        # Extract genotype submatrix for this sample set
        sub = genotype_matrix[:, indices]  # (n_sites x n_subset_samples)
        n = len(indices)

        # Mask missing data (-1)
        valid = sub >= 0  # boolean mask
        n_valid = valid.sum(axis=1)  # valid samples per site

        # Count alt alleles (value == 1) where valid
        alt_counts = (sub == 1).sum(axis=1).astype(np.float64)

        # Allele frequencies (only where n_valid >= 2)
        usable = n_valid >= 2
        p = np.zeros(n_sites)
        q = np.zeros(n_sites)

        n_usable = n_valid[usable].astype(np.float64)
        p[usable] = alt_counts[usable] / n_usable
        q[usable] = 1.0 - p[usable]

        # π = (n/(n-1)) * 2pq per site
        pi_per_site = np.zeros(n_sites)
        pi_per_site[usable] = (n_usable / (n_usable - 1)) * 2 * p[usable] * q[usable]

        # Average π across all usable sites
        n_usable_sites = usable.sum()
        if n_usable_sites > 0:
            avg_pi = pi_per_site[usable].mean()
            std_pi = pi_per_site[usable].std()
        else:
            avg_pi = 0.0
            std_pi = 0.0

        results[ss_str] = {
            'pi_4D': avg_pi,
            'site_count': int(n_usable_sites),
            'pi_std': std_pi
        }

    if verbose:
        print(f"  Computed π for {len(results)} sample sets")

    return results


def main():
    args = parse_arguments()

    if args.verbose:
        print("VCF-based 4-fold degenerate site diversity calculation")
        print("=" * 60)

    # Build introner region index for fast overlap checking
    introner_index = build_introner_region_index(args.genotype_matrix, args.verbose)

    # Extract unique sample sets
    sample_sets = extract_unique_sample_sets(args.genotype_matrix, args.verbose)

    # Process VCF with vectorized computation
    if args.verbose:
        print("\nCalculating π values from VCF...")

    results = process_vcf_vectorized(args.vcf, sample_sets, introner_index, args.verbose)

    # Format output
    output_rows = []
    for ss_str in sorted(results.keys()):
        data = results[ss_str]
        output_rows.append({
            'sample_set': ss_str,
            'pi_4D': data['pi_4D'],
            'site_count': data['site_count'],
            'pi_std': data['pi_std']
        })

    # Save
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, '4fold_avg_pi_present_samples.tsv')
    results_df = pd.DataFrame(output_rows)
    results_df.to_csv(output_file, sep='\t', index=False)

    if args.verbose:
        print(f"\nSaved {len(output_rows)} sample sets to {output_file}")
        top = sorted(output_rows, key=lambda x: x['site_count'], reverse=True)
        print("Top 10 sample sets by site count:")
        for i, row in enumerate(top[:10]):
            print(f"  {i+1:2d}. {row['sample_set']}: π={row['pi_4D']:.6f} (n={row['site_count']})")

    print(f"Generated VCF-based lookup table: {output_file}")
    print("Done!")


if __name__ == "__main__":
    main()

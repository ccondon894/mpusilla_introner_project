#!/usr/bin/env python3
"""
Prepare isoform count data for regression analysis.

Merges isoform counts, introner features (from genotype matrix),
CDS lengths (from GTF), and library-size normalized expression from the
featureCounts count matrix.

Usage (Snakemake):
    python prepare_isoform_data.py \
        --isoform_counts CSV --genotype_matrix TSV --gtf GTF \
        --output CSV --filtered CSV --counts_matrix CSV \
        --mt_scaffold STR --mt_start INT --mt_end INT
"""

import argparse
import pandas as pd
import numpy as np
import re
from collections import defaultdict


def extract_cds_lengths_from_gtf(gtf_path):
    """Extract total CDS length for each gene from GTF file."""
    gene_cds = defaultdict(int)

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'CDS':
                match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', fields[8])
                if match:
                    gene_id = match.group(1)
                    start = int(fields[3])
                    end = int(fields[4])
                    gene_cds[gene_id] += end - start + 1

    return dict(gene_cds)


def extract_mating_type_genes(gtf_path, scaffold, start, end):
    """Extract gene IDs in the mating-type region."""
    mt_genes = set()
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            if scaffold not in chrom:
                continue
            feat_start = int(fields[3])
            feat_end = int(fields[4])
            if feat_start <= end and feat_end >= start:
                match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', fields[8])
                if match:
                    mt_genes.add(match.group(1))
    return mt_genes


def derive_introner_features(genotype_matrix_path, reference='CCMP1545'):
    """Derive introner count, gain, loss features from genotype matrix."""
    gm = pd.read_csv(genotype_matrix_path, sep='\t')

    # Count introners per gene per sample
    introner_counts = (
        gm[gm['presence'] == 1]
        .groupby(['gene', 'sample'])
        .size()
        .reset_index(name='introner_count')
    )

    # Get reference counts
    ref_counts = introner_counts[introner_counts['sample'] == reference][
        ['gene', 'introner_count']
    ].rename(columns={'introner_count': 'baseline_introner_count'})

    # All gene-sample combinations
    all_genes = gm['gene'].unique()
    all_samples = gm['sample'].unique()
    all_combos = pd.DataFrame([
        {'gene': g, 'sample': s}
        for g in all_genes for s in all_samples
    ])

    features = all_combos.merge(introner_counts, on=['gene', 'sample'], how='left')
    features['introner_count'] = features['introner_count'].fillna(0).astype(int)

    features = features.merge(ref_counts, on='gene', how='left')
    features['baseline_introner_count'] = features['baseline_introner_count'].fillna(0).astype(int)

    # Gain = present in sample but not in reference
    # Loss = present in reference but not in sample
    features['introner_gain'] = np.maximum(
        features['introner_count'] - features['baseline_introner_count'], 0
    )
    features['introner_loss'] = np.maximum(
        features['baseline_introner_count'] - features['introner_count'], 0
    )
    features['introner_net_change'] = (
        features['introner_gain'] - features['introner_loss']
    )

    # Map sample names to strain names
    strain_map = {
        'CCMP1545': 'CCMP1545',
        'RCC1614': 'RCC1614',
        'RCC1749': 'RCC1749'
    }
    features['strain'] = features['sample'].map(strain_map)
    features = features.dropna(subset=['strain'])

    return features.rename(columns={'gene': 'gene_id'})[
        ['gene_id', 'strain', 'introner_count', 'baseline_introner_count',
         'introner_gain', 'introner_loss', 'introner_net_change']
    ]


def compute_mean_expression(counts_matrix_path):
    """
    Compute mean expression per gene per strain from featureCounts merged matrix.

    The negative-binomial expression model uses raw counts with log library size
    as an offset. For gene-level summaries and plotting, use the same rate logic
    by converting each replicate to counts per million mapped reads before
    averaging replicates within strain.

    Returns DataFrame with columns: gene_id, strain, mean_expression, log_expression
    """
    counts = pd.read_csv(counts_matrix_path, index_col=0)

    # Drop non-count columns
    drop_cols = [c for c in counts.columns
                 if c in ('Length', 'Chr', 'Start', 'End', 'Strand')]
    sample_cols = [c for c in counts.columns if c not in drop_cols]
    library_sizes = counts[sample_cols].sum()
    normalized_counts = counts[sample_cols].div(library_sizes, axis=1) * 1_000_000

    # Map replicate to strain
    strain_map = {}
    for col in sample_cols:
        if col.startswith('834'):
            strain_map[col] = 'CCMP1545'
        elif col.startswith('1614'):
            strain_map[col] = 'RCC1614'
        elif col.startswith('1749'):
            strain_map[col] = 'RCC1749'

    results = []
    for strain in ['CCMP1545', 'RCC1614', 'RCC1749']:
        strain_cols = [c for c in sample_cols if strain_map.get(c) == strain]
        if not strain_cols:
            continue
        mean_expr = normalized_counts[strain_cols].mean(axis=1)
        gene_expr = pd.DataFrame({
            'gene_id': counts.index,
            'strain': strain,
            'mean_expression': mean_expr.values
        })
        results.append(gene_expr)

    expr_df = pd.concat(results, ignore_index=True)
    expr_df['log_expression'] = np.log(expr_df['mean_expression'] + 1)
    return expr_df


def main():
    parser = argparse.ArgumentParser(
        description='Prepare isoform data for GLM modeling'
    )
    parser.add_argument('--isoform_counts', required=True,
                        help='Isoform counts per gene CSV (from isoform_abundance_analysis)')
    parser.add_argument('--genotype_matrix', required=True,
                        help='Genotype matrix TSV')
    parser.add_argument('--gtf', required=True,
                        help='GTF file for CDS length extraction')
    parser.add_argument('--output', required=True,
                        help='Full output CSV')
    parser.add_argument('--filtered', required=True,
                        help='Filtered output CSV (genes with <=20 isoforms)')
    parser.add_argument('--counts_matrix', default=None,
                        help='Merged featureCounts matrix CSV for expression')
    parser.add_argument('--mt_scaffold', required=True,
                        help='Mating-type scaffold name')
    parser.add_argument('--mt_start', type=int, required=True,
                        help='Mating-type region start')
    parser.add_argument('--mt_end', type=int, required=True,
                        help='Mating-type region end')

    args = parser.parse_args()

    print("=" * 80)
    print("Isoform Data Preparation")
    print("=" * 80)
    print()

    # 1. Load isoform counts
    print("Step 1: Loading isoform counts...")
    isoform_data = pd.read_csv(args.isoform_counts)
    print(f"  {len(isoform_data)} gene-strain combinations")

    # 2. Derive introner features from genotype matrix
    print("Step 2: Deriving introner features from genotype matrix...")
    introner_features = derive_introner_features(args.genotype_matrix)
    print(f"  {len(introner_features)} gene-strain-introner records")

    # 3. Extract CDS lengths from GTF
    print("Step 3: Extracting CDS lengths from GTF...")
    cds_lengths = extract_cds_lengths_from_gtf(args.gtf)
    cds_df = pd.DataFrame([
        {'gene_id': g, 'cds_length': l} for g, l in cds_lengths.items()
    ])
    print(f"  CDS lengths for {len(cds_df)} genes")

    # 4. Compute mean expression from counts matrix (if provided)
    if args.counts_matrix:
        print("Step 4: Computing mean expression from counts matrix...")
        expr_df = compute_mean_expression(args.counts_matrix)
        print(f"  Expression data for {len(expr_df)} gene-strain combinations")
    else:
        print("Step 4: No counts matrix provided, skipping expression")
        expr_df = None

    # 5. Merge everything
    print("Step 5: Merging datasets...")
    merged = isoform_data.merge(introner_features, on=['gene_id', 'strain'], how='left')
    merged = merged.merge(cds_df, on='gene_id', how='inner')

    if expr_df is not None:
        merged = merged.merge(expr_df, on=['gene_id', 'strain'], how='left')

    # Fill missing introner features with 0
    for col in ['introner_count', 'baseline_introner_count',
                'introner_gain', 'introner_loss', 'introner_net_change']:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    # Add log CDS length
    merged['log_cds_length'] = np.log(merged['cds_length'])

    # Add log expression if not already present
    if 'log_expression' not in merged.columns and 'mean_expression' in merged.columns:
        merged['log_expression'] = np.log(merged['mean_expression'] + 1)

    print(f"  Merged: {len(merged)} rows")

    # 6. Filter mating-type region
    print("Step 6: Filtering mating-type region...")
    mt_genes = extract_mating_type_genes(args.gtf, args.mt_scaffold,
                                          args.mt_start, args.mt_end)
    n_before = len(merged)
    merged = merged[~merged['gene_id'].isin(mt_genes)]
    print(f"  Removed {n_before - len(merged)} rows in MT region")

    # 7. Save full output
    print("Step 7: Saving output...")
    merged.to_csv(args.output, index=False)
    print(f"  Full dataset: {args.output} ({len(merged)} rows)")

    # 8. Save filtered output (<=20 isoforms)
    filtered = merged[merged['n_isoforms'] <= 20]
    filtered.to_csv(args.filtered, index=False)
    print(f"  Filtered dataset: {args.filtered} ({len(filtered)} rows)")

    # Summary
    print()
    print("Summary:")
    print(f"  Total observations: {len(merged)}")
    print(f"  Unique genes: {merged['gene_id'].nunique()}")
    print(f"  Strains: {merged['strain'].unique().tolist()}")
    print(f"  Columns: {merged.columns.tolist()}")


if __name__ == '__main__':
    main()

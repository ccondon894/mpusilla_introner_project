#!/usr/bin/env python3
"""
Isoform Abundance Analysis
===========================

Counts isoforms per gene per strain from merged SQANTI3 data.

Usage (Snakemake):
    python isoform_abundance_analysis.py \
        --sqanti TSV --output CSV --mode count
"""

import argparse
import pandas as pd
import sys

# Structural categories considered high-confidence
INCLUDE_CATEGORIES = {'full-splice_match', 'novel_in_catalog', 'novel_not_in_catalog'}


def count_isoforms(sqanti_file, output_file):
    """Count unique isoforms per gene per strain."""
    df = pd.read_csv(sqanti_file, sep='\t')

    # Filter by structural category
    if 'structural_category' in df.columns:
        df = df[df['structural_category'].isin(INCLUDE_CATEGORIES)]

    # Count unique isoforms per gene per strain
    counts = df.groupby(['gene_id', 'strain'])['isoform'].nunique().reset_index(
        name='n_isoforms'
    )

    counts.to_csv(output_file, index=False)
    print(f"Counted isoforms for {len(counts)} gene-strain combinations")
    print(f"Output: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Isoform abundance analysis"
    )
    parser.add_argument('--sqanti', required=True,
                        help='Merged SQANTI3 TSV file')
    parser.add_argument('--output', required=True,
                        help='Output CSV file')
    parser.add_argument('--mode', default='count',
                        choices=['count'],
                        help='Analysis mode (default: count)')

    args = parser.parse_args()

    if args.mode == 'count':
        count_isoforms(args.sqanti, args.output)


if __name__ == '__main__':
    main()

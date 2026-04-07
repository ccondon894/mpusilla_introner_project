#!/usr/bin/env python3
"""
Parse SQANTI3 classification files into a merged TSV for downstream analysis.

Reads per-strain SQANTI3 classification files, adds strain identifiers,
cleans gene IDs, and outputs a single merged TSV.

Usage (Snakemake):
    python parse_sqanti3_data.py \
        --ccmp1545 FILE --rcc1614 FILE --rcc1749 FILE \
        --output TSV --summary TXT
"""

import argparse
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Columns to retain from SQANTI3 classification files
KEEP_COLS = [
    'isoform', 'gene_id', 'strain', 'chrom', 'strand', 'length', 'exons',
    'structural_category', 'associated_gene', 'associated_transcript',
    'coding', 'ORF_length', 'CDS_length', 'predicted_NMD'
]


def load_sqanti3_file(filepath, strain):
    """Load a SQANTI3 classification file and add strain column."""
    df = pd.read_csv(filepath, sep='\t')
    df['strain'] = strain

    # Create clean gene_id from associated_gene
    if 'associated_gene' in df.columns:
        df['gene_id'] = df['associated_gene'].fillna(
            df.get('Gene', df.get('isoform', ''))
        )
    elif 'Gene' in df.columns:
        df['gene_id'] = df['Gene']
    else:
        df['gene_id'] = df['isoform']

    # Remove rows with missing gene information
    df = df[df['gene_id'].notna() & (df['gene_id'] != '')]

    return df


def get_replicate_count_cols(df):
    """Identify replicate count columns (e.g., 834_A1, 1614_B1)."""
    import re
    count_cols = [c for c in df.columns if re.match(r'^\d+_[A-Z]\d+$', c)]
    return count_cols


def main():
    parser = argparse.ArgumentParser(
        description="Parse SQANTI3 classification files into merged TSV"
    )
    parser.add_argument('--ccmp1545', required=True,
                        help='SQANTI3 classification file for CCMP1545')
    parser.add_argument('--rcc1614', required=True,
                        help='SQANTI3 classification file for RCC1614')
    parser.add_argument('--rcc1749', required=True,
                        help='SQANTI3 classification file for RCC1749')
    parser.add_argument('--output', required=True,
                        help='Output merged TSV file')
    parser.add_argument('--summary', required=True,
                        help='Output summary text file')

    args = parser.parse_args()

    # Load each strain
    strain_files = {
        'CCMP1545': args.ccmp1545,
        'RCC1614': args.rcc1614,
        'RCC1749': args.rcc1749
    }

    frames = []
    strain_stats = {}

    for strain, filepath in strain_files.items():
        logger.info(f"Loading {strain}: {filepath}")
        df = load_sqanti3_file(filepath, strain)
        logger.info(f"  {strain}: {len(df)} isoforms, "
                     f"{df['gene_id'].nunique()} genes")

        strain_stats[strain] = {
            'n_isoforms': len(df),
            'n_genes': df['gene_id'].nunique()
        }
        frames.append(df)

    # Concatenate all strains
    merged = pd.concat(frames, ignore_index=True)
    logger.info(f"Merged: {len(merged)} isoforms, "
                f"{merged['gene_id'].nunique()} genes")

    # Identify replicate count columns to keep
    count_cols = get_replicate_count_cols(merged)

    # Select columns to output
    output_cols = []
    for col in KEEP_COLS + count_cols:
        if col in merged.columns:
            output_cols.append(col)

    merged_out = merged[output_cols]

    # Write output
    merged_out.to_csv(args.output, sep='\t', index=False)
    logger.info(f"Wrote {len(merged_out)} rows to {args.output}")

    # Write summary
    with open(args.summary, 'w') as f:
        f.write("SQANTI3 Data Parsing Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total isoforms: {len(merged_out)}\n")
        f.write(f"Total unique genes: {merged_out['gene_id'].nunique()}\n\n")

        f.write("Per-strain counts:\n")
        for strain, stats in strain_stats.items():
            f.write(f"  {strain}: {stats['n_isoforms']} isoforms, "
                    f"{stats['n_genes']} genes\n")

        f.write("\nStructural category breakdown:\n")
        if 'structural_category' in merged_out.columns:
            cats = merged_out['structural_category'].value_counts()
            for cat, count in cats.items():
                pct = count / len(merged_out) * 100
                f.write(f"  {cat}: {count} ({pct:.1f}%)\n")

    logger.info(f"Summary written to {args.summary}")


if __name__ == '__main__':
    main()

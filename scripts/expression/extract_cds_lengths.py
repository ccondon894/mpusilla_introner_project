#!/usr/bin/env python3
"""
Extract CDS Lengths from Strain-Specific GTF Files
===================================================

Parses vg_paths.gtf files for CCMP1545, RCC1614, and RCC1749 to extract
the total CDS length for each gene.

CDS length = sum of all CDS feature lengths for a gene
CDS feature length = end - start + 1
"""

import pandas as pd
import re
from collections import defaultdict

print("=" * 80)
print("Extracting CDS Lengths from Strain-Specific GTF Files")
print("=" * 80)
print()

# Define strains and their GTF files
strains = {
    'CCMP1545': '../gene_annotations/CCMP1545.vg_paths.gtf',
    'RCC1614': '../gene_annotations/RCC1614.vg_paths.gtf',
    'RCC1749': '../gene_annotations/RCC1749.vg_paths.gtf'
}

# Store CDS lengths: {strain: {gene_id: total_cds_length}}
cds_lengths = {}

for strain, gtf_path in strains.items():
    print(f"Processing {strain}...")
    print(f"  GTF file: {gtf_path}")

    # Dictionary to accumulate CDS lengths per gene
    gene_cds = defaultdict(int)

    n_cds_features = 0

    with open(gtf_path, 'r') as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            # Check if this is a CDS feature
            if len(fields) >= 9 and fields[2] == 'CDS':
                # Extract gene_id from attributes (column 9)
                attributes = fields[8]

                # Parse gene_id using regex
                match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', attributes)
                if match:
                    gene_id = match.group(1)

                    # Extract start and end positions
                    start = int(fields[3])
                    end = int(fields[4])

                    # Calculate CDS length
                    cds_len = end - start + 1

                    # Accumulate for this gene
                    gene_cds[gene_id] += cds_len
                    n_cds_features += 1

    # Store results for this strain
    cds_lengths[strain] = dict(gene_cds)

    print(f"  CDS features processed: {n_cds_features:,}")
    print(f"  Unique genes found: {len(gene_cds):,}")
    print()

# Convert to DataFrame
print("Converting to DataFrame...")
records = []
for strain in strains.keys():
    for gene_id, cds_len in cds_lengths[strain].items():
        records.append({
            'strain': strain,
            'gene_id': gene_id,
            'cds_length': cds_len
        })

cds_df = pd.DataFrame(records)

print(f"Total records: {len(cds_df):,}")
print()

# Summary statistics
print("Summary statistics by strain:")
print("-" * 80)
for strain in strains.keys():
    strain_data = cds_df[cds_df['strain'] == strain]['cds_length']
    print(f"{strain}:")
    print(f"  Genes: {len(strain_data):,}")
    print(f"  Mean CDS length: {strain_data.mean():.0f} bp")
    print(f"  Median CDS length: {strain_data.median():.0f} bp")
    print(f"  Min CDS length: {strain_data.min():,} bp")
    print(f"  Max CDS length: {strain_data.max():,} bp")
    print()

# Save to CSV
output_file = 'cds_lengths.csv'
cds_df.to_csv(output_file, index=False)
print(f"Saved to: {output_file}")
print()

print("=" * 80)
print("COMPLETE")
print("=" * 80)

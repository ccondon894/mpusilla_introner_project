#!/usr/bin/env python3

import pandas as pd
import sys

def filter_ccmp1545_present_groups():
    """
    Filter the 864 polymorphic groups to only those where CCMP1545 has introner present.
    This ensures we can use CCMP1545 coordinates for VCF extraction.
    """
    print("=== Phase 1: Data Preparation ===")
    print("Step 1: Filtering for CCMP1545-present polymorphic groups")
    
    # Read the polymorphic groups file
    print("Reading polymorphic groups...")
    polymorphic_df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t')
    print(f"Total polymorphic groups: {len(polymorphic_df)}")
    
    # Read the full genotype matrix
    print("Reading genotype matrix...")
    genotype_df = pd.read_csv('/scratch1/chris/introner_vis/genotype_matrixes/genotype_matrix_updated.tsv', sep='\t')
    
    # Filter for CCMP1545 only
    ccmp1545_df = genotype_df[genotype_df['sample'] == 'CCMP1545'].copy()
    print(f"CCMP1545 rows in genotype matrix: {len(ccmp1545_df)}")
    
    # Find polymorphic groups where CCMP1545 has introner present (presence=1)
    ccmp1545_present = ccmp1545_df[ccmp1545_df['presence'] == 1]['ortholog_id'].tolist()
    print(f"Ortholog groups where CCMP1545 has introner present: {len(ccmp1545_present)}")
    
    # Filter polymorphic groups to only those where CCMP1545 is present
    filtered_polymorphic = polymorphic_df[polymorphic_df['ortholog_id'].isin(ccmp1545_present)].copy()
    print(f"Polymorphic groups with CCMP1545 present: {len(filtered_polymorphic)}")
    
    if len(filtered_polymorphic) == 0:
        print("ERROR: No polymorphic groups found where CCMP1545 has introner present!")
        return None
        
    # Get CCMP1545 coordinates for these groups
    ccmp1545_coords = ccmp1545_df[ccmp1545_df['ortholog_id'].isin(filtered_polymorphic['ortholog_id'])].copy()
    
    # Merge with polymorphic group info
    merged_df = filtered_polymorphic.merge(
        ccmp1545_coords[['ortholog_id', 'contig', 'start', 'end', 'gene']], 
        on='ortholog_id', 
        suffixes=('', '_ccmp1545')
    )
    
    print("\nCCMP1545 coordinate summary:")
    print(f"Unique contigs: {merged_df['contig'].nunique()}")
    print(f"Coordinate range: {merged_df['start'].min()} - {merged_df['end'].max()}")
    
    # Show some examples
    print(f"\nFirst 10 groups with CCMP1545 coordinates:")
    print(merged_df[['ortholog_id', 'contig', 'start', 'end', 'gene', 'present_count', 'absent_count']].head(10))
    
    # Save filtered results
    output_file = '/scratch1/chris/introner_vis/polymorphic_groups_ccmp1545_present.tsv'
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nFiltered polymorphic groups saved to: {output_file}")
    
    return merged_df

if __name__ == "__main__":
    result = filter_ccmp1545_present_groups()
    if result is not None:
        print(f"\n✓ Successfully filtered to {len(result)} polymorphic groups where CCMP1545 has introner present")
    else:
        print("\n✗ Failed to filter polymorphic groups")
        sys.exit(1)
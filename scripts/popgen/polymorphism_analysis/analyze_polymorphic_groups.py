#!/usr/bin/env python3

import pandas as pd
import sys

def analyze_polymorphic_groups():
    """
    Analyze genotype matrix to identify polymorphic ortholog groups in Group1
    with no missing data.
    """
    
    # Define Group1 samples (excluding CCMP490, RCC647, RCC835, RCC1749, RCC3052)
    group1_samples = {
        'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 'RCC373', 
        'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833', 'CCMP1545'
    }
    
    print(f"Group1 samples to analyze: {sorted(group1_samples)}")
    print(f"Total Group1 samples: {len(group1_samples)}")
    
    # Read genotype matrix
    print("\nReading genotype matrix...")
    df = pd.read_csv('/scratch1/chris/introner_vis/genotype_matrixes/genotype_matrix_updated.tsv', sep='\t')
    
    print(f"Total rows in genotype matrix: {len(df)}")
    print(f"Unique ortholog groups: {df['ortholog_id'].nunique()}")
    print(f"Unique samples: {df['sample'].nunique()}")
    
    # Filter for Group1 samples only
    group1_df = df[df['sample'].isin(group1_samples)].copy()
    print(f"\nRows for Group1 samples: {len(group1_df)}")
    
    # Analyze polymorphic groups
    polymorphic_results = []
    
    for ortholog_id in group1_df['ortholog_id'].unique():
        ortholog_data = group1_df[group1_df['ortholog_id'] == ortholog_id].copy()
        
        # Check if we have all Group1 samples for this ortholog
        samples_present = set(ortholog_data['sample'])
        if len(samples_present) != len(group1_samples):
            continue  # Skip if not all samples are present
            
        # Check for missing data (presence == 3)
        if (ortholog_data['presence'] == 3).any():
            continue  # Skip if any missing data
            
        # Check if any sample in this group is on CCMP1545 scaffold_2
        ccmp1545_data = ortholog_data[ortholog_data['sample'] == 'CCMP1545']
        if len(ccmp1545_data) > 0 and (ccmp1545_data['contig'] == 'CCMP1545#0#scaffold_2').any():
            continue  # Skip scaffold_2 groups
            
        # Check if polymorphic (both presence=1 and presence=2)
        presence_values = set(ortholog_data['presence'])
        if len(presence_values) > 1 and presence_values.issubset({1, 2}):
            # This is polymorphic
            present_count = (ortholog_data['presence'] == 1).sum()
            absent_count = (ortholog_data['presence'] == 2).sum()
            
            # Get genomic coordinates (use first row as representative)
            first_row = ortholog_data.iloc[0]
            contig = first_row['contig']
            start = first_row['start']
            end = first_row['end']
            gene = first_row['gene']
            
            polymorphic_results.append({
                'ortholog_id': ortholog_id,
                'present_count': present_count,
                'absent_count': absent_count,
                'contig': contig,
                'start': start,
                'end': end,
                'gene': gene,
                'present_samples': sorted([s for s in ortholog_data[ortholog_data['presence'] == 1]['sample']]),
                'absent_samples': sorted([s for s in ortholog_data[ortholog_data['presence'] == 2]['sample']])
            })
    
    print(f"\nPolymorphic ortholog groups found: {len(polymorphic_results)}")
    
    # Display summary
    if polymorphic_results:
        polymorphic_df = pd.DataFrame(polymorphic_results)
        print("\nPolymorphic group summary:")
        print(polymorphic_df[['ortholog_id', 'present_count', 'absent_count', 'contig', 'start', 'end', 'gene']].head(10))
        
        # Save results
        polymorphic_df.to_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t', index=False)
        print(f"\nResults saved to: polymorphic_groups_group1.tsv")
        
        # Distribution of polymorphism
        print("\nPolymorphism distribution:")
        for _, row in polymorphic_df.iterrows():
            print(f"{row['ortholog_id']}: {row['present_count']} present, {row['absent_count']} absent")
            if len(polymorphic_results) <= 10:  # Only show samples for first few
                print(f"  Present: {row['present_samples']}")
                print(f"  Absent: {row['absent_samples']}")
    
    return polymorphic_results

if __name__ == "__main__":
    results = analyze_polymorphic_groups()
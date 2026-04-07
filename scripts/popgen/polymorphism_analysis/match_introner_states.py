#!/usr/bin/env python3

import pandas as pd
import sys
import re

def match_mutations_to_introner_states():
    """
    Match mutations to introner presence/absence status for each sample
    """
    print("=== Matching Mutations to Introner States ===")
    print("Step 5: Linking variants to introner presence/absence per sample")
    
    # Load data
    print("Loading data...")
    
    # Read variant annotations
    annotations_df = pd.read_csv('/scratch1/chris/introner_vis/variant_annotations.tsv', sep='\t')
    print(f"Loaded {len(annotations_df)} variant annotations")
    
    # Read original variants with genotype information
    variants_df = pd.read_csv('/scratch1/chris/introner_vis/introner_variants.tsv', sep='\t')
    print(f"Loaded {len(variants_df)} variants with genotypes")
    
    # Read polymorphic groups with introner states
    polymorphic_df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_group1.tsv', sep='\t')
    print(f"Loaded {len(polymorphic_df)} polymorphic groups")
    
    # Define Group1 samples
    group1_samples = ['CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
                      'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    # Create introner state lookup: ortholog_id -> {sample: presence_status}
    print("Creating introner state lookup...")
    introner_lookup = {}
    
    for ortholog_id in polymorphic_df['ortholog_id'].unique():
        introner_lookup[ortholog_id] = {}
        
        # Get present and absent samples for this ortholog
        group_data = polymorphic_df[polymorphic_df['ortholog_id'] == ortholog_id].iloc[0]
        present_samples = group_data['present_samples']
        absent_samples = group_data['absent_samples'] 
        
        # Parse sample lists (they're stored as strings that look like lists)
        present_list = eval(present_samples) if isinstance(present_samples, str) else present_samples
        absent_list = eval(absent_samples) if isinstance(absent_samples, str) else absent_samples
        
        # Assign presence status
        for sample in present_list:
            if sample in group1_samples:
                introner_lookup[ortholog_id][sample] = 'present'
        
        for sample in absent_list:
            if sample in group1_samples:
                introner_lookup[ortholog_id][sample] = 'absent'
    
    print(f"Created introner lookup for {len(introner_lookup)} ortholog groups")
    
    # Process each variant and extract sample-specific data
    print("Processing variants and genotypes...")
    
    sample_variant_data = []
    
    for idx, var_row in variants_df.iterrows():
        if idx % 500 == 0:
            print(f"Processed {idx}/{len(variants_df)} variants")
            
        ortholog_id = var_row['ortholog_id']
        chrom = var_row['chrom']
        pos = var_row['pos']
        ref = var_row['ref']
        alt = var_row['alt']
        
        # Get annotations for this variant
        var_annotations = annotations_df[
            (annotations_df['ortholog_id'] == ortholog_id) & 
            (annotations_df['pos'] == pos) & 
            (annotations_df['ref'] == ref) & 
            (annotations_df['alt'] == alt)
        ]
        
        # Get the primary annotation (first one, usually most severe)
        if len(var_annotations) > 0:
            primary_annotation = var_annotations.iloc[0]
            mutation_type = primary_annotation['mutation_type']
            annotation = primary_annotation['annotation']
            impact = primary_annotation['impact']
            gene_name = primary_annotation['gene_name']
        else:
            mutation_type = 'unknown'
            annotation = 'unknown'
            impact = 'unknown'
            gene_name = ''
        
        # Process each sample's genotype
        for sample in group1_samples:
            gt_col = f'{sample}_GT'
            if gt_col in var_row and pd.notna(var_row[gt_col]):
                genotype = var_row[gt_col]
                
                # Parse genotype (format: GT:AD:DP:GQ:PL)
                gt_fields = str(genotype).split(':')
                if len(gt_fields) > 0:
                    gt = gt_fields[0]  # Genotype field
                    
                    # Determine if sample has variant
                    has_variant = determine_variant_status(gt)
                    
                    # Get introner status for this sample and ortholog
                    introner_status = introner_lookup.get(ortholog_id, {}).get(sample, 'unknown')
                    
                    # Only include if we know the introner status
                    if introner_status in ['present', 'absent']:
                        sample_variant_data.append({
                            'ortholog_id': ortholog_id,
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'sample': sample,
                            'genotype': gt,
                            'has_variant': has_variant,
                            'introner_status': introner_status,
                            'mutation_type': mutation_type,
                            'annotation': annotation,
                            'impact': impact,
                            'gene_name': gene_name
                        })
    
    # Convert to DataFrame
    final_df = pd.DataFrame(sample_variant_data)
    
    print(f"\nFinal dataset summary:")
    print(f"Total sample-variant records: {len(final_df)}")
    print(f"Unique variants: {final_df[['chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]}")
    print(f"Unique orthologs: {final_df['ortholog_id'].nunique()}")
    print(f"Samples: {final_df['sample'].nunique()}")
    
    print(f"\nIntroner status distribution:")
    print(final_df['introner_status'].value_counts())
    
    print(f"\nMutation type distribution:")
    print(final_df['mutation_type'].value_counts())
    
    print(f"\nVariant presence distribution:")
    print(final_df['has_variant'].value_counts())
    
    # Cross-tabulation
    print(f"\nCross-tabulation: Introner status vs Mutation type")
    crosstab = pd.crosstab(final_df['introner_status'], final_df['mutation_type'], margins=True)
    print(crosstab)
    
    # Save final dataset
    output_file = '/scratch1/chris/introner_vis/final_selection_dataset.tsv'
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nFinal dataset saved to: {output_file}")
    
    return final_df

def determine_variant_status(genotype):
    """
    Determine if sample has the variant based on genotype
    """
    if genotype in ['./.', '.', '0/0', '0']:
        return False  # No variant (reference)
    elif genotype in ['0/1', '1/0', '0|1', '1|0']:
        return True   # Heterozygous variant
    elif genotype in ['1/1', '1']:
        return True   # Homozygous variant
    else:
        return False  # Unknown/missing, treat as no variant

if __name__ == "__main__":
    result = match_mutations_to_introner_states()
    if result is not None:
        print(f"\n✓ Successfully created selection analysis dataset with {len(result)} records")
    else:
        print("\n✗ Failed to create dataset")
        sys.exit(1)
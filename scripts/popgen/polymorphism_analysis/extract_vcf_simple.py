#!/usr/bin/env python3

import pandas as pd
import sys
import gzip

def extract_vcf_regions_simple():
    """
    Simple VCF extraction without bcftools - parse VCF directly
    CRITICAL FIX: Only extracts flanking regions (200bp total) to avoid coordinate bias
    - Left flanking: start to start+100bp
    - Right flanking: end-100bp to end
    - Excludes introner sequence: start+100bp to end-100bp
    """
    print("=== Simple VCF Region Extraction (Flanking Only) ===")
    
    # Read the regions we want
    print("Reading target regions...")
    regions_df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_ccmp1545_present.tsv', sep='\t')
    
    # Create a lookup for regions: contig -> list of (start, end, ortholog_id)
    region_lookup = {}
    for _, row in regions_df.iterrows():
        contig = row['contig']
        start = row['start']
        end = row['end'] 
        ortholog_id = row['ortholog_id']
        
        if contig not in region_lookup:
            region_lookup[contig] = []
        region_lookup[contig].append((start, end, ortholog_id))
    
    print(f"Target regions: {len(regions_df)} across {len(region_lookup)} contigs")
    
    # Define Group1 samples (excluding unwanted ones)
    group1_samples = ['CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
                      'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    # Read VCF file
    vcf_file = '/scratch1/chris/introner_vis/data_prep/snpEff_run2/mpusilla.snps.snpEff.no_MT.vcf'
    output_file = '/scratch1/chris/introner_vis/introner_variants.tsv'
    
    print(f"Processing VCF: {vcf_file}")
    
    # Parse VCF header to get sample column positions
    sample_indices = {}
    header_line = None
    
    variants_found = []
    
    try:
        with open(vcf_file, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                
                # Progress indicator
                if line_count % 100000 == 0:
                    print(f"Processed {line_count} lines, found {len(variants_found)} variants")
                
                if line.startswith('##'):
                    continue  # Skip metadata
                    
                if line.startswith('#CHROM'):
                    # Parse header to get sample positions
                    header_line = line.strip().split('\t')
                    for i, sample in enumerate(header_line):
                        if sample in group1_samples:
                            sample_indices[sample] = i
                    print(f"Found {len(sample_indices)} Group1 samples in VCF header")
                    continue
                
                # Process variant line
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                        
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    info = fields[7]
                    
                    # Check if this variant is in any of our target flanking regions
                    # Only analyze flanking regions: start to start+100bp and end-100bp to end
                    if chrom in region_lookup:
                        for start, end, ortholog_id in region_lookup[chrom]:
                            # Left flanking region: start to start+100bp
                            # Right flanking region: end-100bp to end  
                            # Exclude introner sequence: start+100bp to end-100bp
                            in_left_flanking = start <= pos <= (start + 100)
                            in_right_flanking = (end - 100) <= pos <= end
                            if in_left_flanking or in_right_flanking:
                                # This variant is in our target region
                                variant_data = {
                                    'ortholog_id': ortholog_id,
                                    'chrom': chrom,
                                    'pos': pos,
                                    'ref': ref,
                                    'alt': alt,
                                    'info': info
                                }
                                
                                # Extract genotypes for Group1 samples
                                for sample in group1_samples:
                                    if sample in sample_indices:
                                        gt_index = sample_indices[sample]
                                        if gt_index < len(fields):
                                            variant_data[f'{sample}_GT'] = fields[gt_index]
                                
                                variants_found.append(variant_data)
                                break  # Found in one region, don't need to check others
                    
    except Exception as e:
        print(f"ERROR reading VCF: {e}")
        return None
    
    print(f"\nFound {len(variants_found)} variants in target regions")
    
    if len(variants_found) > 0:
        # Convert to DataFrame and save
        variants_df = pd.DataFrame(variants_found)
        variants_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Variants saved to: {output_file}")
        print("\nVariant summary:")
        print(f"Unique orthologs with variants: {variants_df['ortholog_id'].nunique()}")
        print(f"Chromosomes: {variants_df['chrom'].nunique()}")
        print(f"Position range: {variants_df['pos'].min()} - {variants_df['pos'].max()}")
        
        # Show first few variants
        print("\nFirst 5 variants:")
        print(variants_df[['ortholog_id', 'chrom', 'pos', 'ref', 'alt']].head())
        
        return output_file
    else:
        print("No variants found in target regions")
        return None

if __name__ == "__main__":
    result = extract_vcf_regions_simple()
    if result:
        print(f"\n✓ VCF extraction successful: {result}")
    else:
        print("\n✗ VCF extraction failed")
        sys.exit(1)
#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import os

'''
To run:
python scripts/calculate_haplotype_diversity_ratios.py --introner_metrics group1_evolution_analysis_no/diversity_metrics/group1_diversity_metrics_200bp.tsv --genotype_matrix genotype_matrixes/genotype_matrix_updated.tsv
  --fourfold_lookup group1_4fold_analysis/diversity_metrics/4fold_avg_pi_present_samples.tsv --output
  group1_evolution_analysis_no_consensus/diversity_metrics/haplotype_diversity_ratios.tsv --min_sites 20
  
Note: Uses a single 4-fold lookup table for both PHDR and AHDR calculations, providing proper sample-set-specific normalization.
Ratios are only calculated when the 4D pi estimate is based on at least --min_sites sites (default: 20).
'''

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate Present and Absent Haplotype Diversity Ratios (PHDR and AHDR) for polymorphic introners')
    parser.add_argument('--introner_metrics', required=True, help='Input TSV file with introner diversity metrics')
    parser.add_argument('--genotype_matrix', required=True, help='Input TSV file with genotype matrix')
    parser.add_argument('--fourfold_lookup', required=True, help='4-fold diversity lookup table (sample_set -> pi_4D)')
    parser.add_argument('--output', required=True, help='Output TSV file for haplotype diversity ratios')
    parser.add_argument('--min_sites', type=int, default=20, help='Minimum number of 4D sites required for ratio calculation (default: 20)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    return parser.parse_args()

def format_sample_set(samples):
    """Convert list of samples to sorted set format with curly braces"""
    if not samples:
        return ''
    sorted_samples = sorted(samples)
    return '{' + ','.join(sorted_samples) + '}'

def load_4fold_lookup_table(fourfold_file, verbose=False):
    """Load 4-fold diversity lookup table with pi_4D and site_count information"""
    if verbose:
        print(f"Loading 4-fold lookup table...")
    
    # Load 4-fold lookup table
    fourfold_df = pd.read_csv(fourfold_file, sep='\t')
    
    # Create dictionaries for both pi_4D and site_count
    pi_lookup = dict(zip(fourfold_df['sample_set'], fourfold_df['pi_4D']))
    site_count_lookup = dict(zip(fourfold_df['sample_set'], fourfold_df['site_count']))
    
    if verbose:
        print(f"  Loaded {len(pi_lookup)} sample sets with 4-fold diversity values")
    
    return pi_lookup, site_count_lookup

def get_sample_presence_for_ortholog(genotype_df, ortholog_id, verbose=False):
    """Get present/absent sample lists for a specific ortholog_id"""
    # Filter for this ortholog_id
    ortholog_data = genotype_df[genotype_df['ortholog_id'] == ortholog_id].copy()
    
    if len(ortholog_data) == 0:
        if verbose:
            print(f"  Warning: No data found for {ortholog_id}")
        return [], []
    
    # Exclude Group2 samples (RCC1749, RCC3052)
    group2_samples = {'RCC1749', 'RCC3052'}
    ortholog_data = ortholog_data[~ortholog_data['sample'].isin(group2_samples)]
    
    # Classify samples by presence
    present_samples = ortholog_data[ortholog_data['presence'] == 1]['sample'].tolist()
    absent_samples = ortholog_data[ortholog_data['presence'] == 2]['sample'].tolist()
    
    if verbose:
        print(f"  {ortholog_id}: {len(present_samples)} present, {len(absent_samples)} absent")
    
    return present_samples, absent_samples

def calculate_haplotype_diversity_ratios(introner_metrics_file, genotype_matrix_file, 
                                       fourfold_lookup_file, output_file, min_sites=20, verbose=False):
    """Main function to calculate haplotype diversity ratios"""
    
    # Load input data
    if verbose:
        print("Loading input data...")
    
    introner_df = pd.read_csv(introner_metrics_file, sep='\t')
    genotype_df = pd.read_csv(genotype_matrix_file, sep='\t')
    fourfold_pi_lookup, fourfold_site_lookup = load_4fold_lookup_table(fourfold_lookup_file, verbose)

    # Filter to only truly polymorphic loci: must have both presence==1 AND
    # presence==2 calls within Group 1. Loci where "absent" samples are all
    # presence==3 (not callable) are excluded.
    group2_samples = {'RCC1749', 'RCC3052'}
    g1_gm = genotype_df[~genotype_df['sample'].isin(group2_samples)]

    polymorphic_oids = set()
    for oid, grp in g1_gm.groupby('ortholog_id'):
        presences = set(grp['presence'])
        if 1 in presences and 2 in presences:
            polymorphic_oids.add(oid)

    n_before = len(introner_df)
    introner_df = introner_df[introner_df['ortholog_id'].isin(polymorphic_oids)]
    n_after = len(introner_df)

    if verbose:
        print(f"Filtered to truly polymorphic loci (both present and absent calls in Group 1):")
        print(f"  Before: {n_before}, After: {n_after}, Removed: {n_before - n_after}")

    results = []
    missing_fourfold_count = 0
    insufficient_sites_count = 0
    singleton_excluded_count = 0

    introner_ortholog_ids = set(introner_df['ortholog_id'])

    if verbose:
        print(f"Processing {len(introner_ortholog_ids)} polymorphic introners...")
    
    # Process each introner from the diversity metrics
    for idx, row in introner_df.iterrows():
        ortholog_id = row['ortholog_id']
        pi_present_introner = row['pi_present']
        pi_absent_introner = row['pi_absent']
        
        if verbose and (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(introner_df)} introners...")
        
        # Get sample presence for this ortholog
        present_samples, absent_samples = get_sample_presence_for_ortholog(genotype_df, ortholog_id, verbose)
        
        # Check for singleton sample sets (pi cannot be calculated with <2 samples)
        can_calculate_phdr = len(present_samples) >= 2
        can_calculate_ahdr = len(absent_samples) >= 2
        
        # Skip introners where no meaningful pi calculations can be made
        if not can_calculate_phdr and not can_calculate_ahdr:
            singleton_excluded_count += 1
            if verbose:
                print(f"  Skipping {ortholog_id}: Both present ({len(present_samples)}) and absent ({len(absent_samples)}) are singletons")
            continue
        
        # Format sample sets for lookup
        present_sample_set = format_sample_set(present_samples)
        absent_sample_set = format_sample_set(absent_samples)
        
        # Look up 4-fold pi values and site counts using the same lookup table for both
        avg_4D_pi_present = fourfold_pi_lookup.get(present_sample_set, None) if can_calculate_phdr else None
        avg_4D_pi_absent = fourfold_pi_lookup.get(absent_sample_set, None) if can_calculate_ahdr else None
        sites_4D_present = fourfold_site_lookup.get(present_sample_set, 0) if can_calculate_phdr else 0
        sites_4D_absent = fourfold_site_lookup.get(absent_sample_set, 0) if can_calculate_ahdr else 0
        
        # Track missing lookups (only for sample sets that can have pi calculated)
        if can_calculate_phdr and present_sample_set and avg_4D_pi_present is None:
            missing_fourfold_count += 1
            if verbose:
                print(f"  Warning: No 4-fold data for present set {present_sample_set}")
        
        if can_calculate_ahdr and absent_sample_set and avg_4D_pi_absent is None:
            missing_fourfold_count += 1
            if verbose:
                print(f"  Warning: No 4-fold data for absent set {absent_sample_set}")
        
        # Check minimum site thresholds
        phdr_sufficient_sites = sites_4D_present >= min_sites if can_calculate_phdr else False
        ahdr_sufficient_sites = sites_4D_absent >= min_sites if can_calculate_ahdr else False
        
        if can_calculate_phdr and avg_4D_pi_present is not None and not phdr_sufficient_sites:
            insufficient_sites_count += 1
            if verbose:
                print(f"  Warning: Insufficient 4D sites for present set {present_sample_set}: {sites_4D_present} < {min_sites}")
        
        if can_calculate_ahdr and avg_4D_pi_absent is not None and not ahdr_sufficient_sites:
            insufficient_sites_count += 1
            if verbose:
                print(f"  Warning: Insufficient 4D sites for absent set {absent_sample_set}: {sites_4D_absent} < {min_sites}")
        
        # Calculate ratios (only for non-singleton sample sets)
        phdr = None
        ahdr = None
        
        if (can_calculate_phdr and not pd.isna(pi_present_introner) and 
            avg_4D_pi_present is not None and avg_4D_pi_present > 0 and phdr_sufficient_sites):
            phdr = pi_present_introner / avg_4D_pi_present
        elif can_calculate_phdr and verbose:
            if not phdr_sufficient_sites and avg_4D_pi_present is not None:
                print(f"  Skipping PHDR for {ortholog_id}: insufficient 4D sites ({sites_4D_present} < {min_sites})")
            elif not avg_4D_pi_present:
                print(f"  Skipping PHDR for {ortholog_id}: no 4D data available")
        elif not can_calculate_phdr and verbose:
            print(f"  Skipping PHDR for {ortholog_id}: present samples ({len(present_samples)}) insufficient for pi calculation")
        
        if (can_calculate_ahdr and not pd.isna(pi_absent_introner) and 
            avg_4D_pi_absent is not None and avg_4D_pi_absent > 0 and ahdr_sufficient_sites):
            ahdr = pi_absent_introner / avg_4D_pi_absent
        elif can_calculate_ahdr and verbose:
            if not ahdr_sufficient_sites and avg_4D_pi_absent is not None:
                print(f"  Skipping AHDR for {ortholog_id}: insufficient 4D sites ({sites_4D_absent} < {min_sites})")
            elif not avg_4D_pi_absent:
                print(f"  Skipping AHDR for {ortholog_id}: no 4D data available")
        elif not can_calculate_ahdr and verbose:
            print(f"  Skipping AHDR for {ortholog_id}: absent samples ({len(absent_samples)}) insufficient for pi calculation")
        
        # Store results (only for introners in diversity metrics)
        results.append({
            'ortholog_id': ortholog_id,
            'present_samples': ','.join(present_samples),
            'absent_samples': ','.join(absent_samples),
            'pi_present_introner': pi_present_introner,
            'pi_absent_introner': pi_absent_introner,
            'avg_4D_pi_present': avg_4D_pi_present,
            'avg_4D_pi_absent': avg_4D_pi_absent,
            'sites_4D_present': sites_4D_present,
            'sites_4D_absent': sites_4D_absent,
            'PHDR': phdr,
            'AHDR': ahdr
        })
    
    # Create output DataFrame and save
    results_df = pd.DataFrame(results)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save results
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary statistics
    if verbose:
        print(f"\nSummary:")
        print(f"  Total introners in input: {len(introner_df)}")
        print(f"  Introners excluded (both present and absent singletons): {singleton_excluded_count}")
        print(f"  Total introners included in output: {len(results)}")
        print(f"  Missing sample sets in 4-fold lookup: {missing_fourfold_count}")
        print(f"  Ratios skipped due to insufficient 4D sites (< {min_sites}): {insufficient_sites_count}")
        
        # Calculate ratio statistics
        phdr_values = [r['PHDR'] for r in results if r['PHDR'] is not None]
        ahdr_values = [r['AHDR'] for r in results if r['AHDR'] is not None]
        
        print(f"  PHDR calculated for: {len(phdr_values)} introners")
        if phdr_values:
            print(f"    PHDR range: {min(phdr_values):.4f} - {max(phdr_values):.4f}")
            print(f"    PHDR median: {np.median(phdr_values):.4f}")
        
        print(f"  AHDR calculated for: {len(ahdr_values)} introners")
        if ahdr_values:
            print(f"    AHDR range: {min(ahdr_values):.4f} - {max(ahdr_values):.4f}")
            print(f"    AHDR median: {np.median(ahdr_values):.4f}")
    
    print(f"Saved haplotype diversity ratios to {output_file}")
    print("Done!")

def main():
    args = parse_arguments()
    
    calculate_haplotype_diversity_ratios(
        args.introner_metrics,
        args.genotype_matrix,
        args.fourfold_lookup,
        args.output,
        args.min_sites,
        args.verbose
    )

if __name__ == "__main__":
    main()
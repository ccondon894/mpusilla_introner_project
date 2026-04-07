#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import gzip
from collections import defaultdict, Counter

def calculate_introner_sfs(genotype_file, group1_samples):
    """Calculate site frequency spectrum for introners"""
    
    df = pd.read_csv(genotype_file, sep='\t')
    
    # Filter to Group1 samples only
    df_group1 = df[df['sample'].isin(group1_samples)].copy()
    
    # Find ortholog groups with missing data (presence == 3)
    missing_data_groups = df_group1[df_group1['presence'] == 3]['ortholog_id'].unique()
    print(f"Excluding {len(missing_data_groups)} ortholog groups with missing data")
    
    # Find ortholog groups on CCMP1545 scaffold_2 contig
    scaffold2_groups = df_group1[(df_group1['sample'] == 'CCMP1545') & 
                                (df_group1['contig'] == 'CCMP1545#0#scaffold_2')]['ortholog_id'].unique()
    print(f"Excluding {len(scaffold2_groups)} ortholog groups on CCMP1545 scaffold_2")
    
    # Combine exclusion sets
    excluded_groups = set(missing_data_groups) | set(scaffold2_groups)
    print(f"Total excluded groups: {len(excluded_groups)}")
    
    # Exclude ortholog groups with any missing data or on scaffold_2
    df_complete = df_group1[~df_group1['ortholog_id'].isin(excluded_groups)]
    
    # Count presence (1) for each ortholog in complete data only
    introner_counts = df_complete[df_complete['presence'] == 1].groupby('ortholog_id').size()
    
    # Create frequency spectrum (count how many introners appear in 1, 2, 3, ... samples)
    freq_counts = Counter(introner_counts.values)
    
    # Fill in missing frequencies with 0, excluding the fixed frequency (max_samples)
    max_samples = len(group1_samples)
    sfs = [freq_counts.get(i, 0) for i in range(1, max_samples)]
    
    print(f"Used {len(introner_counts)} complete ortholog groups for SFS calculation")
    
    return sfs

def parse_snpeff_vcf(vcf_file, group1_samples, annotation_type):
    """Extract mutations of specific type from SnpEff VCF and calculate SFS"""
    
    # Map annotation types to SnpEff impact categories
    impact_map = {
        'synonymous': ['synonymous_variant'],
        'nonsynonymous': ['missense_variant', 'stop_gained', 'stop_lost', 'start_lost']
    }
    
    target_impacts = impact_map[annotation_type]
    mutation_counts = Counter()
    
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(vcf_file, mode) as f:
        header_samples = None
        group1_indices = None
        
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                # Parse sample names from header
                fields = line.strip().split('\t')
                header_samples = fields[9:]  # Sample columns start at index 9
                
                # Find indices of Group1 samples
                group1_indices = []
                for sample in group1_samples:
                    if sample in header_samples:
                        group1_indices.append(header_samples.index(sample) + 9)  # +9 for VCF offset
                
                print(f"Found {len(group1_indices)} Group1 samples in VCF for {annotation_type} analysis")
                continue
            
            # Process variant lines
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            info_field = fields[7]
            
            # Check if this variant has the target annotation
            has_target_annotation = False
            for impact in target_impacts:
                if impact in info_field:
                    has_target_annotation = True
                    break
            
            if not has_target_annotation:
                continue
            
            # Count alleles in Group1 samples
            alt_count = 0
            total_count = 0
            
            for idx in group1_indices:
                if idx < len(fields):
                    genotype = fields[idx].split(':')[0]
                    if genotype not in ['./.', '.']:
                        alleles = genotype.replace('|', '/').split('/')
                        for allele in alleles:
                            if allele != '.' and allele != '0':  # Non-reference allele
                                alt_count += 1
                            total_count += 1
            
            if alt_count > 0:
                mutation_counts[alt_count] += 1
    
    # Convert to SFS format (frequencies 1 to max_samples-1, excluding fixed)
    max_samples = len(group1_samples)
    sfs = [mutation_counts.get(i, 0) for i in range(1, max_samples)]
    
    return sfs

def plot_comparative_sfs(introner_sfs, syn_sfs, nonsyn_sfs, output_file):
    """Create comparative site frequency spectrum plot"""
    
    frequencies = list(range(1, len(introner_sfs) + 1))
    
    plt.figure(figsize=(12, 8))
    
    # Convert zeros to small values for log scale (log of 0 is undefined)
    introner_sfs_log = [max(x, 0.1) for x in introner_sfs]
    syn_sfs_log = [max(x, 0.1) for x in syn_sfs]
    nonsyn_sfs_log = [max(x, 0.1) for x in nonsyn_sfs]
    
    # Plot all three spectra
    plt.bar([f - 0.25 for f in frequencies], introner_sfs_log, width=0.25, 
            label='Introners', alpha=0.8, color='blue')
    plt.bar(frequencies, syn_sfs_log, width=0.25, 
            label='Synonymous mutations', alpha=0.8, color='green')
    plt.bar([f + 0.25 for f in frequencies], nonsyn_sfs_log, width=0.25, 
            label='Nonsynonymous mutations', alpha=0.8, color='red')
    
    plt.xlabel('Allele frequency (number of samples)')
    plt.ylabel('Number of variants/introners (log scale)')
    plt.title('Site Frequency Spectrum Comparison:\nIntroners vs. Genome-wide Mutations')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Set log scale for y-axis
    plt.yscale('log')
    
    # Set x-axis to show all frequencies
    plt.xticks(frequencies)
    plt.xlim(0.5, len(frequencies) + 0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

def main():
    # Define Group1 samples (excluding Group2: RCC1749, RCC3052)
    group1_samples = [
        'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
        'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
    ]
    
    print("Calculating comparative site frequency spectra...")
    print(f"Group1 samples: {group1_samples}")
    print(f"Total Group1 samples: {len(group1_samples)}")
    
    # File paths
    genotype_file = 'genotype_matrixes/genotype_matrix_updated.tsv'
    vcf_file = 'data_prep/snpEff_run2/mpusilla.snps.snpEff.no_MT.vcf' # exclude mating type chromosome
    
    print("\n1. Calculating introner SFS...")
    introner_sfs = calculate_introner_sfs(genotype_file, group1_samples)
    print(f"Introner SFS: {introner_sfs}")
    
    print("\n2. Calculating synonymous mutation SFS...")
    syn_sfs = parse_snpeff_vcf(vcf_file, group1_samples, 'synonymous')
    print(f"Synonymous SFS: {syn_sfs}")
    
    print("\n3. Calculating nonsynonymous mutation SFS...")
    nonsyn_sfs = parse_snpeff_vcf(vcf_file, group1_samples, 'nonsynonymous')
    print(f"Nonsynonymous SFS: {nonsyn_sfs}")
    
    print("\n4. Creating comparative plot...")
    plot_comparative_sfs(introner_sfs, syn_sfs, nonsyn_sfs, 'comparative_sfs_analysis.png')
    
    # Summary statistics
    print("\nSummary Statistics:")
    print(f"Total introners: {sum(introner_sfs)}")
    print(f"Total synonymous mutations: {sum(syn_sfs)}")
    print(f"Total nonsynonymous mutations: {sum(nonsyn_sfs)}")
    
    # Save results to file
    results_df = pd.DataFrame({
        'allele_frequency': list(range(1, len(introner_sfs) + 1)),
        'introners': introner_sfs,
        'synonymous_mutations': syn_sfs,
        'nonsynonymous_mutations': nonsyn_sfs
    })
    
    results_df.to_csv('comparative_sfs_results.tsv', sep='\t', index=False)
    print("\nResults saved to comparative_sfs_results.tsv")

if __name__ == "__main__":
    main()
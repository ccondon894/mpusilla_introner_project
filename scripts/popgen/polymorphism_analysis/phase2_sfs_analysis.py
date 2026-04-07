#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import Counter
import sys
import argparse
import pysam

def fold_frequency(freq, total_samples):
    """
    Fold allele frequency to remove ancestral/derived polarity.
    For a frequency f in n samples, the folded frequency is min(f, n-f).
    """
    return min(freq, total_samples - freq)

def calculate_unfolded_allele_frequencies(vcf_file):
    """
    Calculate unfolded allele frequencies using Group2 (RCC1749, RCC3052) as outgroup for polarization
    Assumes haploid genome data.
    """
    print("=== Phase 2: Site Frequency Spectrum Analysis (Unfolded) ===")
    
    # Sample definitions
    outgroup_samples = ['RCC1749', 'RCC3052']  # Group2 - outgroup for ancestral state inference
    group1_samples = [
        'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
        'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
    ]
    obsolete_samples = {'CCMP490', 'RCC647', 'RCC835'}
    
    print(f"Outgroup samples (Group2): {outgroup_samples}")
    print(f"Group1 samples: {group1_samples}")
    
    # Process variants directly from SnpEff VCF
    
    # Open VCF for genotype data
    print(f"Opening VCF file: {vcf_file}")
    try:
        vcf = pysam.VariantFile(vcf_file)
        sample_names = list(vcf.header.samples)
        print(f"VCF opened successfully with {len(sample_names)} samples")
    except Exception as e:
        print(f"Error opening VCF: {e}")
        return pd.DataFrame()
    
    # Get sample indices
    outgroup_indices = [sample_names.index(sample) for sample in outgroup_samples if sample in sample_names]
    group1_indices = [sample_names.index(sample) for sample in group1_samples if sample in sample_names]
    
    print(f"Found {len(outgroup_indices)} outgroup samples in VCF")
    print(f"Found {len(group1_indices)} Group1 samples in VCF")
    
    unfolded_freq_data = []
    sites_processed = 0
    sites_excluded_polymorphic_outgroup = 0
    sites_excluded_missing_outgroup = 0
    sites_excluded_monomorphic = 0
    sites_excluded_missing_group1 = 0
    sites_used = 0
    
    # Process variants directly from VCF
    print("Processing variants from VCF...")
    
    # Process each variant record
    for record in vcf:
        if record.alts and len(record.alts) == 1:  # Biallelic SNP only
            sites_processed += 1
            
            if sites_processed % 10000 == 0:
                print(f"Processed {sites_processed} sites...")
            
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]
            
            # Get outgroup genotypes for polarization
            outgroup_genotypes = []
            outgroup_missing = False
            
            for i in outgroup_indices:
                gt = record.samples[i]['GT']
                if None in gt:  # Missing data
                    outgroup_missing = True
                    break
                outgroup_genotypes.extend([allele for allele in gt if allele is not None])
            
            # Skip if missing data in outgroup
            if outgroup_missing:
                sites_excluded_missing_outgroup += 1
                continue
            
            # Check if outgroup is monomorphic
            unique_alleles = set(outgroup_genotypes)
            if len(unique_alleles) > 1:
                sites_excluded_polymorphic_outgroup += 1
                continue
            
            # Determine ancestral allele (consensus in outgroup)
            ancestral_allele = list(unique_alleles)[0]
            
            # Count derived alleles in Group1
            derived_count = 0
            group1_missing = False
            
            for i in group1_indices:
                gt = record.samples[i]['GT']
                if None in gt:  # Missing data
                    group1_missing = True
                    break
                
                # Count alleles different from ancestral
                for allele in gt:
                    if allele != ancestral_allele:
                        derived_count += 1
            
            # Skip if missing data in Group1
            if group1_missing:
                sites_excluded_missing_group1 += 1
                continue
            
            # Skip monomorphic sites (fixed in Group1, haploid)
            n_group1_chromosomes = len(group1_indices)
            if derived_count == 0 or derived_count == n_group1_chromosomes:
                sites_excluded_monomorphic += 1
                continue
            
            # Calculate derived allele frequency
            derived_frequency = derived_count / n_group1_chromosomes
            
            # Extract SnpEff annotations from INFO field
            mutation_type = 'unknown'
            annotation = 'unknown'
            impact = 'unknown'
            
            if 'ANN' in record.info:
                # SnpEff annotation format: Allele|Annotation|Impact|Gene_Name|...
                ann_data = record.info['ANN']
                if isinstance(ann_data, (list, tuple)):
                    ann_string = ann_data[0]
                else:
                    ann_string = ann_data
                ann_fields = ann_string.split('|')
                if len(ann_fields) >= 3:
                    annotation = ann_fields[1]
                    impact = ann_fields[2]
                    
                    # Map SnpEff annotations to mutation types
                    if 'synonymous_variant' in annotation:
                        mutation_type = 'synonymous'
                    elif any(x in annotation for x in ['missense_variant', 'stop_gained', 'stop_lost']):
                        mutation_type = 'nonsynonymous'
                        if 'stop_gained' in annotation or 'stop_lost' in annotation:
                            mutation_type = 'nonsense'
                        else:
                            mutation_type = 'missense'
            
            # Add to unfolded frequency data
            unfolded_freq_data.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'derived_count': derived_count,
                'total_chromosomes': n_group1_chromosomes,
                'derived_frequency': derived_frequency,
                'mutation_type': mutation_type,
                'annotation': annotation,
                'impact': impact,
                'ancestral_allele': ancestral_allele
            })
            
            sites_used += 1
    
    print(f"\nSite filtering summary:")
    print(f"Total variants processed: {sites_processed}")
    print(f"Sites excluded (polymorphic outgroup): {sites_excluded_polymorphic_outgroup}")
    print(f"Sites excluded (missing data in outgroup): {sites_excluded_missing_outgroup}")
    print(f"Sites excluded (missing data in Group1): {sites_excluded_missing_group1}")
    print(f"Sites excluded (monomorphic in Group1): {sites_excluded_monomorphic}")
    print(f"Sites used for unfolded AFS: {sites_used}")
    
    unfolded_freq_df = pd.DataFrame(unfolded_freq_data)
    print(f"Calculated unfolded allele frequencies for {len(unfolded_freq_df)} variant-introner state combinations")
    
    return unfolded_freq_df

def calculate_allele_frequencies():
    """
    Calculate allele frequencies for variants in the dataset with folded frequencies
    """
    print("=== Phase 2: Site Frequency Spectrum Analysis (Folded) ===")
    
    # Load the final dataset
    print("Loading selection dataset...")
    df = pd.read_csv('/scratch1/chris/introner_vis/polymorphism_analysis/final_selection_dataset.tsv', sep='\t')
    print(f"Loaded {len(df)} sample-variant records")
    
    # Calculate allele frequencies for each variant position
    print("Calculating folded allele frequencies...")
    
    allele_freq_data = []
    
    # Group by variant position (chrom, pos, ref, alt) and ortholog
    for (ortholog_id, chrom, pos, ref, alt), variant_group in df.groupby(['ortholog_id', 'chrom', 'pos', 'ref', 'alt']):
        
        # Count samples with and without variant for each introner state
        present_samples = variant_group[variant_group['introner_status'] == 'present']
        absent_samples = variant_group[variant_group['introner_status'] == 'absent']
        
        # Calculate frequencies for introner-present samples
        if len(present_samples) > 0:
            present_with_variant = len(present_samples[present_samples['has_variant'] == True])
            present_total = len(present_samples)
            present_freq = present_with_variant / present_total
            present_folded_count = fold_frequency(present_with_variant, present_total)
            present_folded_freq = present_folded_count / present_total
            
            # Get mutation info
            mutation_type = present_samples['mutation_type'].iloc[0]
            annotation = present_samples['annotation'].iloc[0]
            impact = present_samples['impact'].iloc[0]
            
            allele_freq_data.append({
                'ortholog_id': ortholog_id,
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'introner_status': 'present',
                'with_variant': present_with_variant,
                'total_samples': present_total,
                'allele_frequency': present_freq,
                'folded_count': present_folded_count,
                'folded_frequency': present_folded_freq,
                'mutation_type': mutation_type,
                'annotation': annotation,
                'impact': impact
            })
        
        # Calculate frequencies for introner-absent samples
        if len(absent_samples) > 0:
            absent_with_variant = len(absent_samples[absent_samples['has_variant'] == True])
            absent_total = len(absent_samples)
            absent_freq = absent_with_variant / absent_total
            absent_folded_count = fold_frequency(absent_with_variant, absent_total)
            absent_folded_freq = absent_folded_count / absent_total
            
            # Get mutation info
            mutation_type = absent_samples['mutation_type'].iloc[0]
            annotation = absent_samples['annotation'].iloc[0]
            impact = absent_samples['impact'].iloc[0]
            
            allele_freq_data.append({
                'ortholog_id': ortholog_id,
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'introner_status': 'absent',
                'with_variant': absent_with_variant,
                'total_samples': absent_total,
                'allele_frequency': absent_freq,
                'folded_count': absent_folded_count,
                'folded_frequency': absent_folded_freq,
                'mutation_type': mutation_type,
                'annotation': annotation,
                'impact': impact
            })
    
    allele_freq_df = pd.DataFrame(allele_freq_data)
    print(f"Calculated folded allele frequencies for {len(allele_freq_df)} variant-introner state combinations")
    
    return allele_freq_df

def calculate_unfolded_introner_frequencies():
    """
    Calculate unfolded introner frequencies assuming absent state is ancestral
    """
    print("\n=== Calculating Unfolded Introner Frequencies ===")
    
    # Import required modules for introner analysis
    import pandas as pd
    from collections import Counter
    
    # Define Group1 samples (from the comparative analysis)
    group1_samples = [
        'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
        'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
    ]
    
    genotype_file = '/scratch1/chris/introner_vis/genotype_matrixes/genotype_matrix_updated_12112025.tsv'
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

    # Count presence (1) for each ortholog in complete data only - this is the DERIVED count
    # since we assume absent (0) is ancestral
    introner_derived_counts = df_complete[df_complete['presence'] == 1].groupby('ortholog_id').size()
    
    # Create unfolded frequency spectrum (1 to n-1, excluding 0 and n)
    max_samples = len(group1_samples)
    unfolded_introner_counts = []
    
    for count in introner_derived_counts.values:
        # Exclude monomorphic sites (all absent = 0, all present = max_samples)
        if count > 0 and count < max_samples:
            unfolded_introner_counts.append(count)
    
    # Create unfolded frequency spectrum
    unfolded_freq_counts = Counter(unfolded_introner_counts)
    
    # Fill in missing frequencies with 0 (frequencies from 1 to max_samples-1)
    unfolded_sfs = [unfolded_freq_counts.get(i, 0) for i in range(1, max_samples)]
    
    print(f"Used {len(introner_derived_counts)} complete ortholog groups")
    print(f"Polymorphic introner groups (1 to {max_samples-1}): {len(unfolded_introner_counts)}")
    print(f"Unfolded introner SFS: {unfolded_sfs}")
    
    return unfolded_sfs, list(range(1, max_samples))

def calculate_folded_introner_frequencies():
    """
    Calculate folded introner frequencies using the comparative SFS approach
    """
    print("\n=== Calculating Folded Introner Frequencies ===")
    
    # Import required modules for introner analysis
    import pandas as pd
    from collections import Counter
    
    # Define Group1 samples (from the comparative analysis)
    group1_samples = [
        'CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
        'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833'
    ]
    
    genotype_file = '/scratch1/chris/introner_vis/genotype_matrixes/genotype_matrix_updated_12112025.tsv'
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
    
    # Apply folding to introner counts
    max_samples = len(group1_samples)
    folded_introner_counts = []
    
    for count in introner_counts.values:
        folded_count = fold_frequency(count, max_samples)
        folded_introner_counts.append(folded_count)
    
    # Create folded frequency spectrum
    folded_freq_counts = Counter(folded_introner_counts)
    
    # Fill in missing frequencies with 0 (frequencies from 1 to max_samples//2)
    max_folded_freq = max_samples // 2
    folded_sfs = [folded_freq_counts.get(i, 0) for i in range(1, max_folded_freq + 1)]
    
    print(f"Used {len(introner_counts)} complete ortholog groups for folded introner SFS")
    print(f"Folded introner SFS: {folded_sfs}")
    
    return folded_sfs, list(range(1, max_folded_freq + 1))

def perform_unfolded_sfs_analysis(allele_freq_df):
    """
    Perform unfolded Site Frequency Spectrum analysis
    """
    print("\n=== Unfolded Site Frequency Spectrum Analysis ===")
    
    # Filter for coding variants (synonymous, missense, nonsense)
    coding_variants = allele_freq_df[allele_freq_df['mutation_type'].isin(['synonymous', 'missense', 'nonsense'])].copy()
    print(f"Coding variants for SFS: {len(coding_variants)}")
    
    # Calculate unfolded introner frequencies
    unfolded_introner_sfs, introner_freq_bins = calculate_unfolded_introner_frequencies()
    
    # Calculate unfolded synonymous and nonsynonymous site frequencies
    syn_variants = coding_variants[coding_variants['mutation_type'] == 'synonymous']
    nonsyn_variants = coding_variants[coding_variants['mutation_type'].isin(['missense', 'nonsense'])]
    
    # Create unfolded frequency counts for synonymous variants
    syn_unfolded_counts = Counter()
    for _, row in syn_variants.iterrows():
        derived_count = int(row['derived_count'])
        if derived_count > 0:  # Exclude monomorphic sites (already excluded in calculation)
            syn_unfolded_counts[derived_count] += 1
    
    # Create unfolded frequency counts for nonsynonymous variants
    nonsyn_unfolded_counts = Counter()
    for _, row in nonsyn_variants.iterrows():
        derived_count = int(row['derived_count'])
        if derived_count > 0:  # Exclude monomorphic sites
            nonsyn_unfolded_counts[derived_count] += 1
    
    # Create unfolded SFS arrays for synonymous and nonsynonymous  
    # Note: introner_freq_bins goes from 1 to n, but we use 1 to n-1 for the actual SFS
    max_unfolded_freq = len(unfolded_introner_sfs) + 1  # This should match the introner SFS length
    unfolded_syn_sfs = [syn_unfolded_counts.get(i, 0) for i in range(1, max_unfolded_freq)]
    unfolded_nonsyn_sfs = [nonsyn_unfolded_counts.get(i, 0) for i in range(1, max_unfolded_freq)]
    
    print(f"Unfolded synonymous SFS: {unfolded_syn_sfs}")
    print(f"Unfolded nonsynonymous SFS: {unfolded_nonsyn_sfs}")
    
    # Analyze unfolded frequency distributions
    print("\n--- Unfolded Frequency Distribution Summary ---")
    print(f"Total polymorphic introners: {sum(unfolded_introner_sfs)}")
    print(f"Total synonymous mutations: {sum(unfolded_syn_sfs)}")
    print(f"Total nonsynonymous mutations: {sum(unfolded_nonsyn_sfs)}")
    
    # Create comprehensive unfolded SFS visualization
    # Make sure all arrays have the same length
    freq_bins_plot = list(range(1, len(unfolded_introner_sfs) + 1))
    create_unfolded_sfs_plots(unfolded_introner_sfs, unfolded_syn_sfs, unfolded_nonsyn_sfs, freq_bins_plot)
    
    # Save unfolded results
    unfolded_results_df = pd.DataFrame({
        'derived_frequency_count': freq_bins_plot,
        'introners': unfolded_introner_sfs,
        'synonymous_mutations': unfolded_syn_sfs,
        'nonsynonymous_mutations': unfolded_nonsyn_sfs
    })
    
    output_file = '/scratch1/chris/introner_vis/polymorphism_analysis/unfolded_sfs_analysis_results.tsv'
    unfolded_results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nUnfolded SFS results saved to: {output_file}")
    
    return coding_variants, unfolded_results_df

def perform_sfs_analysis(allele_freq_df):
    """
    Perform Site Frequency Spectrum analysis using folded frequencies
    """
    print("\n=== Folded Site Frequency Spectrum Analysis ===")
    
    # Filter for coding variants (synonymous, missense, nonsense)
    coding_variants = allele_freq_df[allele_freq_df['mutation_type'].isin(['synonymous', 'missense', 'nonsense'])].copy()
    print(f"Coding variants for SFS: {len(coding_variants)}")
    
    # Calculate folded introner frequencies
    folded_introner_sfs, introner_freq_bins = calculate_folded_introner_frequencies()
    
    # Calculate folded synonymous and nonsynonymous site frequencies
    syn_variants = coding_variants[coding_variants['mutation_type'] == 'synonymous']
    nonsyn_variants = coding_variants[coding_variants['mutation_type'].isin(['missense', 'nonsense'])]
    
    # Create folded frequency counts for synonymous variants
    syn_folded_counts = Counter()
    for _, row in syn_variants.iterrows():
        folded_count = int(row['folded_count'])
        if folded_count > 0:  # Exclude monomorphic sites
            syn_folded_counts[folded_count] += 1
    
    # Create folded frequency counts for nonsynonymous variants
    nonsyn_folded_counts = Counter()
    for _, row in nonsyn_variants.iterrows():
        folded_count = int(row['folded_count'])
        if folded_count > 0:  # Exclude monomorphic sites
            nonsyn_folded_counts[folded_count] += 1
    
    # Create folded SFS arrays for synonymous and nonsynonymous
    max_folded_freq = max(introner_freq_bins)
    folded_syn_sfs = [syn_folded_counts.get(i, 0) for i in range(1, max_folded_freq + 1)]
    folded_nonsyn_sfs = [nonsyn_folded_counts.get(i, 0) for i in range(1, max_folded_freq + 1)]
    
    print(f"Folded synonymous SFS: {folded_syn_sfs}")
    print(f"Folded nonsynonymous SFS: {folded_nonsyn_sfs}")
    
    # Analyze folded frequency distributions
    print("\n--- Folded Frequency Distribution Summary ---")
    print(f"Total introners: {sum(folded_introner_sfs)}")
    print(f"Total synonymous mutations: {sum(folded_syn_sfs)}")
    print(f"Total nonsynonymous mutations: {sum(folded_nonsyn_sfs)}")
    
    # Create comprehensive folded SFS visualization
    create_folded_sfs_plots(folded_introner_sfs, folded_syn_sfs, folded_nonsyn_sfs, introner_freq_bins)
    
    # Save folded results
    folded_results_df = pd.DataFrame({
        'folded_frequency': introner_freq_bins,
        'introners': folded_introner_sfs,
        'synonymous_mutations': folded_syn_sfs,
        'nonsynonymous_mutations': folded_nonsyn_sfs
    })
    
    output_file = '/scratch1/chris/introner_vis/polymorphism_analysis/folded_sfs_analysis_results.tsv'
    folded_results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nFolded SFS results saved to: {output_file}")
    
    return coding_variants, folded_results_df

def create_unfolded_sfs_plots(introner_sfs, syn_sfs, nonsyn_sfs, freq_bins):
    """
    Create unfolded Site Frequency Spectrum visualizations with grouped bars
    """
    print("\nCreating unfolded SFS visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, orange, green
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Unfolded Site Frequency Spectrum Analysis', fontsize=16)
    
    # 1. Grouped bar plot (main plot requested by user)
    x = np.array(freq_bins)
    width = 0.25
    
    bars1 = ax1.bar(x - width, introner_sfs, width, label='Introners (derived=present)', 
                   color=colors[0], alpha=0.8)
    bars2 = ax1.bar(x, syn_sfs, width, label='Synonymous mutations', 
                   color=colors[1], alpha=0.8)
    bars3 = ax1.bar(x + width, nonsyn_sfs, width, label='Nonsynonymous mutations', 
                   color=colors[2], alpha=0.8)
    
    ax1.set_xlabel('Derived Allele Frequency (number of samples)')
    ax1.set_ylabel('Number of variants/introners')
    ax1.set_title('Unfolded SFS: Grouped Bar Plot')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(x)
    
    # Add value labels on bars
    def add_value_labels(ax, bars):
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{int(height)}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),  # 3 points vertical offset
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)
    
    add_value_labels(ax1, bars1)
    add_value_labels(ax1, bars2)
    add_value_labels(ax1, bars3)
    
    # 2. Log-scale comparison plot
    # Convert zeros to small values for log scale
    introner_sfs_log = [max(x, 0.1) for x in introner_sfs]
    syn_sfs_log = [max(x, 0.1) for x in syn_sfs]
    nonsyn_sfs_log = [max(x, 0.1) for x in nonsyn_sfs]
    
    ax2.bar(x - width, introner_sfs_log, width, label='Introners (derived=present)', 
           color=colors[0], alpha=0.8)
    ax2.bar(x, syn_sfs_log, width, label='Synonymous mutations', 
           color=colors[1], alpha=0.8)
    ax2.bar(x + width, nonsyn_sfs_log, width, label='Nonsynonymous mutations', 
           color=colors[2], alpha=0.8)
    
    ax2.set_xlabel('Derived Allele Frequency (number of samples)')
    ax2.set_ylabel('Number of variants/introners (log scale)')
    ax2.set_title('Unfolded SFS: Log Scale Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(x)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = '/scratch1/chris/introner_vis/polymorphism_analysis/unfolded_sfs_analysis_plots.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Unfolded SFS plots saved to: {plot_file}")
    
    plt.close()

    # Create density plot for comparison
    # Half-page width for publication (single column width ~3.5-4 inches)
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Calculate densities (normalize by total counts)
    total_introners = sum(introner_sfs)
    total_syn = sum(syn_sfs)
    total_nonsyn = sum(nonsyn_sfs)

    if total_introners > 0:
        introner_density = [x / total_introners for x in introner_sfs]
    else:
        introner_density = introner_sfs

    if total_syn > 0:
        syn_density = [x / total_syn for x in syn_sfs]
    else:
        syn_density = syn_sfs

    if total_nonsyn > 0:
        nonsyn_density = [x / total_nonsyn for x in nonsyn_sfs]
    else:
        nonsyn_density = nonsyn_sfs

    # Use seaborn color palette for cleaner colors
    clean_colors = sns.color_palette("Set2", 3)

    ax.bar(x - width, introner_density, width, label=f'Introners (n={total_introners})',
           color=clean_colors[0], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.bar(x, syn_density, width, label=f'Synonymous (n={total_syn})',
           color=clean_colors[1], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.bar(x + width, nonsyn_density, width, label=f'Nonsynonymous (n={total_nonsyn})',
           color=clean_colors[2], alpha=0.8, edgecolor='black', linewidth=1.2)

    ax.set_xlabel('Derived Allele Frequency (number of samples)', fontsize=10)
    ax.set_ylabel('Density (proportion of total)', fontsize=10)
    ax.set_title('Unfolded SFS: Density Comparison', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, frameon=True, edgecolor='black')
    ax.set_xticks(x)
    ax.tick_params(labelsize=9)

    # Remove background grid
    ax.grid(False)

    # Clean up spines for publication quality
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    
    # Save density plot
    density_plot_file = '/scratch1/chris/introner_vis/polymorphism_analysis/unfolded_sfs_density_plots.png'
    plt.savefig(density_plot_file, dpi=300, bbox_inches='tight')
    print(f"Unfolded SFS density plots saved to: {density_plot_file}")
    
    plt.close()

def create_folded_sfs_plots(introner_sfs, syn_sfs, nonsyn_sfs, freq_bins):
    """
    Create folded Site Frequency Spectrum visualizations with grouped bars
    """
    print("\nCreating folded SFS visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, orange, green
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Folded Site Frequency Spectrum Analysis', fontsize=16)
    
    # 1. Grouped bar plot (main plot requested by user)
    x = np.array(freq_bins)
    width = 0.25
    
    bars1 = ax1.bar(x - width, introner_sfs, width, label='Introners', 
                   color=colors[0], alpha=0.8)
    bars2 = ax1.bar(x, syn_sfs, width, label='Synonymous mutations', 
                   color=colors[1], alpha=0.8)
    bars3 = ax1.bar(x + width, nonsyn_sfs, width, label='Nonsynonymous mutations', 
                   color=colors[2], alpha=0.8)
    
    ax1.set_xlabel('Folded Allele Frequency (number of samples)')
    ax1.set_ylabel('Number of variants/introners')
    ax1.set_title('Folded SFS: Grouped Bar Plot')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(x)
    
    # Add value labels on bars
    def add_value_labels(ax, bars):
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{int(height)}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),  # 3 points vertical offset
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)
    
    add_value_labels(ax1, bars1)
    add_value_labels(ax1, bars2)
    add_value_labels(ax1, bars3)
    
    # 2. Log-scale comparison plot
    # Convert zeros to small values for log scale
    introner_sfs_log = [max(x, 0.1) for x in introner_sfs]
    syn_sfs_log = [max(x, 0.1) for x in syn_sfs]
    nonsyn_sfs_log = [max(x, 0.1) for x in nonsyn_sfs]
    
    ax2.bar(x - width, introner_sfs_log, width, label='Introners', 
           color=colors[0], alpha=0.8)
    ax2.bar(x, syn_sfs_log, width, label='Synonymous mutations', 
           color=colors[1], alpha=0.8)
    ax2.bar(x + width, nonsyn_sfs_log, width, label='Nonsynonymous mutations', 
           color=colors[2], alpha=0.8)
    
    ax2.set_xlabel('Folded Allele Frequency (number of samples)')
    ax2.set_ylabel('Number of variants/introners (log scale)')
    ax2.set_title('Folded SFS: Log Scale Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(x)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = '/scratch1/chris/introner_vis/polymorphism_analysis/folded_sfs_analysis_plots.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Folded SFS plots saved to: {plot_file}")
    
    plt.close()

    # Create density plot for comparison
    # Half-page width for publication (single column width ~3.5-4 inches)
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Calculate densities (normalize by total counts)
    total_introners = sum(introner_sfs)
    total_syn = sum(syn_sfs)
    total_nonsyn = sum(nonsyn_sfs)

    if total_introners > 0:
        introner_density = [x / total_introners for x in introner_sfs]
    else:
        introner_density = introner_sfs

    if total_syn > 0:
        syn_density = [x / total_syn for x in syn_sfs]
    else:
        syn_density = syn_sfs

    if total_nonsyn > 0:
        nonsyn_density = [x / total_nonsyn for x in nonsyn_sfs]
    else:
        nonsyn_density = nonsyn_sfs

    # Use seaborn color palette for cleaner colors
    clean_colors = sns.color_palette("Set2", 3)

    ax.bar(x - width, introner_density, width, label=f'Introners (n={total_introners})',
           color=clean_colors[0], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.bar(x, syn_density, width, label=f'Synonymous (n={total_syn})',
           color=clean_colors[1], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.bar(x + width, nonsyn_density, width, label=f'Nonsynonymous (n={total_nonsyn})',
           color=clean_colors[2], alpha=0.8, edgecolor='black', linewidth=1.2)

    ax.set_xlabel('Folded Allele Frequency (number of samples)', fontsize=10)
    ax.set_ylabel('Density (proportion of total)', fontsize=10)
    ax.set_title('Folded SFS: Density Comparison', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, frameon=True, edgecolor='black')
    ax.set_xticks(x)
    ax.tick_params(labelsize=9)

    # Remove background grid
    ax.grid(False)

    # Clean up spines for publication quality
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save density plot
    density_plot_file = '/scratch1/chris/introner_vis/polymorphism_analysis/folded_sfs_density_plots.png'
    plt.savefig(density_plot_file, dpi=300, bbox_inches='tight')
    print(f"Folded SFS density plots saved to: {density_plot_file}")
    
    plt.close()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Site Frequency Spectrum Analysis with unfolded (default) or folded options")
    parser.add_argument('--folded', action='store_true', help='Use folded SFS analysis (default is unfolded)')
    parser.add_argument('--vcf', type=str, default='/scratch1/chris/introner_vis/data_prep/snpEff_run2/mpusilla.snps.snpEff.vcf',
                       help='Path to SnpEff annotated VCF file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    if args.folded:
        print("=== Running FOLDED Site Frequency Spectrum Analysis ===")
        # Calculate allele frequencies (folded)
        allele_freq_df = calculate_allele_frequencies()
        
        # Perform SFS analysis
        coding_variants, results_df = perform_sfs_analysis(allele_freq_df)
        
        if coding_variants is not None and results_df is not None:
            print(f"\n✓ Folded Site Frequency Spectrum analysis completed successfully")
            print(f"✓ Results saved to folded_sfs_analysis_results.tsv")
            print(f"✓ Plots saved to folded_sfs_analysis_plots.png and folded_sfs_density_plots.png")
        else:
            print("\n✗ SFS analysis failed")
            sys.exit(1)
    else:
        print("=== Running UNFOLDED Site Frequency Spectrum Analysis ===")
        # Calculate unfolded allele frequencies (default)
        allele_freq_df = calculate_unfolded_allele_frequencies(args.vcf)
        
        # Perform unfolded SFS analysis
        coding_variants, results_df = perform_unfolded_sfs_analysis(allele_freq_df)
        
        if coding_variants is not None and results_df is not None:
            print(f"\n✓ Unfolded Site Frequency Spectrum analysis completed successfully")
            print(f"✓ Results saved to unfolded_sfs_analysis_results.tsv")
            print(f"✓ Plots saved to unfolded_sfs_analysis_plots.png and unfolded_sfs_density_plots.png")
        else:
            print("\n✗ SFS analysis failed")
            sys.exit(1)
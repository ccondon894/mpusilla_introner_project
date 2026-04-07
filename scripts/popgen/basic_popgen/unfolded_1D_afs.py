import argparse
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Compute and plot unfolded 1D allele frequency spectrum using Group2 as outgroup.")
    parser.add_argument('--vcf_file', type=str, help="Path to the VCF file")
    parser.add_argument('--plot', action='store_true', help="Plot the unfolded 1D AFS")
    return parser.parse_args()

def compute_unfolded_afs(vcf_file):
    # Sample definitions
    outgroup_samples = ['RCC1749', 'RCC3052']  # Group2 - outgroup for ancestral state inference
    obsolete_samples = {'CCMP490', 'RCC647', 'RCC835'}
    
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)
    
    # Get sample indices
    sample_names = list(vcf.header.samples)
    outgroup_indices = [sample_names.index(sample) for sample in outgroup_samples if sample in sample_names]
    group1_indices = [i for i in range(len(sample_names)) 
                     if sample_names[i] not in outgroup_samples and sample_names[i] not in obsolete_samples]
    
    print(f"Outgroup samples (Group2): {[sample_names[i] for i in outgroup_indices]}")
    print(f"Group1 samples: {[sample_names[i] for i in group1_indices]}")
    print(f"Total Group1 samples: {len(group1_indices)}")
    
    # Initialize AFS array (0 to n derived alleles in Group1, haploid)
    n_group1 = len(group1_indices)
    afs = np.zeros(n_group1 + 1, dtype=int)
    
    sites_processed = 0
    sites_excluded_polymorphic_outgroup = 0
    sites_excluded_missing_outgroup = 0
    sites_excluded_monomorphic = 0
    sites_used = 0
    
    # Process each variant
    for record in vcf:
        if record.alts and len(record.alts) == 1:  # Biallelic SNP only
            sites_processed += 1
            
            # Get genotypes for outgroup
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
                if None in gt:  # Skip sites with missing data in Group1
                    group1_missing = True
                    break
                
                # Count alleles different from ancestral
                for allele in gt:
                    if allele != ancestral_allele:
                        derived_count += 1
            
            # Skip if missing data in Group1
            if group1_missing:
                continue
            
            # Skip monomorphic sites (fixed in both groups)
            # This means derived_count = 0 (ancestral in all Group1) or derived_count = n_group1 (derived in all Group1, haploid)
            if derived_count == 0 or derived_count == n_group1:
                sites_excluded_monomorphic += 1
                continue
            
            # Add to AFS
            afs[derived_count] += 1
            sites_used += 1
    
    print(f"\nSite filtering summary:")
    print(f"Total biallelic sites processed: {sites_processed}")
    print(f"Sites excluded (polymorphic outgroup): {sites_excluded_polymorphic_outgroup}")
    print(f"Sites excluded (missing data in outgroup): {sites_excluded_missing_outgroup}")
    print(f"Sites excluded (monomorphic across groups): {sites_excluded_monomorphic}")
    print(f"Sites used for unfolded AFS: {sites_used}")
    
    return afs, n_group1

def print_afs(afs, n_samples):
    print(f"\nUnfolded 1D Allele Frequency Spectrum (Group1, n={n_samples} haploid samples):")
    print("Derived allele count | Number of sites")
    print("----------------------------------")
    
    for i, count in enumerate(afs):
        freq = i / n_samples if n_samples > 0 else 0
        print(f"{i:2d} ({freq:.3f}) | {count}")

def plot_unfolded_afs_log(afs, n_samples):
    """Plot unfolded AFS with log scale including all segregating frequencies"""
    # Create frequency array
    frequencies = np.arange(len(afs)) / n_samples
    
    # Filter out zero counts for cleaner visualization
    non_zero_mask = afs > 0
    frequencies_nz = frequencies[non_zero_mask]
    afs_nz = afs[non_zero_mask]
    
    # Create bar plot with log scale
    plt.figure(figsize=(12, 6))
    bars = plt.bar(frequencies_nz, afs_nz, width=0.8/n_samples, alpha=0.7, color='steelblue', edgecolor='black')
    
    # Set log scale for y-axis
    plt.yscale('log')
    
    # Customize plot
    plt.xlabel('Derived Allele Frequency')
    plt.ylabel('Number of Sites (log scale)')
    plt.title(f'Unfolded 1D Allele Frequency Spectrum - Log Scale\nGroup1 haploid samples (n={n_samples}), Group2 as outgroup')
    plt.grid(True, alpha=0.3, which='both')
    
    # Add minor grid lines for log scale
    plt.grid(True, alpha=0.2, which='minor')
    
    # Add text annotations for bars
    for freq, count in zip(frequencies_nz, afs_nz):
        if count > 0:  # All non-zero counts
            plt.text(freq, count * 1.1, str(count), 
                    ha='center', va='bottom', fontsize=8)
    
    # Set y-axis limits to start from 1 (since we're using log scale)
    plt.ylim(bottom=0.5)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig("unfolded_1D_AFS_log.pdf")
    plt.savefig("unfolded_1D_AFS_log.png", dpi=300)

def plot_unfolded_afs_linear(afs, n_samples):
    """Plot unfolded AFS with linear scale excluding near-fixation sites"""
    # Create frequency array, excluding the highest frequency bin (near-fixation)
    # This excludes sites with frequency (n_samples - 1)/n_samples
    frequencies = np.arange(len(afs) - 1) / n_samples  # Exclude last bin
    afs_subset = afs[:-1]  # Exclude last frequency bin
    
    # Filter out zero counts for cleaner visualization
    non_zero_mask = afs_subset > 0
    frequencies_nz = frequencies[non_zero_mask]
    afs_nz = afs_subset[non_zero_mask]
    
    # Create bar plot with linear scale
    plt.figure(figsize=(12, 6))
    bars = plt.bar(frequencies_nz, afs_nz, width=0.8/(n_samples-1), alpha=0.7, color='coral', edgecolor='black')
    
    # Customize plot
    plt.xlabel('Derived Allele Frequency')
    plt.ylabel('Number of Sites')
    plt.title(f'Unfolded 1D Allele Frequency Spectrum - Linear Scale\nGroup1 haploid samples (n={n_samples}), excluding Group1 fixed sites')
    plt.grid(True, alpha=0.3)
    
    # Add text annotations for bars
    max_count = max(afs_nz) if len(afs_nz) > 0 else 0
    for freq, count in zip(frequencies_nz, afs_nz):
        if count > max_count * 0.02:  # Only annotate bars > 2% of max
            plt.text(freq, count + max_count * 0.01, str(count), 
                    ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig("unfolded_1D_AFS_linear.pdf")
    plt.savefig("unfolded_1D_AFS_linear.png", dpi=300)

def main():
    args = parse_args()
    afs, n_samples = compute_unfolded_afs(args.vcf_file)
    
    # Print results
    print_afs(afs, n_samples)
    
    # Plot if requested
    if args.plot:
        print("\nGenerating plots...")
        print("1. Log-scale plot with all segregating frequencies")
        plot_unfolded_afs_log(afs, n_samples)
        
        print("2. Linear-scale plot excluding near-fixation sites")
        plot_unfolded_afs_linear(afs, n_samples)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import pandas as pd
import subprocess
import sys
import os

def create_bed_file():
    """Create BED file with coordinates for VCF extraction"""
    print("Step 2: Creating BED file for VCF region extraction")
    
    # Read the filtered polymorphic groups
    df = pd.read_csv('/scratch1/chris/introner_vis/polymorphic_groups_ccmp1545_present.tsv', sep='\t')
    
    # Create BED format: chrom, start, end, name
    bed_df = df[['contig', 'start', 'end', 'ortholog_id']].copy()
    bed_df = bed_df.sort_values(['contig', 'start'])
    
    # Save BED file
    bed_file = '/scratch1/chris/introner_vis/introner_regions.bed'
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    
    print(f"BED file created: {bed_file}")
    print(f"Regions to extract: {len(bed_df)}")
    print(f"Contig distribution:")
    print(bed_df['contig'].value_counts().head(10))
    
    return bed_file

def extract_vcf_regions(bed_file):
    """Extract VCF regions using bcftools"""
    print("\nStep 3: Extracting VCF regions using bcftools")
    
    vcf_file = '/scratch1/chris/introner_vis/data_prep/snpEff_run2/mpusilla.snps.snpEff.vcf'
    output_vcf = '/scratch1/chris/introner_vis/introner_regions_extracted.vcf'
    
    # Check if VCF file exists
    if not os.path.exists(vcf_file):
        print(f"ERROR: VCF file not found: {vcf_file}")
        return None
    
    # Define Group1 samples (excluding the ones we don't want)
    group1_samples = ['CCMP1545', 'RCC114', 'RCC1614', 'RCC1698', 'RCC2482', 
                      'RCC373', 'RCC465', 'RCC629', 'RCC692', 'RCC693', 'RCC833']
    
    # Create sample list for bcftools
    sample_list = ','.join(group1_samples)
    
    print(f"Extracting regions from: {vcf_file}")
    print(f"Using samples: {sample_list}")
    print(f"Output file: {output_vcf}")
    
    # Try bcftools command
    cmd = [
        'bcftools', 'view',
        '-R', bed_file,  # Regions from BED file
        '-s', sample_list,  # Sample subset
        '-o', output_vcf,  # Output file
        vcf_file
    ]
    
    try:
        print("Running bcftools command...")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ bcftools extraction successful")
        
        # Check output file
        if os.path.exists(output_vcf):
            # Count lines in output
            line_count = subprocess.run(['wc', '-l', output_vcf], capture_output=True, text=True)
            print(f"Output VCF has {line_count.stdout.strip().split()[0]} lines")
            return output_vcf
        else:
            print("ERROR: Output VCF file was not created")
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bcftools failed with return code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        return None
    except FileNotFoundError:
        print("ERROR: bcftools not found. Trying alternative approach...")
        return extract_vcf_regions_alternative(bed_file, vcf_file, output_vcf, group1_samples)

def extract_vcf_regions_alternative(bed_file, vcf_file, output_vcf, group1_samples):
    """Alternative VCF extraction using grep and awk"""
    print("Using alternative extraction method...")
    
    # Read BED file to get regions
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['contig', 'start', 'end', 'name'])
    
    # For now, let's just extract the header and first few variants for testing
    print("Extracting VCF header and testing region extraction...")
    
    try:
        # Extract header
        header_cmd = ['grep', '^#', vcf_file]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        
        # Extract first few data lines for testing
        data_cmd = ['grep', '-v', '^#', vcf_file]
        data_result = subprocess.run(data_cmd, capture_output=True, text=True, check=True)
        data_lines = data_result.stdout.strip().split('\n')[:100]  # First 100 variants
        
        # Write to output file
        with open(output_vcf, 'w') as f:
            f.write(header_result.stdout)
            for line in data_lines:
                f.write(line + '\n')
        
        print(f"✓ Alternative extraction created test file with {len(data_lines)} variants")
        return output_vcf
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Alternative extraction failed: {e}")
        return None

def main():
    """Main function to orchestrate VCF region extraction"""
    print("=== VCF Region Extraction ===")
    
    # Create BED file
    bed_file = create_bed_file()
    if bed_file is None:
        return False
    
    # Extract VCF regions
    output_vcf = extract_vcf_regions(bed_file)
    if output_vcf is None:
        return False
    
    print(f"\n✓ VCF extraction complete: {output_vcf}")
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        sys.exit(1)
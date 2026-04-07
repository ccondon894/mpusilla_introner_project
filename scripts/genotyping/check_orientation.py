#!/usr/bin/env python3

import pandas as pd
import subprocess
import tempfile
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import multiprocessing
from functools import partial
import argparse

def extract_sequence_from_fasta(fasta_file, contig, start, end):
    """Extract sequence from a FASTA file given coordinates"""
    try:
        # Use samtools faidx to extract the sequence
        region = f"{contig}:{start}-{end}"
        result = subprocess.run(['samtools', 'faidx', fasta_file, region], 
                              capture_output=True, text=True, check=True)
        
        # Parse the output to get just the sequence (skip header line)
        lines = result.stdout.strip().split('\n')
        sequence = ''.join(lines[1:])  # Skip the FASTA header
        return sequence
    except subprocess.CalledProcessError as e:
        print(f"Error extracting sequence from {fasta_file} at {region}: {e}")
        return None

def run_mafft_with_direction(sequences, ortholog_id):
    """Run MAFFT with --adjustdirection and return which sequences were reversed"""
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_input:
        temp_input_path = temp_input.name
        for sample, seq in sequences.items():
            temp_input.write(f">{sample}\n{seq}\n")
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_output:
        temp_output_path = temp_output.name
    
    try:
        # Run MAFFT with --adjustdirection
        with open(temp_output_path, 'w') as outfile:
            result = subprocess.run([
                'mafft', '--adjustdirection', '--quiet', '--maxiterate', '1000', 
                '--globalpair', temp_input_path
            ], stdout=outfile, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            print(f"MAFFT failed for {ortholog_id}: {result.stderr}")
            return {}
        
        # Parse the output to check for _R_ sequences (reversed by MAFFT)
        reversed_samples = {}
        with open(temp_output_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sample_name = record.id
                # MAFFT adds _R_ to sequence names that were reverse complemented
                if '_R_' in sample_name:
                    # Extract original sample name (remove _R_ suffix)
                    original_name = sample_name.replace('_R_', '')
                    reversed_samples[original_name] = True
                else:
                    reversed_samples[sample_name] = False
        
        return reversed_samples
        
    finally:
        # Clean up temporary files
        if os.path.exists(temp_input_path):
            os.unlink(temp_input_path)
        if os.path.exists(temp_output_path):
            os.unlink(temp_output_path)

def process_ortholog_group(ortholog_data, fasta_dir, group1_samples, group2_samples):
    """Process a single ortholog group to determine orientations"""
    
    ortholog_id, ortholog_df = ortholog_data
    
    # Filter out samples with missing data (presence == 3)
    valid_samples_df = ortholog_df[ortholog_df['presence'] != 3]
    
    if len(valid_samples_df) < 2:
        print(f"Skipping {ortholog_id}: fewer than 2 valid samples")
        return []
    
    # Choose reference sample: prefer CCMP1545, then other Group1, then Group2
    reference_sample = None
    if 'CCMP1545' in valid_samples_df['sample'].values:
        reference_sample = 'CCMP1545'
    else:
        # Look for any Group1 sample
        group1_present = [s for s in group1_samples if s in valid_samples_df['sample'].values]
        if group1_present:
            reference_sample = group1_present[0]  # Pick first available Group1 sample
        else:
            # Fall back to Group2 sample
            group2_present = [s for s in group2_samples if s in valid_samples_df['sample'].values]
            if group2_present:
                reference_sample = group2_present[0]
            else:
                print(f"Skipping {ortholog_id}: no valid reference sample found")
                return []
    
    print(f"Processing {ortholog_id} with reference sample {reference_sample}...")
    
    # Extract sequences for all valid samples (including absent ones)
    sequences = {}
    skip_group = False
    
    for _, row in valid_samples_df.iterrows():
        sample = row['sample']
        contig = row['contig']
        start = int(row['start'])
        end = int(row['end'])
        presence = row['presence']
        
        # Build path to FASTA file
        fasta_file = os.path.join(fasta_dir, f"{sample}.vg_paths.fa")
        
        if not os.path.exists(fasta_file):
            print(f"Warning: FASTA file not found: {fasta_file}")
            skip_group = True
            break
        
        # Extract sequence for both present (1) and absent (2) introners
        # Absent introners will have flanking sequence with gap in middle
        sequence = extract_sequence_from_fasta(fasta_file, contig, start, end)
        if sequence is None:
            print(f"Failed to extract sequence for {sample} in {ortholog_id}")
            skip_group = True
            break
        
        sequences[sample] = sequence
    
    if skip_group or len(sequences) < 2:
        print(f"Skipping {ortholog_id}: insufficient sequences for alignment")
        return []
    
    # Ensure reference sample has a sequence
    if reference_sample not in sequences:
        print(f"Skipping {ortholog_id}: reference sample {reference_sample} has no sequence")
        return []
    
    # Run MAFFT to determine orientations
    orientations = run_mafft_with_direction(sequences, ortholog_id)
    
    if not orientations:
        print(f"Failed to determine orientations for {ortholog_id}")
        return []
    
    # Check reference sample orientation and adjust relative to it
    reference_reversed = orientations.get(reference_sample, False)
    
    # Record results for ALL samples in the ortholog group (including those with missing data)
    results = []
    for _, row in ortholog_df.iterrows():
        sample = row['sample']
        presence = row['presence']
        
        if presence == 3:  # Missing data
            # Cannot determine orientation for missing data
            orientation = None
            left_reverse = None
            right_reverse = None
        else:  # Present (1) or absent (2) introner - both have sequences
            sample_reversed = orientations.get(sample, False)
            # Determine if sample should be reversed relative to reference
            reverse_relative_to_reference = sample_reversed != reference_reversed
            
            # Set orientation columns based on reverse_relative_to_reference
            if reverse_relative_to_reference:
                orientation = "reverse"
                left_reverse = True
                right_reverse = True
            else:
                orientation = "forward"
                left_reverse = False
                right_reverse = False
        
        # Copy all original columns and add new orientation columns
        result = row.to_dict()
        result['orientation'] = orientation
        result['left_reverse'] = left_reverse
        result['right_reverse'] = right_reverse
        
        results.append(result)
    
    return results

def determine_orientations(genotype_matrix_file, fasta_dir, output_matrix_file, n_threads=12):
    """Main function to determine orientations relative to a reference sample using multiprocessing"""
    
    # Read genotype matrix
    df = pd.read_csv(genotype_matrix_file, sep='\t')
    
    # Define sample groups
    group1_samples = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482", 
                      "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
    group2_samples = ["RCC1749", "RCC3052"]
    
    # Group data by ortholog_id
    ortholog_groups = [(ortholog_id, ortholog_df) for ortholog_id, ortholog_df in df.groupby('ortholog_id')]
    
    print(f"Processing {len(ortholog_groups)} ortholog groups using {n_threads} threads...")
    
    # Create partial function with fixed arguments
    process_func = partial(process_ortholog_group, 
                          fasta_dir=fasta_dir, 
                          group1_samples=group1_samples, 
                          group2_samples=group2_samples)
    
    # Process ortholog groups in parallel
    with multiprocessing.Pool(processes=n_threads) as pool:
        results_nested = pool.map(process_func, ortholog_groups)
    
    # Flatten results
    results = []
    for result_list in results_nested:
        results.extend(result_list)
    
    if not results:
        print("No results generated!")
        return
    
    # Create results DataFrame with all original columns plus new orientation columns
    results_df = pd.DataFrame(results)
    
    # Ensure proper column order (original columns first, then orientation columns)
    original_columns = df.columns.tolist()
    new_columns = ['orientation', 'left_reverse', 'right_reverse']
    
    # Reorder columns to match original + new columns
    column_order = original_columns + new_columns
    results_df = results_df[column_order]
    
    # Fix data types for start, end, and family columns to avoid .0 suffixes
    numeric_columns = ['start', 'end', 'family']
    for col in numeric_columns:
        if col in results_df.columns:
            # Convert to string, handling NaN values and removing .0 suffixes
            results_df[col] = results_df[col].apply(
                lambda x: str(int(x)) if pd.notna(x) and x != '' else ''
            )
    
    # Write results to file
    results_df.to_csv(output_matrix_file, sep='\t', index=False)
    
    print(f"\nOrientation results written to {output_matrix_file}")
    print(f"Processed {len(results_df['ortholog_id'].unique())} ortholog groups")
    print(f"Total entries: {len(results_df)}")
    
    # Summary statistics
    valid_orientations = results_df[results_df['orientation'].notna()]
    if len(valid_orientations) > 0:
        orientation_counts = valid_orientations['orientation'].value_counts()
        print(f"\nOrientation summary:")
        for orientation, count in orientation_counts.items():
            print(f"  {orientation}: {count} entries ({count/len(valid_orientations)*100:.1f}%)")
    
    # Reference sample usage (based on first non-null entry per ortholog group)
    sample_ref_map = {}
    for ortholog_id, group_df in results_df.groupby('ortholog_id'):
        # Find a sample with valid orientation data to get reference info
        valid_samples = group_df[group_df['orientation'].notna()]
        if not valid_samples.empty:
            # All samples in an ortholog group should have the same reference
            # We can infer the reference from the processing logic
            if 'CCMP1545' in group_df['sample'].values:
                sample_ref_map[ortholog_id] = 'CCMP1545'
            else:
                # Find first Group1 sample, then Group2
                group1_present = [s for s in group1_samples if s in group_df['sample'].values]
                if group1_present:
                    sample_ref_map[ortholog_id] = group1_present[0]
                else:
                    group2_present = [s for s in group2_samples if s in group_df['sample'].values]
                    if group2_present:
                        sample_ref_map[ortholog_id] = group2_present[0]
    
    ref_counts = pd.Series(sample_ref_map.values()).value_counts()
    print(f"\nReference sample usage:")
    for ref_sample, count in ref_counts.items():
        print(f"  {ref_sample}: {count} ortholog groups")
    
    # Presence status summary
    presence_counts = results_df['presence'].value_counts().sort_index()
    print(f"\nPresence status summary:")
    for presence, count in presence_counts.items():
        status = {1: "Present", 2: "Absent", 3: "Missing"}[presence]
        print(f"  {status} (presence={presence}): {count} entries")

def main():
    parser = argparse.ArgumentParser(description='Determine sequence orientations and create corrected genotype matrix')
    parser.add_argument('genotype_matrix', help='Path to input genotype matrix TSV file')
    parser.add_argument('fasta_dir', help='Directory containing sample FASTA files')
    parser.add_argument('output_matrix', help='Output corrected genotype matrix TSV file')
    parser.add_argument('--threads', '-t', type=int, default=12, 
                       help='Number of threads to use (default: 12)')
    
    args = parser.parse_args()
    
    # Call function with corrected arguments
    determine_orientations(args.genotype_matrix, args.fasta_dir, 
                          args.output_matrix, args.threads)

if __name__ == "__main__":
    main()
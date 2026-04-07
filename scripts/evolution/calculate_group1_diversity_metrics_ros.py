#!/usr/bin/env python3

import argparse
import os
import glob
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate Group1 diversity metrics using ratio of sums approach for bottleneck model analysis')
    parser.add_argument('--alignment_dir', help='Directory containing alignment files')
    parser.add_argument('--classification', help='Path to Group1 classification JSON file')
    parser.add_argument('--flank_length', help='Flanking sequence length (100 or 200)')
    parser.add_argument('--output', help='Output file for diversity metrics')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    return parser.parse_args()

def calculate_pairwise_counts(seq1, seq2):
    """Calculate pairwise differences and valid sites between two sequences"""
    def is_valid_base(base):
        return base not in ['-', 'N', 'n']
    
    differences = sum(1 for a, b in zip(seq1, seq2) 
                     if a != b and is_valid_base(a) and is_valid_base(b))
    sites = sum(1 for a, b in zip(seq1, seq2) 
                if is_valid_base(a) and is_valid_base(b))
    return differences, sites

def calculate_group_pi_ros(sequences):
    """Calculate pi within a group using ratio of sums approach"""
    if len(sequences) < 2:
        return None  # Cannot calculate pi with fewer than 2 sequences
    
    n = len(sequences)
    total_differences = 0
    total_sites = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            differences, sites = calculate_pairwise_counts(sequences[i], sequences[j])
            total_differences += differences
            total_sites += sites
    
    return total_differences / total_sites if total_sites > 0 else None

def calculate_between_groups_pi_ros(sequences1, sequences2):
    """Calculate pi between two groups using ratio of sums approach (equivalent to dXY)"""
    if not sequences1 or not sequences2:
        return None
    
    total_differences = 0
    total_sites = 0
    
    for seq1 in sequences1:
        for seq2 in sequences2:
            differences, sites = calculate_pairwise_counts(seq1, seq2)
            total_differences += differences
            total_sites += sites
    
    return total_differences / total_sites if total_sites > 0 else None

def process_group1_orthologs(alignment_dir, classification, flank_length, verbose=False):
    """Process Group1 ortholog groups for frequency-based diversity analysis"""
    results = []
    
    for ortholog_id, info in classification.items():
        frequency = info["frequency"]
        present_count = info["present_count"]
        absent_count = info["absent_count"]
        freq_category = info["frequency_category"]
        
        # Initialize side-specific results
        side_results = {}
        
        for side in ["left", "right"]:
            # Get alignment files
            present_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.freq_{frequency}.present.{side}_flank_{flank_length}bp.mafft.fa")
            absent_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.freq_{frequency}.absent.{side}_flank_{flank_length}bp.mafft.fa")
            combined_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.freq_{frequency}.combined.{side}_flank_{flank_length}bp.mafft.fa")
            
            # Check file existence
            present_exists = os.path.exists(present_file)
            absent_exists = os.path.exists(absent_file)
            combined_exists = os.path.exists(combined_file)
            
            if not present_exists and not absent_exists:
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: No alignment files found")
                continue
            
            try:
                # Initialize metrics
                pi_present = None
                pi_absent = None  
                pi_between = None
                
                # Load sequences
                present_seqs = []
                absent_seqs = []
                
                if present_exists:
                    present_seqs = [str(record.seq) for record in SeqIO.parse(present_file, "fasta")]
                
                if absent_exists:
                    absent_seqs = [str(record.seq) for record in SeqIO.parse(absent_file, "fasta")]
                
                # Calculate pi within present group (if >1 sequence)
                if len(present_seqs) > 1:
                    pi_present = calculate_group_pi_ros(present_seqs)
                
                # Calculate pi within absent group (if >1 sequence)
                if len(absent_seqs) > 1:
                    pi_absent = calculate_group_pi_ros(absent_seqs)
                
                # Calculate pi between present and absent groups
                if present_seqs and absent_seqs and combined_exists:
                    # Use combined alignment for accurate between-group calculations
                    combined_records = list(SeqIO.parse(combined_file, "fasta"))
                    
                    # Get sequence IDs to separate groups
                    present_ids = set(record.id for record in SeqIO.parse(present_file, "fasta")) if present_exists else set()
                    
                    # Separate sequences from combined alignment
                    present_combined_seqs = [str(record.seq) for record in combined_records if record.id in present_ids]
                    absent_combined_seqs = [str(record.seq) for record in combined_records if record.id not in present_ids]
                    
                    if present_combined_seqs and absent_combined_seqs:
                        pi_between = calculate_between_groups_pi_ros(present_combined_seqs, absent_combined_seqs)
                
                elif present_seqs and absent_seqs:
                    # Fallback: direct calculation without combined alignment
                    pi_between = calculate_between_groups_pi_ros(present_seqs, absent_seqs)
                
                # Store results for this side
                side_results[side] = {
                    "pi_present": pi_present,
                    "pi_absent": pi_absent,
                    "pi_between": pi_between
                }
                
            except Exception as e:
                if verbose:
                    print(f"Error processing {ortholog_id} {side}: {e}")
                continue
        
        # Calculate average metrics across sides if we have results
        if side_results:
            pi_present_values = [data["pi_present"] for data in side_results.values() if data["pi_present"] is not None]
            pi_absent_values = [data["pi_absent"] for data in side_results.values() if data["pi_absent"] is not None]
            pi_between_values = [data["pi_between"] for data in side_results.values() if data["pi_between"] is not None]
            
            # Determine singleton analysis categories
            is_singleton_present = (present_count == 1)
            is_singleton_absent = (absent_count == 1)
            
            # For singleton analysis, calculate pi of non-singleton state
            pi_non_singleton = None
            singleton_type = None
            
            if is_singleton_present and absent_count > 1:
                pi_non_singleton = np.mean(pi_absent_values) if pi_absent_values else None
                singleton_type = "singleton_present"
            elif is_singleton_absent and present_count > 1:
                pi_non_singleton = np.mean(pi_present_values) if pi_present_values else None
                singleton_type = "singleton_absent"
            
            results.append({
                "ortholog_id": ortholog_id,
                "frequency": frequency,
                "frequency_category": freq_category,
                "present_count": present_count,
                "absent_count": absent_count,
                "pi_present": np.mean(pi_present_values) if pi_present_values else None,
                "pi_absent": np.mean(pi_absent_values) if pi_absent_values else None,
                "pi_between": np.mean(pi_between_values) if pi_between_values else None,
                "pi_present_std": np.std(pi_present_values) if len(pi_present_values) > 1 else None,
                "pi_absent_std": np.std(pi_absent_values) if len(pi_absent_values) > 1 else None,
                "pi_between_std": np.std(pi_between_values) if len(pi_between_values) > 1 else None,
                "sides_processed": ",".join(side_results.keys()),
                "is_singleton": is_singleton_present or is_singleton_absent,
                "singleton_type": singleton_type,
                "pi_non_singleton": pi_non_singleton,
                "pi_present_absent_diff": (np.mean(pi_present_values) - np.mean(pi_absent_values)) if pi_present_values and pi_absent_values else None,
                "flank_length": flank_length
            })
    
    return results

def main():
    args = parse_arguments()
    
    # Load Group1 classification
    print(f"Loading Group1 classification from {args.classification}")
    with open(args.classification, "r") as f:
        classification = json.load(f)
    
    print(f"Processing {len(classification)} classified Group1 ortholog groups using ratio of sums approach")
    
    # Process ortholog groups
    results = process_group1_orthologs(args.alignment_dir, classification, args.flank_length, args.verbose)
    
    print(f"Processed {len(results)} ortholog groups")
    
    # Print summary statistics
    freq_counts = defaultdict(int) 
    for result in results:
        freq_counts[result["frequency"]] += 1
    
    print("\nFrequency distribution of processed orthologs:")
    for freq in sorted(freq_counts.keys()):
        print(f"  {freq}/11 present: {freq_counts[freq]} ortholog groups")
    
    # Convert results to DataFrame and save
    df = pd.DataFrame(results)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    df.to_csv(args.output, sep="\t", index=False)
    
    print(f"Saved diversity metrics (ratio of sums) for {len(results)} ortholog groups to {args.output}")
    print("Done!")

if __name__ == "__main__":
    main()
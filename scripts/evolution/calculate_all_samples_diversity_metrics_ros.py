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
    parser = argparse.ArgumentParser(description='Calculate all-samples diversity metrics using ratio of sums approach for group comparison analysis')
    parser.add_argument('--alignment_dir', help='Directory containing alignment files')
    parser.add_argument('--classification', help='Path to all-samples classification JSON file')
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

def calculate_between_groups_dxy_ros(sequences1, sequences2):
    """Calculate dXY between two groups using ratio of sums approach"""
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

def process_all_samples_orthologs(alignment_dir, classification, flank_length, verbose=False):
    """Process all-samples ortholog groups for fixation-based diversity analysis"""
    results = []
    
    for ortholog_id, info in classification.items():
        category = info["category"]
        cross_group_status = info.get("cross_group_status", "NA")
        within_group_status = info.get("within_group_status", "")
        group1_present_count = info["group1_present_count"]
        group1_absent_count = info["group1_absent_count"]
        group2_present_count = info["group2_present_count"]
        group2_absent_count = info["group2_absent_count"]
        
        # Initialize side-specific results
        side_results = {}
        
        for side in ["left", "right"]:
            # Get alignment files
            group1_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.{category}.group1.{side}_flank_{flank_length}bp.mafft.fa")
            group2_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.{category}.group2.{side}_flank_{flank_length}bp.mafft.fa")
            combined_file = os.path.join(alignment_dir, f"{flank_length}bp", f"{ortholog_id}.{category}.combined.{side}_flank_{flank_length}bp.mafft.fa")
            
            # Check file existence
            group1_exists = os.path.exists(group1_file)
            group2_exists = os.path.exists(group2_file)
            combined_exists = os.path.exists(combined_file)
            
            if not group1_exists and not group2_exists:
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: No alignment files found")
                continue
            
            try:
                # Initialize metrics
                pi_group1 = None
                pi_group2 = None  
                dxy_group1_group2 = None
                
                # Load sequences
                group1_seqs = []
                group2_seqs = []
                
                if group1_exists:
                    group1_seqs = [str(record.seq) for record in SeqIO.parse(group1_file, "fasta")]
                
                if group2_exists:
                    group2_seqs = [str(record.seq) for record in SeqIO.parse(group2_file, "fasta")]
                
                # Calculate pi within Group1 (if >1 sequence)
                if len(group1_seqs) > 1:
                    pi_group1 = calculate_group_pi_ros(group1_seqs)
                
                # Calculate pi within Group2 (if >1 sequence)
                if len(group2_seqs) > 1:
                    pi_group2 = calculate_group_pi_ros(group2_seqs)
                
                # Calculate dXY between Group1 and Group2
                if group1_seqs and group2_seqs and combined_exists:
                    # Use combined alignment for accurate between-group calculations
                    combined_records = list(SeqIO.parse(combined_file, "fasta"))
                    
                    # Get sequence IDs to separate groups
                    group1_ids = set(record.id for record in SeqIO.parse(group1_file, "fasta")) if group1_exists else set()
                    
                    # Separate sequences from combined alignment
                    group1_combined_seqs = [str(record.seq) for record in combined_records if record.id in group1_ids]
                    group2_combined_seqs = [str(record.seq) for record in combined_records if record.id not in group1_ids]
                    
                    if group1_combined_seqs and group2_combined_seqs:
                        dxy_group1_group2 = calculate_between_groups_dxy_ros(group1_combined_seqs, group2_combined_seqs)
                
                elif group1_seqs and group2_seqs:
                    # Fallback: direct calculation without combined alignment
                    dxy_group1_group2 = calculate_between_groups_dxy_ros(group1_seqs, group2_seqs)
                
                # Store results for this side
                side_results[side] = {
                    "pi_group1": pi_group1,
                    "pi_group2": pi_group2,
                    "dxy_group1_group2": dxy_group1_group2
                }
                
            except Exception as e:
                if verbose:
                    print(f"Error processing {ortholog_id} {side}: {e}")
                continue
        
        # Calculate average metrics across sides if we have results
        if side_results:
            pi_group1_values = [data["pi_group1"] for data in side_results.values() if data["pi_group1"] is not None]
            pi_group2_values = [data["pi_group2"] for data in side_results.values() if data["pi_group2"] is not None]
            dxy_values = [data["dxy_group1_group2"] for data in side_results.values() if data["dxy_group1_group2"] is not None]
            
            results.append({
                "ortholog_id": ortholog_id,
                "category": category,
                "cross_group_status": cross_group_status,
                "within_group_status": within_group_status,
                "group1_present_count": group1_present_count,
                "group1_absent_count": group1_absent_count,
                "group2_present_count": group2_present_count,
                "group2_absent_count": group2_absent_count,
                "pi_group1": np.mean(pi_group1_values) if pi_group1_values else None,
                "pi_group2": np.mean(pi_group2_values) if pi_group2_values else None,
                "dxy_group1_group2": np.mean(dxy_values) if dxy_values else None,
                "pi_group1_std": np.std(pi_group1_values) if len(pi_group1_values) > 1 else None,
                "pi_group2_std": np.std(pi_group2_values) if len(pi_group2_values) > 1 else None,
                "dxy_group1_group2_std": np.std(dxy_values) if len(dxy_values) > 1 else None,
                "sides_processed": ",".join(side_results.keys()),
                "pi_group1_group2_diff": (np.mean(pi_group1_values) - np.mean(pi_group2_values)) if pi_group1_values and pi_group2_values else None,
                "flank_length": flank_length
            })
    
    return results

def main():
    args = parse_arguments()
    
    # Load all-samples classification
    print(f"Loading all-samples classification from {args.classification}")
    with open(args.classification, "r") as f:
        classification = json.load(f)
    
    print(f"Processing {len(classification)} classified all-samples ortholog groups using ratio of sums approach")
    
    # Process ortholog groups
    results = process_all_samples_orthologs(args.alignment_dir, classification, args.flank_length, args.verbose)
    
    print(f"Processed {len(results)} ortholog groups")
    
    # Print summary statistics
    cat_counts = defaultdict(int) 
    for result in results:
        cat_counts[result["category"]] += 1
    
    print("\nFixation pattern distribution of processed orthologs:")
    for category in sorted(cat_counts.keys()):
        print(f"  {category}: {cat_counts[category]} ortholog groups")
    
    # Convert results to DataFrame and save
    df = pd.DataFrame(results)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    df.to_csv(args.output, sep="\t", index=False)
    
    print(f"Saved diversity metrics (ratio of sums) for {len(results)} ortholog groups to {args.output}")
    print("Done!")

if __name__ == "__main__":
    main()
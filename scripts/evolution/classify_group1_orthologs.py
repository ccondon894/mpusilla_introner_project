#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
import json
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify Group1 ortholog groups by frequency for bottleneck model analysis')
    parser.add_argument('genotype_matrix', help='Path to genotype matrix file')
    parser.add_argument('fasta_dir', help='Directory containing consensus FASTA files')
    parser.add_argument('output_dir', help='Output directory for prepared FASTA files')
    parser.add_argument('--group1', help='Comma-separated list of samples in group 1')
    return parser.parse_args()

def load_fasta_sequences(fasta_dir, flank_length):
    """Load all FASTA sequences from the directory for a specific flanking length"""
    fasta_dict = {}
    pattern = f"{fasta_dir}/*flank_{flank_length}bp.consensus.fa"
    
    for file in glob.glob(pattern):
        with open(file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                fasta_dict[record.id] = str(record.seq)
    return fasta_dict

ACCEPTED_WITHIN_STATUS = {'consistent', 'singleton'}


def classify_group1_orthologs(df, group1):
    """
    Classify Group1 ortholog groups by allele frequency.

    Applies within-group classification filter: only groups where
    within_group_status is in ACCEPTED_WITHIN_STATUS are considered. Drops
    low_identity groups (suspicious ortholog grouping / partial deletions).
    Singletons are kept because they are legitimate freq=1 Group1
    polymorphisms — the downstream diversity computation handles the case
    where pi(present) can't be computed from a single sequence.
    """
    classification = {}
    skipped_low_identity = 0

    # Get unique ortholog IDs
    ortholog_ids = df['ortholog_id'].unique()

    for ortholog_id in ortholog_ids:
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]

        # Drop groups with suspect within-group classification
        if 'within_group_status' in ortholog_df.columns:
            within_vals = ortholog_df['within_group_status'].dropna().unique()
            if len(within_vals) > 0 and str(within_vals[0]) not in ACCEPTED_WITHIN_STATUS:
                skipped_low_identity += 1
                continue

        # Get Group1 samples only
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]

        # Skip if no Group1 data
        if g1_rows.empty:
            continue

        # Get presence status (1=present, 2=absent, 3=missing)
        g1_presence = g1_rows['presence'].to_list()

        # Skip if any missing data (3) - we need complete data for frequency analysis
        if 3 in g1_presence:
            continue

        # Count present and absent
        present_count = g1_presence.count(1)
        absent_count = g1_presence.count(2)
        total_count = present_count + absent_count

        # Skip if not all samples represented
        if total_count != len(group1):
            continue

        # Skip monomorphic cases (all present or all absent) - no polymorphism to analyze
        if present_count == 0 or absent_count == 0:
            continue

        # Calculate frequency (number of present alleles)
        frequency = present_count

        # Store classification information
        classification[ortholog_id] = {
            "frequency": frequency,
            "present_count": present_count,
            "absent_count": absent_count,
            "total_count": total_count,
            "frequency_category": get_frequency_category(frequency, total_count)
        }

    print(f"  Skipped {skipped_low_identity} groups with low_identity within-group status")
    return classification

def get_frequency_category(present_count, total_count):
    """Get descriptive frequency category"""
    if present_count == 1:
        return "singleton_present"
    elif present_count == total_count - 1:
        return "singleton_absent"
    elif present_count == 2:
        return "doubleton_present"
    elif present_count == total_count - 2:
        return "doubleton_absent"
    elif present_count <= total_count // 3:
        return "low_frequency"
    elif present_count >= 2 * total_count // 3:
        return "high_frequency"
    else:
        return "intermediate_frequency"

def prepare_group1_fasta_files(df, fasta_dict, classification, output_dir, group1, flank_length):
    """Prepare FASTA files for each ortholog group based on frequency classification"""
    
    subdir = os.path.join(output_dir, f"{flank_length}bp")
    os.makedirs(subdir, exist_ok=True)
    
    for ortholog_id, info in classification.items():
        frequency = info["frequency"]
        
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]
        
        # Process for both left and right flanking regions
        for side in ["left", "right"]:
            # Create FASTA for present sequences
            present_file = os.path.join(subdir, f"{ortholog_id}.freq_{frequency}.present.{side}_flank_{flank_length}bp.fa")
            present_sequences = []
            
            for _, row in g1_rows[g1_rows['presence'] == 1].iterrows():
                flank_id = get_flank_id(row, side, flank_length)
                if flank_id and flank_id in fasta_dict:
                    present_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            
            if present_sequences:
                with open(present_file, 'w') as f:
                    for seq_id, seq in present_sequences:
                        f.write(f">{seq_id}\n{seq}\n")
            
            # Create FASTA for absent sequences  
            absent_file = os.path.join(subdir, f"{ortholog_id}.freq_{frequency}.absent.{side}_flank_{flank_length}bp.fa")
            absent_sequences = []
            
            for _, row in g1_rows[g1_rows['presence'] == 2].iterrows():
                # For absent sequences, check if they have flanking sequences
                if row['start'] != -1 and row['end'] != -1:
                    flank_id = get_flank_id(row, side, flank_length)
                    if flank_id and flank_id in fasta_dict:
                        absent_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            
            if absent_sequences:
                with open(absent_file, 'w') as f:
                    for seq_id, seq in absent_sequences:
                        f.write(f">{seq_id}\n{seq}\n")

def get_flank_id(row, side, flank_length):
    """Get the flank ID for a given row, side, and flank length"""
    if row['start'] == -1 or row['end'] == -1:
        return None
        
    start = int(row['start'])
    end = int(row['end'])
    contig = row['contig']
    orientation = row['orientation']
    
    # Handle reverse orientation by swapping biological left/right flanks
    if orientation == 'reverse':
        if side == "left":  # Biological 5' flank extends upstream from sequence end
            flank_start, flank_end = end, end + flank_length
        else:  # Biological 3' flank extends downstream from sequence start
            flank_start, flank_end = start - flank_length, start
    else:  # Forward orientation
        if side == "left":  # Biological 5' flank extends upstream from sequence start
            flank_start, flank_end = start - flank_length, start
        else:  # Biological 3' flank extends downstream from sequence end
            flank_start, flank_end = end, end + flank_length
        
    return f"{contig}:{flank_start}-{flank_end}"

def main():
    args = parse_arguments()
    
    # Parse group1 list
    group1 = args.group1.split(',')
    
    print(f"Reading genotype matrix from {args.genotype_matrix}")
    df = pd.read_csv(args.genotype_matrix, sep="\t")
    
    # Process both flanking lengths
    flank_lengths = [100, 200]
    
    for flank_length in flank_lengths:
        print(f"\nProcessing {flank_length}bp flanking sequences...")
        
        print(f"Loading FASTA sequences for {flank_length}bp flanks")
        fasta_dict = load_fasta_sequences(args.fasta_dir, flank_length)
        
        print("Classifying Group1 ortholog groups by frequency")
        classification = classify_group1_orthologs(df, group1)
        
        print(f"Found {len(classification)} polymorphic ortholog groups in Group1")
        
        # Print frequency distribution
        freq_dist = defaultdict(int)
        for info in classification.values():
            freq_dist[info["frequency"]] += 1
        
        print("Frequency distribution:")
        for freq in sorted(freq_dist.keys()):
            print(f"  {freq}/11 present: {freq_dist[freq]} ortholog groups")
        
        print(f"Preparing FASTA files for {len(classification)} ortholog groups")
        prepare_group1_fasta_files(df, fasta_dict, classification, args.output_dir, group1, flank_length)
        
        # Save classification dictionary
        classification_file = os.path.join(args.output_dir, f"group1_classification_{flank_length}bp.json")
        with open(classification_file, "w") as f:
            json.dump(classification, f, indent=2)
        
        print(f"Classification saved to {classification_file}")
    
    print("Done!")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
import json
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify all-samples ortholog groups by fixation patterns between Group1 and Group2')
    parser.add_argument('genotype_matrix', help='Path to genotype matrix file')
    parser.add_argument('fasta_dir', help='Directory containing consensus FASTA files')
    parser.add_argument('output_dir', help='Output directory for prepared FASTA files')
    parser.add_argument('--group1', help='Comma-separated list of samples in group 1')
    parser.add_argument('--group2', help='Comma-separated list of samples in group 2')
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

ANCESTRAL_CROSS_GROUP = {
    'ancestral', 'likely_ancestral', 'ancestral_low_identity',
}


def _cross_group_value(ortholog_df):
    """Get the cross_group_status for this ortholog (first non-null row)."""
    if 'cross_group_status' not in ortholog_df.columns:
        return ''
    vals = ortholog_df['cross_group_status'].dropna().unique()
    if len(vals) == 0:
        return ''
    v = vals[0]
    if pd.isna(v):
        return ''
    return str(v)


def _within_group_value(ortholog_df):
    """Get the within_group_status for this ortholog (first non-null row)."""
    if 'within_group_status' not in ortholog_df.columns:
        return ''
    vals = ortholog_df['within_group_status'].dropna().unique()
    if len(vals) == 0:
        return ''
    return str(vals[0])


def classify_all_samples_orthologs(df, group1, group2):
    """
    Classify ortholog groups by fixation patterns between Group1 and Group2.

    Applies ancestry-aware filtering:
    - Drops any group where within_group_status == 'low_identity' (suspicious
      ortholog grouping, likely paralog mismerge or partial deletion).
    - For the group1_fixed_group2_fixed category, additionally requires
      cross_group_status in ANCESTRAL_CROSS_GROUP so Dxy is computed only on
      orthologs with evidence of shared ancestry (not convergent insertions
      at the same locus).
    """
    classification = {}
    skipped_low_identity = 0
    skipped_non_ancestral = 0

    # Get unique ortholog IDs
    ortholog_ids = df['ortholog_id'].unique()

    for ortholog_id in ortholog_ids:
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]

        # Drop low_identity groups entirely (suspicious ortholog grouping)
        within_status = _within_group_value(ortholog_df)
        if within_status == 'low_identity':
            skipped_low_identity += 1
            continue

        # Get Group1 and Group2 samples
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]
        g2_rows = ortholog_df[ortholog_df['sample'].isin(group2)]

        # Skip if no data for either group
        if g1_rows.empty or g2_rows.empty:
            continue

        # Get presence status (1=present, 2=absent, 3=missing)
        g1_presence = g1_rows['presence'].to_list()
        g2_presence = g2_rows['presence'].to_list()

        # Skip if any missing data (3) - we need complete data for fixation analysis
        if 3 in g1_presence or 3 in g2_presence:
            continue

        # Skip if not all samples represented in each group
        if len(g1_presence) != len(group1) or len(g2_presence) != len(group2):
            continue

        # Count present and absent in each group
        g1_present_count = g1_presence.count(1)
        g1_absent_count = g1_presence.count(2)
        g2_present_count = g2_presence.count(1)
        g2_absent_count = g2_presence.count(2)

        # Determine fixation category based on strict fixation patterns
        category = None

        if g1_present_count == len(group1) and g2_absent_count == len(group2):
            category = "group1_fixed_group2_absent"
        elif g1_absent_count == len(group1) and g2_present_count == len(group2):
            category = "group1_absent_group2_fixed"
        elif g1_present_count == len(group1) and g2_present_count == len(group2):
            category = "group1_fixed_group2_fixed"
        else:
            # Skip cases that don't match strict fixation patterns
            continue

        # For fixed_fixed (both clades have the introner), require cross-group
        # evidence of shared ancestry. Clade-specific categories don't have
        # cross-group members so this filter doesn't apply to them.
        cross_status = _cross_group_value(ortholog_df)
        if category == "group1_fixed_group2_fixed":
            if cross_status not in ANCESTRAL_CROSS_GROUP:
                skipped_non_ancestral += 1
                continue
        else:
            # Clade-specific: no cross-group comparison possible
            cross_status = 'NA'

        # Store classification information
        classification[ortholog_id] = {
            "category": category,
            "cross_group_status": cross_status,
            "within_group_status": within_status,
            "group1_present_count": g1_present_count,
            "group1_absent_count": g1_absent_count,
            "group2_present_count": g2_present_count,
            "group2_absent_count": g2_absent_count,
            "group1_total": len(group1),
            "group2_total": len(group2)
        }

    print(f"  Skipped {skipped_low_identity} low_identity groups")
    print(f"  Skipped {skipped_non_ancestral} non-ancestral fixed_fixed groups")
    return classification

def prepare_all_samples_fasta_files(df, fasta_dict, classification, output_dir, group1, group2, flank_length):
    """Prepare FASTA files for each ortholog group based on fixation classification"""
    
    subdir = os.path.join(output_dir, f"{flank_length}bp")
    os.makedirs(subdir, exist_ok=True)
    
    for ortholog_id, info in classification.items():
        category = info["category"]
        
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]
        g2_rows = ortholog_df[ortholog_df['sample'].isin(group2)]
        
        # Process for both left and right flanking regions
        for side in ["left", "right"]:
            # Create FASTA for Group1 sequences
            group1_file = os.path.join(subdir, f"{ortholog_id}.{category}.group1.{side}_flank_{flank_length}bp.fa")
            group1_sequences = []
            
            # For Group1, get sequences based on category
            if category == "group1_fixed_group2_absent" or category == "group1_fixed_group2_fixed":
                # Group1 has present sequences
                for _, row in g1_rows[g1_rows['presence'] == 1].iterrows():
                    flank_id = get_flank_id(row, side, flank_length)
                    if flank_id and flank_id in fasta_dict:
                        group1_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            elif category == "group1_absent_group2_fixed":
                # Group1 has absent sequences (flanking regions where introner would be)
                for _, row in g1_rows[g1_rows['presence'] == 2].iterrows():
                    # For absent sequences, check if they have flanking sequences
                    if row['start'] != -1 and row['end'] != -1:
                        flank_id = get_flank_id(row, side, flank_length)
                        if flank_id and flank_id in fasta_dict:
                            group1_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            
            if group1_sequences:
                with open(group1_file, 'w') as f:
                    for seq_id, seq in group1_sequences:
                        f.write(f">{seq_id}\n{seq}\n")
            
            # Create FASTA for Group2 sequences  
            group2_file = os.path.join(subdir, f"{ortholog_id}.{category}.group2.{side}_flank_{flank_length}bp.fa")
            group2_sequences = []
            
            # For Group2, get sequences based on category
            if category == "group1_absent_group2_fixed" or category == "group1_fixed_group2_fixed":
                # Group2 has present sequences
                for _, row in g2_rows[g2_rows['presence'] == 1].iterrows():
                    flank_id = get_flank_id(row, side, flank_length)
                    if flank_id and flank_id in fasta_dict:
                        group2_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            elif category == "group1_fixed_group2_absent":
                # Group2 has absent sequences (flanking regions where introner would be)
                for _, row in g2_rows[g2_rows['presence'] == 2].iterrows():
                    # For absent sequences, check if they have flanking sequences
                    if row['start'] != -1 and row['end'] != -1:
                        flank_id = get_flank_id(row, side, flank_length)
                        if flank_id and flank_id in fasta_dict:
                            group2_sequences.append((f"{row['sample']}_{flank_id}", fasta_dict[flank_id]))
            
            if group2_sequences:
                with open(group2_file, 'w') as f:
                    for seq_id, seq in group2_sequences:
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
    
    # Parse group lists
    group1 = args.group1.split(',')
    group2 = args.group2.split(',')
    
    print(f"Reading genotype matrix from {args.genotype_matrix}")
    df = pd.read_csv(args.genotype_matrix, sep="\t")
    
    # Process both flanking lengths
    flank_lengths = [100, 200]
    
    for flank_length in flank_lengths:
        print(f"\nProcessing {flank_length}bp flanking sequences...")
        
        print(f"Loading FASTA sequences for {flank_length}bp flanks")
        fasta_dict = load_fasta_sequences(args.fasta_dir, flank_length)
        
        print("Classifying all-samples ortholog groups by fixation patterns")
        classification = classify_all_samples_orthologs(df, group1, group2)
        
        print(f"Found {len(classification)} ortholog groups with strict fixation patterns")
        
        # Print category distribution
        cat_dist = defaultdict(int)
        for info in classification.values():
            cat_dist[info["category"]] += 1
        
        print("Fixation pattern distribution:")
        for category in sorted(cat_dist.keys()):
            print(f"  {category}: {cat_dist[category]} ortholog groups")
        
        print(f"Preparing FASTA files for {len(classification)} ortholog groups")
        prepare_all_samples_fasta_files(df, fasta_dict, classification, args.output_dir, group1, group2, flank_length)
        
        # Save classification dictionary
        classification_file = os.path.join(args.output_dir, f"all_samples_classification_{flank_length}bp.json")
        with open(classification_file, "w") as f:
            json.dump(classification, f, indent=2)
        
        print(f"Classification saved to {classification_file}")
    
    print("Done!")

if __name__ == "__main__":
    main()
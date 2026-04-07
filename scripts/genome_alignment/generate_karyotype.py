#!/usr/bin/env python3
"""
Generate a circos karyotype file from two FASTA index (.fai) files.

Usage:
    python generate_karyotype.py --strain1-fai strain1.fasta.fai --strain2-fai strain2.fasta.fai 
                                 --strain1-name Strain1 --strain2-name Strain2 
                                 --output karyotype.txt
"""

import argparse
import pandas as pd
import os
import re
import sys

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate a circos karyotype file from FASTA index files')
    parser.add_argument('--strain1-fai', required=True, help='FASTA index file for strain 1')
    parser.add_argument('--strain2-fai', required=True, help='FASTA index file for strain 2')
    parser.add_argument('--strain1-name', default='Strain1', help='Name of strain 1')
    parser.add_argument('--strain2-name', default='Strain2', help='Name of strain 2')
    parser.add_argument('--output', default='karyotype.txt', help='Output karyotype file')
    parser.add_argument('--colors', help='File with predefined colors (optional)')
    parser.add_argument('--min-length', type=int, default=0, 
                        help='Minimum chromosome length to include (default: 0)')
    return parser.parse_args()

def parse_fai(fai_file, min_length=0):
    """Parse a FASTA index file to get chromosome names and lengths."""
    if not os.path.isfile(fai_file):
        print(f"Error: {fai_file} not found.", file=sys.stderr)
        sys.exit(1)
    
    # Read the fai file
    columns = ['name', 'length', 'offset', 'linebases', 'linewidth']
    try:
        fai_data = pd.read_csv(fai_file, sep='\t', header=None, names=columns)
    except pd.errors.EmptyDataError:
        print(f"Error: {fai_file} is empty.", file=sys.stderr)
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Error: {fai_file} has an invalid format.", file=sys.stderr)
        sys.exit(1)
    
    # Filter by minimum length
    if min_length > 0:
        fai_data = fai_data[fai_data['length'] >= min_length]
    
    return fai_data[['name', 'length']]

def load_color_map(color_file):
    """Load predefined colors from a file."""
    color_map = {}
    if not color_file or not os.path.isfile(color_file):
        return color_map
    
    with open(color_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('=')
            if len(parts) == 2:
                chr_name = parts[0].strip()
                color = parts[1].strip()
                color_map[chr_name] = color
    
    return color_map

def generate_hue_colors(n, reverse=False):
    """Generate a list of n distinct hue colors for chromosomes.

    Uses hueNNN color names defined in circos.conf's <colors> block.
    For forward (strain1): evenly spaces n hue values across the color wheel
    and snaps each to the nearest defined hue name.
    For reverse (strain2): uses all defined hues in reverse order, cycling if needed.
    """
    # Defined hue values from circos.conf <colors> block
    DEFINED_HUE_VALUES = [
        0, 20, 32, 40, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192,
        208, 224, 240, 256, 272, 288, 304, 308, 312, 316, 320, 324, 328, 332, 340
    ]

    def snap_to_nearest(value):
        return min(DEFINED_HUE_VALUES, key=lambda x: abs(x - value))

    def hue_to_name(value):
        return f"hue{value:03d}"

    if reverse:
        # Use all defined hues in reverse order, cycling if needed
        reversed_hues = list(reversed(DEFINED_HUE_VALUES))
        return [hue_to_name(reversed_hues[i % len(reversed_hues)]) for i in range(n)]
    else:
        # Evenly space n values across 0-320 and snap to nearest defined hue
        if n == 1:
            return [hue_to_name(0)]
        max_hue = 320
        colors = []
        for i in range(n):
            hue_value = i * max_hue / (n - 1)
            snapped = snap_to_nearest(hue_value)
            colors.append(hue_to_name(snapped))
        return colors


def make_label(strain_name, chr_name):
    """Generate a label like 'CCMP1545-Chr1' from the chromosome name.

    Extracts the trailing number from scaffold/contig names.
    """
    # Try to extract trailing number from the chromosome name
    match = re.search(r'(\d+)$', chr_name)
    if match:
        num = match.group(1)
        return f"{strain_name}-Chr{num}"
    else:
        # Fallback: use the name with '#' replaced
        return chr_name.replace('#', '_')

def sort_key_by_number(name):
    """Extract trailing number from chromosome name for numeric sorting."""
    match = re.search(r'(\d+)$', name)
    return int(match.group(1)) if match else 0


def write_karyotype(strain1_data, strain2_data, strain1_name, strain2_name, color_map, output_file):
    """Write the karyotype file."""
    # Sort both strains by numeric suffix to ensure consistent color assignment
    # (scaffold_1 before scaffold_2, contig_1 before contig_2, etc.)
    strain1_sorted = strain1_data.copy()
    strain1_sorted['sort_key'] = strain1_sorted['name'].apply(sort_key_by_number)
    strain1_sorted = strain1_sorted.sort_values('sort_key').drop(columns='sort_key')

    strain2_sorted = strain2_data.copy()
    strain2_sorted['sort_key'] = strain2_sorted['name'].apply(sort_key_by_number)
    strain2_sorted = strain2_sorted.sort_values('sort_key').drop(columns='sort_key')

    # Generate hue colors: forward for strain1, reversed for strain2
    strain1_colors = generate_hue_colors(len(strain1_sorted), reverse=False)
    strain2_colors = generate_hue_colors(len(strain2_sorted), reverse=True)

    with open(output_file, 'w') as f:
        # Write strain1 chromosomes
        for i, (_, row) in enumerate(strain1_sorted.iterrows()):
            chr_name = row['name']

            # Build chr_id: ensure strain prefix and replace '#' with '_'
            if chr_name.startswith(f"{strain1_name}_") or chr_name.startswith(f"{strain1_name}#"):
                chr_id = chr_name.replace('#', '_')
            else:
                chr_id = f"{strain1_name}_{chr_name}".replace('#', '_')

            chr_label = make_label(strain1_name, chr_name)
            start = 0
            end = row['length']

            if chr_name in color_map:
                color = color_map[chr_name]
            else:
                color = strain1_colors[i]

            f.write(f"chr - {chr_id} {chr_label} {start} {end} {color}\n")

        # Write strain2 chromosomes
        for i, (_, row) in enumerate(strain2_sorted.iterrows()):
            chr_name = row['name']

            if chr_name.startswith(f"{strain2_name}_") or chr_name.startswith(f"{strain2_name}#"):
                chr_id = chr_name.replace('#', '_')
            else:
                chr_id = f"{strain2_name}_{chr_name}".replace('#', '_')

            chr_label = make_label(strain2_name, chr_name)
            start = 0
            end = row['length']

            if chr_name in color_map:
                color = color_map[chr_name]
            else:
                color = strain2_colors[i]

            f.write(f"chr - {chr_id} {chr_label} {start} {end} {color}\n")

def main():
    args = parse_args()
    
    print(f"Parsing FASTA index for {args.strain1_name}...")
    strain1_data = parse_fai(args.strain1_fai, args.min_length)
    print(f"Found {len(strain1_data)} chromosomes in {args.strain1_name}.")
    
    print(f"Parsing FASTA index for {args.strain2_name}...")
    strain2_data = parse_fai(args.strain2_fai, args.min_length)
    print(f"Found {len(strain2_data)} chromosomes in {args.strain2_name}.")
    
    # Load color map if provided
    color_map = load_color_map(args.colors)
    
    # Write the karyotype file
    print(f"Generating karyotype file {args.output}...")
    write_karyotype(strain1_data, strain2_data, args.strain1_name, args.strain2_name, color_map, args.output)
    print("Done!")
    
    # Print a reminder about the order, checking for existing strain prefixes
    print("\nReminder for your circos.conf file:")
    strain1_chroms = []
    for name in strain1_data['name']:
        if name.startswith(f"{args.strain1_name}_") or name.startswith(f"{args.strain1_name}#"):
            strain1_chroms.append(name)
        else:
            strain1_chroms.append(f"{args.strain1_name}_{name}")
            
    strain2_chroms = []
    for name in strain2_data['name']:
        if name.startswith(f"{args.strain2_name}_") or name.startswith(f"{args.strain2_name}#"):
            strain2_chroms.append(name)
        else:
            strain2_chroms.append(f"{args.strain2_name}_{name}")
    
    print(f"chromosomes_order = {','.join(strain1_chroms)},{','.join(strain2_chroms)}")
    print(f"chromosomes_reverse = {','.join(strain2_chroms)}")

if __name__ == "__main__":
    main()
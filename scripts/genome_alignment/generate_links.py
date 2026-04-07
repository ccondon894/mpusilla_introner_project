#!/usr/bin/env python3
"""
Generate a circos links file from a MUMmer coords file.

Usage:
    python generate_links.py --coords comparison_output.coords 
                             --strain1-name Strain1 --strain2-name Strain2 
                             --output links.txt
"""

import argparse
import sys
import os
import re
import math

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate a circos links file from MUMmer coords')
    parser.add_argument('--coords', required=True, help='MUMmer coords file')
    parser.add_argument('--strain1-name', default='Strain1', help='Name of strain 1')
    parser.add_argument('--strain2-name', default='Strain2', help='Name of strain 2')
    parser.add_argument('--output', default='links.txt', help='Output links file')
    parser.add_argument('--min-length', type=int, default=1000, 
                        help='Minimum alignment length to include (default: 1000bp)')
    parser.add_argument('--min-identity', type=float, default=80.0,
                        help='Minimum percent identity (default: 80.0)')
    parser.add_argument('--highlight', nargs='+', default=[],
                        help='Highlight specific chromosomes (space-separated)')
    parser.add_argument('--highlight-color', default='red',
                        help='Color for highlighted chromosomes (default: red)')
    return parser.parse_args()

def parse_coords_file(coords_file, min_length, min_identity):
    """Parse a MUMmer coords file."""
    alignments = []
    
    if not os.path.isfile(coords_file):
        print(f"Error: {coords_file} not found.", file=sys.stderr)
        sys.exit(1)
    
    with open(coords_file, 'r') as f:
        # Skip header lines
        line = f.readline()
        while line and not re.match(r'^\s*\d', line):
            line = f.readline()
        
        # Process alignment lines
        while line:
            # Handle different formats of show-coords output
            parts = line.strip().split()
            
            # Try to parse the line based on format
            try:
                if len(parts) >= 13:  # Standard format with tags
                    s1, e1, s2, e2 = map(int, parts[0:4])
                    len1, len2 = map(int, parts[4:6])
                    identity = float(parts[6])
                    ref_name = parts[11]
                    query_name = parts[12]
                elif len(parts) >= 10:  # Format without tags
                    s1, e1, s2, e2 = map(int, parts[0:4])
                    len1 = abs(e1 - s1) + 1
                    len2 = abs(e2 - s2) + 1
                    identity = float(parts[9])
                    ref_name = parts[7]
                    query_name = parts[8]
                else:
                    print(f"Warning: Couldn't parse line: {line.strip()}", file=sys.stderr)
                    line = f.readline()
                    continue
                
                # Filter by length and identity
                alignment_length = min(len1, len2)
                if alignment_length >= min_length and identity >= min_identity:
                    # Ensure start < end for circos
                    if s1 > e1:
                        s1, e1 = e1, s1
                    if s2 > e2:
                        s2, e2 = e2, s2
                    
                    alignments.append({
                        'ref_chr': ref_name,
                        'ref_start': s1,
                        'ref_end': e1,
                        'query_chr': query_name,
                        'query_start': s2,
                        'query_end': e2,
                        'identity': identity,
                        'length': alignment_length
                    })
            except (ValueError, IndexError) as e:
                print(f"Warning: Error parsing line: {line.strip()}", file=sys.stderr)
                print(f"Error details: {str(e)}", file=sys.stderr)
            
            line = f.readline()
    
    return alignments

def get_color_by_identity(identity, min_identity=80.0, max_identity=100.0):
    """Generate a color based on alignment identity."""
    # Normalize identity to a 0-1 scale
    normalized = (identity - min_identity) / (max_identity - min_identity)
    normalized = max(0, min(1, normalized))  # Clamp to 0-1 range
    
    # Create a blue-to-red gradient
    if normalized < 0.5:
        # Blue to purple
        r = int(255 * (normalized * 2))
        g = 0
        b = 255
    else:
        # Purple to red
        r = 255
        g = 0
        b = int(255 * (1 - (normalized - 0.5) * 2))
    
    return f"({r},{g},{b})"

def write_links_file(alignments, strain1_name, strain2_name, output_file, highlight=None, highlight_color='red'):
    """Write the links file for circos."""
    if highlight is None:
        highlight = []
    
    with open(output_file, 'w') as f:
        f.write("# Links file for circos\n")
        f.write("# Format: ID chr1 start1 end1 chr2 start2 end2 [options]\n\n")
        
        for i, aln in enumerate(alignments):
            # Convert chromosome names to circos format
            # Replace '#' with '_' for circos compatibility (# is comment char)
            # If strain name is already in the name, use as-is; otherwise prepend
            ref_name = aln['ref_chr'].replace('#', '_')
            query_name = aln['query_chr'].replace('#', '_')
            ref_chr = ref_name if strain1_name in ref_name else f"{strain1_name}_{ref_name}"
            query_chr = query_name if strain2_name in query_name else f"{strain2_name}_{query_name}"
            
            # Determine if this link should be highlighted
            should_highlight = (aln['ref_chr'] in highlight) or (aln['query_chr'] in highlight)
            
            # Set color based on highlight status or identity
            if should_highlight:
                color = highlight_color
            else:
                color = get_color_by_identity(aln['identity'])
            
            # Set thickness based on alignment length (optional)
            thickness = 1 + int(math.log10(aln['length']) - 2)
            thickness = max(1, min(5, thickness))  # Clamp to 1-5 range
            
            f.write(f"link{i} {ref_chr} {aln['ref_start']} {aln['ref_end']} "
                   f"{query_chr} {aln['query_start']} {aln['query_end']} "
                   f"color={color},thickness={thickness}p\n")

def main():
    args = parse_args()
    
    print(f"Parsing MUMmer coords file: {args.coords}")
    alignments = parse_coords_file(args.coords, args.min_length, args.min_identity)
    print(f"Found {len(alignments)} alignments meeting criteria (length >= {args.min_length}bp, identity >= {args.min_identity}%)")
    
    # Write the links file
    print(f"Generating links file: {args.output}")
    write_links_file(alignments, args.strain1_name, args.strain2_name, args.output, 
                    highlight=args.highlight, highlight_color=args.highlight_color)
    print("Done!")
    
    # Print summary of highlighted features if any
    if args.highlight:
        print(f"\nHighlighted chromosomes ({args.highlight_color}):")
        for chrom in args.highlight:
            print(f"  - {chrom}")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Generate a circos configuration file for comparing two genomes.

Usage:
    python generate_config.py --karyotype karyotype.txt --links links.txt
                             --strain1-name Strain1 --strain2-name Strain2
                             --output circos.conf
"""

import argparse
import os
import sys

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate a circos configuration file')
    parser.add_argument('--karyotype', required=True, help='Path to karyotype file')
    parser.add_argument('--links', required=True, help='Path to links file')
    parser.add_argument('--strain1-name', default='Strain1', help='Name of strain 1')
    parser.add_argument('--strain2-name', default='Strain2', help='Name of strain 2')
    parser.add_argument('--output', default='circos.conf', help='Output configuration file')
    parser.add_argument('--image-file', default='circos_comparison.png', help='Output image filename')
    parser.add_argument('--ticks', choices=['none', 'basic', 'detailed'], default='basic',
                        help='Tick style (none, basic, detailed)')
    parser.add_argument('--ideogram-thickness', type=int, default=30,
                        help='Ideogram thickness in pixels (default: 30)')
    parser.add_argument('--link-thickness', type=int, default=2,
                        help='Link thickness in pixels (default: 2)')
    return parser.parse_args()

def load_chromosome_ids(karyotype_file, strain1_name, strain2_name):
    """Load chromosome IDs from karyotype file."""
    strain1_chroms = []
    strain2_chroms = []
    
    if not os.path.isfile(karyotype_file):
        print(f"Error: {karyotype_file} not found.", file=sys.stderr)
        sys.exit(1)
    
    with open(karyotype_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 6 and parts[0] == 'chr' and parts[1] == '-':
                chr_id = parts[2]  # This is the ID exactly as it appears in the karyotype file
                
                # Don't try to prepend strain names - they should already be in the ID
                # Just check which strain this chromosome belongs to
                if strain1_name in chr_id and not chr_id.startswith(strain2_name):
                    strain1_chroms.append(chr_id)
                elif strain2_name in chr_id:
                    strain2_chroms.append(chr_id)
                # If neither strain name is in the ID, this is a problem
                elif not (strain1_name in chr_id or strain2_name in chr_id):
                    print(f"Warning: Chromosome ID {chr_id} doesn't contain either strain name.")
                    # Make a best guess based on position in the file
                    if len(strain1_chroms) == 0 or len(strain1_chroms) > len(strain2_chroms):
                        strain2_chroms.append(chr_id)
                    else:
                        strain1_chroms.append(chr_id)
    
    # Print what we found to help with debugging
    print(f"Found these chromosomes for {strain1_name}: {strain1_chroms}")
    print(f"Found these chromosomes for {strain2_name}: {strain2_chroms}")
    
    return strain1_chroms, strain2_chroms

def generate_config(karyotype_file, links_file, strain1_chroms, strain2_chroms, 
                   strain1_name, strain2_name, image_file, ticks_style, 
                   ideogram_thickness, link_thickness, output_file):
    """Generate the circos configuration file."""
    
    # Start with the basic configuration
    config = f"""# Circos configuration for genome comparison
# Comparing {strain1_name} vs {strain2_name}

# Include default configuration
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>

# Specify karyotype file
karyotype = {karyotype_file}

# Define ideogram
<ideogram>
    <spacing>
        default = 0.005r
    </spacing>
    
    # Ideogram position
    radius    = 0.9r
    thickness = {ideogram_thickness}p
    fill      = yes
    
    # Labels
    show_label     = yes
    label_font     = default
    label_radius   = dims(ideogram,radius_outer) + 0.05r
    label_size     = 30
    label_parallel = yes
</ideogram>

# Define how links are displayed
<links>
    <link>
        file          = {links_file}
        radius        = 0.6r
        bezier_radius = 0r
        thickness     = {link_thickness}p
        ribbon        = yes
    </link>
</links>
"""

    # Add tick configuration based on style
    if ticks_style != 'none':
        config += """
# Configure ticks
show_ticks          = yes
show_tick_labels    = yes

<ticks>
    radius           = dims(ideogram,radius_outer)
    color            = black
    thickness        = 2p
"""
        
        if ticks_style == 'basic':
            config += """    
    # 10kb ticks
    <tick>
        spacing        = 10u
        size           = 5p
        show_label     = no
    </tick>
    
    # 100kb ticks
    <tick>
        spacing        = 100u
        size           = 10p
        show_label     = yes
        label_size     = 20p
        format         = %d
    </tick>
</ticks>
"""
        elif ticks_style == 'detailed':
            config += """    
    # 10kb ticks
    <tick>
        spacing        = 10u
        size           = 5p
        show_label     = no
    </tick>
    
    # 100kb ticks
    <tick>
        spacing        = 100u
        size           = 8p
        show_label     = no
    </tick>
    
    # 1Mb ticks
    <tick>
        spacing        = 1000u
        size           = 12p
        show_label     = yes
        label_size     = 24p
        format         = %dMb
        suffix         = " "
    </tick>
</ticks>
"""

    # Add chromosome order configuration
    strain1_chroms_str = ';'.join(strain1_chroms)
    strain2_chroms_str = ';'.join(strain2_chroms)
    
    config += f"""
# Split the chromosomes into two halves of the circle
# {strain1_name} on top, {strain2_name} on bottom
chromosomes = {strain1_chroms_str};{strain2_chroms_str}
"""

    # Add image configuration
    config += f"""
# Image configuration
<image>
    dir   = .
    file  = {image_file}
    png   = yes
    svg   = yes
    radius = 1500p
    background = white
    # Enables antialiasing
    anti_aliasing = yes
</image>
"""


    # Write the configuration to file
    with open(output_file, 'w') as f:
        f.write(config)

def main():
    args = parse_args()
    
    # Load chromosome IDs from karyotype file
    print(f"Loading chromosome information from {args.karyotype}...")
    strain1_chroms, strain2_chroms = load_chromosome_ids(args.karyotype, args.strain1_name, args.strain2_name)
    print(f"Found {len(strain1_chroms)} chromosomes for {args.strain1_name} and {len(strain2_chroms)} for {args.strain2_name}.")
    
    # Generate circos configuration
    print(f"Generating circos configuration file: {args.output}")
    generate_config(args.karyotype, args.links, strain1_chroms, strain2_chroms,
                   args.strain1_name, args.strain2_name, args.image_file, args.ticks,
                   args.ideogram_thickness, args.link_thickness, args.output)
    print("Done!")
    
    # Print run command
    print("\nTo run circos with this configuration:")
    print(f"circos -conf {args.output}")

if __name__ == "__main__":
    main()
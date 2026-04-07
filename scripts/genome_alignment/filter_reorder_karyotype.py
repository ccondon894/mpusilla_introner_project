#!/usr/bin/env python3
"""
Filter and reorder RCC1749 chromosomes in karyotype file.

This script:
1. Removes specified small RCC1749 contigs
2. Reorders RCC1749-Chr20 to a specific position

Usage:
    python3 filter_reorder_karyotype.py --input input.kar --output output.kar \
                                        --remove-contigs 2,1,33,30,32,40,43,24,42 \
                                        --move-contig 20 --insert-after 12 --insert-before 25
"""

import argparse
import sys


def parse_karyotype_line(line):
    """Parse a karyotype line into components."""
    if not line.strip() or line.startswith('#'):
        return None

    parts = line.strip().split()
    if len(parts) >= 7 and parts[0] == 'chr' and parts[1] == '-':
        return {
            'full_line': line.rstrip('\n'),
            'chr_id': parts[2],  # e.g., RCC1749_0_intronerless_contig_20
            'label': parts[3],   # e.g., RCC1749-Chr20
            'parts': parts
        }
    return None


def extract_contig_number(chr_id):
    """Extract contig number from chromosome ID."""
    # Format: RCC1749_0_intronerless_contig_20 -> 20
    try:
        return int(chr_id.split('_')[-1])
    except (ValueError, IndexError):
        return None


def main():
    parser = argparse.ArgumentParser(description='Filter and reorder RCC1749 chromosomes in karyotype')
    parser.add_argument('--input', default='circos-files/CCMP1545_vs_RCC1749_reordered_flipped.kar',
                       help='Input karyotype file')
    parser.add_argument('--output', default='circos-files/CCMP1545_vs_RCC1749_reordered_flipped_filtered.kar',
                       help='Output karyotype file')
    parser.add_argument('--remove-contigs', default='2,1,33,30,32,40,43,24,42',
                       help='Comma-separated list of RCC1749 contig numbers to remove')
    parser.add_argument('--move-contig', type=int, default=20,
                       help='Contig number to move (default: 20)')
    parser.add_argument('--insert-after', type=int, default=12,
                       help='Insert after this contig number (default: 12)')
    parser.add_argument('--insert-before', type=int, default=25,
                       help='Insert before this contig number (default: 25)')
    parser.add_argument('--move-contig-color', default='',
                       help='Override color for the moved contig (e.g., hue080)')

    args = parser.parse_args()

    # Parse remove list
    remove_contigs = set(int(x.strip()) for x in args.remove_contigs.split(','))

    print(f"Reading karyotype from {args.input}")
    with open(args.input, 'r') as f:
        lines = f.readlines()

    # Separate CCMP1545 and RCC1749 entries
    header_lines = []
    ccmp_lines = []
    rcc_lines = []

    for line in lines:
        parsed = parse_karyotype_line(line)
        if parsed is None:
            # Header or comment
            header_lines.append(line)
        elif 'CCMP1545' in parsed['chr_id']:
            ccmp_lines.append(parsed)
        elif 'RCC1749' in parsed['chr_id']:
            rcc_lines.append(parsed)

    # Filter RCC1749 lines
    print(f"\nFiltering RCC1749 contigs...")
    print(f"Removing contigs: {', '.join(map(str, sorted(remove_contigs)))}")
    filtered_rcc = []
    removed_count = 0

    for entry in rcc_lines:
        contig_num = extract_contig_number(entry['chr_id'])
        if contig_num in remove_contigs:
            print(f"  Removed: {entry['chr_id']} ({entry['label']})")
            removed_count += 1
        else:
            filtered_rcc.append(entry)

    print(f"Removed {removed_count} contigs")

    # Reorder: move specified contig to new position
    print(f"\nReordering RCC1749 chromosomes...")
    print(f"Moving contig {args.move_contig} to be between contig {args.insert_after} and contig {args.insert_before}")

    move_entry = None
    reordered_rcc = []

    # First pass: find and remove the contig to move
    for entry in filtered_rcc:
        contig_num = extract_contig_number(entry['chr_id'])
        if contig_num == args.move_contig:
            move_entry = entry
            # Apply color override if specified
            if args.move_contig_color:
                old_parts = entry['full_line'].split()
                old_parts[6] = args.move_contig_color
                move_entry['full_line'] = ' '.join(old_parts)
                print(f"  Found contig {args.move_contig}, recolored to {args.move_contig_color}")
            else:
                print(f"  Found contig {args.move_contig} at position {len(reordered_rcc)}")
        else:
            reordered_rcc.append(entry)

    # Second pass: insert the contig at the correct position
    final_rcc = []
    insert_position = None

    for i, entry in enumerate(reordered_rcc):
        contig_num = extract_contig_number(entry['chr_id'])

        # Insert before the insert_before contig
        if contig_num == args.insert_before and insert_position is None:
            if move_entry:
                final_rcc.append(move_entry)
                insert_position = len(final_rcc) - 1
                print(f"  Inserted contig {args.move_contig} at position {insert_position}")

        final_rcc.append(entry)

    # If we never found the insert_before contig, append at end
    if insert_position is None and move_entry:
        final_rcc.append(move_entry)
        insert_position = len(final_rcc) - 1
        print(f"  Inserted contig {args.move_contig} at end (position {insert_position})")

    # Write output
    print(f"\nWriting filtered and reordered karyotype to {args.output}")
    with open(args.output, 'w') as f:
        # Write headers
        for line in header_lines:
            f.write(line)

        # Write CCMP1545 chromosomes
        for entry in ccmp_lines:
            f.write(entry['full_line'] + '\n')

        # Write filtered/reordered RCC1749 chromosomes
        for entry in final_rcc:
            f.write(entry['full_line'] + '\n')

    print(f"\nSummary:")
    print(f"  CCMP1545 chromosomes: {len(ccmp_lines)}")
    print(f"  RCC1749 chromosomes (original): {len(rcc_lines)}")
    print(f"  RCC1749 chromosomes (filtered): {len(final_rcc)}")
    print(f"  Reduction: {len(rcc_lines) - len(final_rcc)} removed")


if __name__ == '__main__':
    main()

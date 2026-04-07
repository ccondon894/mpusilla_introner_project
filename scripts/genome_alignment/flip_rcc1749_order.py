#!/usr/bin/env python3
"""
Script to flip the order of RCC1749 chromosomes in the karyotype file
to correctly align with CCMP1545 chromosomes in the circos plot.
"""
import argparse

def flip_rcc1749_chromosomes(input_file, output_file):
    """
    Read the reordered karyotype file and flip the order of RCC1749 chromosomes
    while keeping CCMP1545 chromosomes in the same order.
    """
    # Read all lines
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    # Separate CCMP1545 and RCC1749 chromosomes
    ccmp1545_lines = []
    rcc1749_lines = []

    for line in lines:
        if 'CCMP1545' in line and line.startswith('chr'):
            ccmp1545_lines.append(line)
        elif 'RCC1749' in line and line.startswith('chr'):
            rcc1749_lines.append(line)

    # Keep CCMP1545 in original order, reverse RCC1749 order
    rcc1749_lines.reverse()

    # Write new karyotype file
    with open(output_file, 'w') as f:
        # Write CCMP1545 chromosomes first
        for line in ccmp1545_lines:
            f.write(line + '\n')

        # Write RCC1749 chromosomes in reversed order
        for line in rcc1749_lines:
            f.write(line + '\n')

    print(f"Created flipped karyotype file: {output_file}")
    print(f"CCMP1545 chromosomes: {len(ccmp1545_lines)} (original order)")
    print(f"RCC1749 chromosomes: {len(rcc1749_lines)} (reversed order)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Flip RCC1749 chromosome order in karyotype file")
    parser.add_argument('--input', required=True, help="Input karyotype file")
    parser.add_argument('--output', required=True, help="Output karyotype file")
    args = parser.parse_args()

    flip_rcc1749_chromosomes(args.input, args.output)
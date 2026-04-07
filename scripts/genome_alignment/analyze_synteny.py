#!/usr/bin/env python3
"""
Synteny Analysis Script for CCMP1545 vs RCC1749 Genome Comparison

This script analyzes syntenic relationships between CCMP1545 scaffolds and RCC1749 contigs
based on genome alignment data from MUMmer/NUCmer output converted to circos links format.
"""

import pandas as pd
import argparse
from collections import defaultdict
import sys

def parse_links_file(links_file):
    """Parse the circos links file and calculate synteny statistics."""
    print(f"Reading links file: {links_file}")
    
    # Read the links file (space-separated: linkID chr1 start1 end1 chr2 start2 end2 options)
    # Note: cannot use comment='#' because chromosome IDs contain '#' (vg_paths format)
    rows = []
    with open(links_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 7:
                rows.append(parts[:8] if len(parts) >= 8 else parts[:7] + [''])
    df = pd.DataFrame(rows, columns=['link_id', 'CCMP1545_scaffold', 'CCMP1545_start', 'CCMP1545_end',
                                     'RCC1749_contig', 'RCC1749_start', 'RCC1749_end', 'options'])
    df['CCMP1545_start'] = pd.to_numeric(df['CCMP1545_start'])
    df['CCMP1545_end'] = pd.to_numeric(df['CCMP1545_end'])
    df['RCC1749_start'] = pd.to_numeric(df['RCC1749_start'])
    df['RCC1749_end'] = pd.to_numeric(df['RCC1749_end'])
    
    print(f"Total syntenic blocks: {len(df):,}")
    
    # Calculate alignment lengths for each block
    df['CCMP1545_length'] = abs(df['CCMP1545_end'] - df['CCMP1545_start'])
    df['RCC1749_length'] = abs(df['RCC1749_end'] - df['RCC1749_start'])
    df['alignment_length'] = df[['CCMP1545_length', 'RCC1749_length']].min(axis=1)
    
    return df

def calculate_synteny_stats(df):
    """Calculate synteny statistics for each scaffold-contig pair."""
    synteny_stats = defaultdict(lambda: {'total_length': 0, 'block_count': 0})
    
    for _, row in df.iterrows():
        key = (row['CCMP1545_scaffold'], row['RCC1749_contig'])
        synteny_stats[key]['total_length'] += row['alignment_length']
        synteny_stats[key]['block_count'] += 1
    
    return synteny_stats

def find_best_syntenic_partners(synteny_stats):
    """Find the best syntenic partner for each CCMP1545 scaffold."""
    scaffold_best = {}
    
    # Group by scaffold and find the contig with maximum total alignment length
    scaffold_contigs = defaultdict(list)
    for (scaffold, contig), stats in synteny_stats.items():
        scaffold_contigs[scaffold].append((contig, stats['total_length'], stats['block_count']))
    
    for scaffold, contigs in scaffold_contigs.items():
        # Sort by total alignment length (descending)
        contigs.sort(key=lambda x: x[1], reverse=True)
        best_contig, best_length, best_blocks = contigs[0]
        scaffold_best[scaffold] = {
            'best_contig': best_contig,
            'total_length': best_length,
            'block_count': best_blocks,
            'all_partners': contigs
        }
    
    return scaffold_best

def parse_karyotype_file(karyotype_file):
    """Parse the karyotype file to get scaffold and contig information."""
    print(f"Reading karyotype file: {karyotype_file}")
    
    ccmp1545_scaffolds = []
    rcc1749_contigs = []
    
    with open(karyotype_file, 'r') as f:
        for line in f:
            if line.strip() and not line.strip().startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 7:
                    chr_id = parts[2]
                    chr_name = parts[3]
                    chr_length = int(parts[5])
                    color = parts[6]
                    
                    if 'CCMP1545' in chr_id:
                        scaffold_num = int(chr_id.split('_')[-1])
                        ccmp1545_scaffolds.append({
                            'scaffold_num': scaffold_num,
                            'scaffold_id': chr_id,
                            'display_name': chr_name,
                            'length': chr_length,
                            'color': color
                        })
                    elif 'RCC1749' in chr_id:
                        rcc1749_contigs.append({
                            'contig_id': chr_id,
                            'display_name': chr_name,
                            'length': chr_length,
                            'color': color
                        })
    
    # Sort CCMP1545 scaffolds by number
    ccmp1545_scaffolds.sort(key=lambda x: x['scaffold_num'])
    
    return ccmp1545_scaffolds, rcc1749_contigs

def generate_synteny_report(scaffold_best, ccmp1545_scaffolds):
    """Generate a comprehensive synteny report."""
    print("\n" + "="*80)
    print("SYNTENY ANALYSIS REPORT: CCMP1545 vs RCC1749")
    print("="*80)
    
    print(f"\nAnalyzing {len(ccmp1545_scaffolds)} CCMP1545 scaffolds:")
    print("-" * 80)
    print(f"{'Scaffold':<20} {'Best Partner':<25} {'Alignment (bp)':<15} {'Blocks':<8} {'% of Scaffold'}")
    print("-" * 80)
    
    total_aligned_bp = 0
    scaffolds_with_synteny = 0
    
    for scaffold_info in ccmp1545_scaffolds:
        scaffold_id = scaffold_info['scaffold_id']
        scaffold_length = scaffold_info['length']
        scaffold_display = f"Chr{scaffold_info['scaffold_num']}"
        
        if scaffold_id in scaffold_best:
            scaffolds_with_synteny += 1
            stats = scaffold_best[scaffold_id]
            contig_name = stats['best_contig'].replace('RCC1749_0_intronerless_contig_', 'Chr')
            alignment_bp = stats['total_length']
            blocks = stats['block_count']
            coverage_pct = (alignment_bp / scaffold_length) * 100
            total_aligned_bp += alignment_bp
            
            print(f"{scaffold_display:<20} {contig_name:<25} {alignment_bp:>12,} {blocks:>7} {coverage_pct:>8.1f}%")
        else:
            print(f"{scaffold_display:<20} {'No synteny':<25} {0:>12,} {0:>7} {0:>8.1f}%")
    
    print("-" * 80)
    print(f"Total scaffolds with synteny: {scaffolds_with_synteny}/{len(ccmp1545_scaffolds)}")
    print(f"Total aligned base pairs: {total_aligned_bp:,}")
    
    return total_aligned_bp, scaffolds_with_synteny

def create_reordered_karyotype(scaffold_best, ccmp1545_scaffolds, rcc1749_contigs, output_file):
    """Create a new karyotype file with RCC1749 contigs ordered by synteny."""
    print(f"\nCreating reordered karyotype file: {output_file}")
    
    # Create mapping of contig_id to contig info
    contig_lookup = {c['contig_id']: c for c in rcc1749_contigs}
    
    # Track which contigs have been used
    used_contigs = set()
    
    lines = []
    
    # Add CCMP1545 scaffolds (maintain original order)
    for scaffold_info in ccmp1545_scaffolds:
        line = f"chr - {scaffold_info['scaffold_id']} {scaffold_info['display_name']} 0 {scaffold_info['length']} {scaffold_info['color']}"
        lines.append(line)
    
    # Add RCC1749 contigs in synteny order, colored to match their CCMP1545 partner
    for scaffold_info in ccmp1545_scaffolds:
        scaffold_id = scaffold_info['scaffold_id']
        if scaffold_id in scaffold_best:
            best_contig_id = scaffold_best[scaffold_id]['best_contig']
            if best_contig_id in contig_lookup and best_contig_id not in used_contigs:
                contig_info = contig_lookup[best_contig_id]
                # Use the CCMP1545 scaffold's color so contig matches its syntenic partner
                partner_color = scaffold_info['color']
                line = f"chr - {contig_info['contig_id']} {contig_info['display_name']} 0 {contig_info['length']} {partner_color}"
                lines.append(line)
                used_contigs.add(best_contig_id)
    
    # Add remaining RCC1749 contigs (those without strong synteny) sorted by size
    remaining_contigs = [c for c in rcc1749_contigs if c['contig_id'] not in used_contigs]
    remaining_contigs.sort(key=lambda x: x['length'], reverse=True)
    
    for contig_info in remaining_contigs:
        line = f"chr - {contig_info['contig_id']} {contig_info['display_name']} 0 {contig_info['length']} {contig_info['color']}"
        lines.append(line)
    
    # Write the new karyotype file
    with open(output_file, 'w') as f:
        for line in lines:
            f.write(line + '\n')
    
    print(f"Reordered karyotype written with {len(lines)} chromosomes")
    print(f"  - {len(ccmp1545_scaffolds)} CCMP1545 scaffolds")
    print(f"  - {len(used_contigs)} RCC1749 contigs with synteny (ordered by partner)")
    print(f"  - {len(remaining_contigs)} RCC1749 contigs without strong synteny (ordered by size)")

def main():
    parser = argparse.ArgumentParser(description='Analyze synteny between CCMP1545 and RCC1749 genomes')
    parser.add_argument('--links', required=True, help='Path to circos links file (.tsv)')
    parser.add_argument('--karyotype', required=True, help='Path to circos karyotype file')
    parser.add_argument('--output-karyotype', help='Output path for reordered karyotype file')
    parser.add_argument('--detailed-report', help='Output path for detailed synteny report (CSV)')
    
    args = parser.parse_args()
    
    # Parse input files
    df = parse_links_file(args.links)
    synteny_stats = calculate_synteny_stats(df)
    scaffold_best = find_best_syntenic_partners(synteny_stats)
    ccmp1545_scaffolds, rcc1749_contigs = parse_karyotype_file(args.karyotype)
    
    # Generate report
    total_aligned, scaffolds_with_synteny = generate_synteny_report(scaffold_best, ccmp1545_scaffolds)
    
    # Create reordered karyotype if requested
    if args.output_karyotype:
        create_reordered_karyotype(scaffold_best, ccmp1545_scaffolds, rcc1749_contigs, args.output_karyotype)
    
    # Save detailed report if requested
    if args.detailed_report:
        detailed_data = []
        for scaffold_info in ccmp1545_scaffolds:
            scaffold_id = scaffold_info['scaffold_id']
            row = {
                'CCMP1545_Scaffold': f"Chr{scaffold_info['scaffold_num']}",
                'Scaffold_Length': scaffold_info['length'],
                'Best_RCC1749_Partner': 'No synteny',
                'Total_Alignment_Length': 0,
                'Number_of_Blocks': 0,
                'Coverage_Percent': 0.0
            }
            
            if scaffold_id in scaffold_best:
                stats = scaffold_best[scaffold_id]
                row.update({
                    'Best_RCC1749_Partner': stats['best_contig'].replace('RCC1749_0_intronerless_contig_', 'Chr'),
                    'Total_Alignment_Length': stats['total_length'],
                    'Number_of_Blocks': stats['block_count'],
                    'Coverage_Percent': (stats['total_length'] / scaffold_info['length']) * 100
                })
            
            detailed_data.append(row)
        
        pd.DataFrame(detailed_data).to_csv(args.detailed_report, index=False)
        print(f"\nDetailed report saved to: {args.detailed_report}")

if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
Generate ribbon-style links for circos visualization with colors matching the karyotype file.
This ensures ribbons use the same color scheme as their corresponding chromosomes.
"""

import pandas as pd
from collections import defaultdict
import argparse

def parse_karyotype_colors(karyotype_file):
    """Extract chromosome colors from karyotype file."""
    colors = {}
    with open(karyotype_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and line.startswith('chr'):
                parts = line.split()
                if len(parts) >= 7:
                    chr_id = parts[2]  # e.g., CCMP1545_0_scaffold_1
                    color = parts[6]   # e.g., hue000 (7th field, index 6)
                    colors[chr_id] = color
                    
    print(f"Loaded colors for {len(colors)} chromosomes from karyotype")
    # Debug: print first few color mappings
    sample_colors = dict(list(colors.items())[:3])
    print(f"Sample color mappings: {sample_colors}")
    return colors

def parse_links_file(links_file, min_link_length=1000):
    """Parse the original links file and filter by minimum length."""
    print(f"Reading links from {links_file}")
    
    all_links = []
    filtered_links = []
    short_links_count = 0
    
    with open(links_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            # Links format: linkID chr1 start1 end1 chr2 start2 end2 options
            # Skip the link ID prefix if present
            if parts[0].startswith('link'):
                parts = parts[1:]
            if len(parts) >= 6:
                chr1, start1, end1, chr2, start2, end2 = parts[:6]
                
                try:
                    start1, end1 = int(start1), int(end1)
                    orig_start2, orig_end2 = int(start2), int(end2)
                    start2, end2 = orig_start2, orig_end2

                    if start1 > end1:
                        start1, end1 = end1, start1
                    if start2 > end2:
                        start2, end2 = end2, start2

                    orientation = "same" if orig_start2 < orig_end2 else "opposite"
                    
                    link_length = end1 - start1
                    
                    link = {
                        'chr1': chr1,
                        'start1': start1,
                        'end1': end1,
                        'chr2': chr2,
                        'start2': start2,
                        'end2': end2,
                        'length': link_length,
                        'orientation': orientation
                    }
                    
                    all_links.append(link)
                    
                    # Filter by minimum length
                    if link_length >= min_link_length:
                        filtered_links.append(link)
                    else:
                        short_links_count += 1
                        
                except ValueError:
                    continue
    
    print(f"Parsed {len(all_links)} total valid links")
    print(f"Filtered out {short_links_count} links shorter than {min_link_length:,} bp")
    print(f"Retained {len(filtered_links)} links for ribbon generation")
    
    return filtered_links

def find_connected_components(links, merge_threshold):
    """
    Find connected components of links that should be merged.
    Two links are connected if they are within merge_threshold on both chromosomes.
    Returns a list of link groups (connected components).
    """
    from collections import defaultdict, deque
    
    # Build adjacency graph
    adjacency = defaultdict(list)
    n = len(links)
    
    for i in range(n):
        for j in range(i + 1, n):
            link1, link2 = links[i], links[j]
            
            # Only merge links from same chromosome pair
            if link1['chr1'] != link2['chr1'] or link1['chr2'] != link2['chr2']:
                continue
            
            # Check if orientations are compatible (same or can be reconciled)
            if link1['orientation'] != link2['orientation']:
                continue
            
            # Calculate gaps on both chromosomes
            gap_chr1 = min(abs(link1['start1'] - link2['end1']), abs(link1['end1'] - link2['start1']))
            gap_chr2 = min(abs(link1['start2'] - link2['end2']), abs(link1['end2'] - link2['start2']))
            
            # If both gaps are within threshold, links are adjacent
            if gap_chr1 <= merge_threshold and gap_chr2 <= merge_threshold:
                adjacency[i].append(j)
                adjacency[j].append(i)
    
    # Find connected components using BFS
    visited = [False] * n
    components = []
    
    for i in range(n):
        if not visited[i]:
            # Start new component
            component = []
            queue = deque([i])
            visited[i] = True
            
            while queue:
                node = queue.popleft()
                component.append(node)
                
                for neighbor in adjacency[node]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            
            components.append(component)
    
    return components

def merge_nearby_links(links, merge_threshold):
    """
    Merge links that are within merge_threshold of each other using connected components.
    Returns a list of merged links with statistics.
    """
    if not links:
        return []
    
    print(f"Starting link merging with threshold: {merge_threshold:,} bp")
    
    # Group links by chromosome pair first for efficiency
    chr_pair_groups = defaultdict(list)
    for i, link in enumerate(links):
        pair_key = (link['chr1'], link['chr2'])
        chr_pair_groups[pair_key].append((i, link))
    
    merged_links = []
    total_original = len(links)
    total_merged = 0
    merge_events = 0
    
    # Process each chromosome pair separately
    for pair_key, pair_items in chr_pair_groups.items():
        pair_links = [item[1] for item in pair_items]
        
        if len(pair_links) == 1:
            # Single link, no merging needed
            merged_links.extend(pair_links)
            continue
        
        # Find connected components for this chromosome pair
        components = find_connected_components(pair_links, merge_threshold)
        
        # Merge each connected component
        for component_indices in components:
            component_links = [pair_links[i] for i in component_indices]
            
            if len(component_links) == 1:
                # Single link in component
                merged_links.append(component_links[0])
            else:
                # Merge multiple links in component
                merge_events += 1
                merged_link = merge_link_group(component_links, len(component_links))
                merged_links.append(merged_link)
                
                # Debug info for large merges
                if len(component_links) > 3:
                    chr1, chr2 = pair_key
                    print(f"  Large merge: {chr1} ↔ {chr2}: {len(component_links)} links merged")
    
    total_merged = len(merged_links)
    reduction = total_original - total_merged
    
    print(f"Link merging summary:")
    print(f"  Original links: {total_original:,}")
    print(f"  Merged links: {total_merged:,}")
    print(f"  Links eliminated: {reduction:,}")
    print(f"  Merge events: {merge_events}")
    print(f"  Reduction factor: {total_original/max(total_merged,1):.1f}x")
    
    return merged_links

def merge_link_group(link_group, original_count):
    """
    Merge a group of links into a single representative link.
    Preserves metadata about the merge operation.
    """
    if not link_group:
        return None
    
    # Use bounds from all links in group
    min_start1 = min(link['start1'] for link in link_group)
    max_end1 = max(link['end1'] for link in link_group)
    min_start2 = min(link['start2'] for link in link_group)
    max_end2 = max(link['end2'] for link in link_group)
    
    # Calculate total length and dominant orientation
    total_length = sum(link['length'] for link in link_group)
    orientations = [link['orientation'] for link in link_group]
    same_count = sum(1 for o in orientations if o == 'same')
    dominant_orientation = 'same' if same_count >= len(orientations)/2 else 'opposite'
    
    # Create merged link using the first link as template
    merged_link = {
        'chr1': link_group[0]['chr1'],
        'start1': min_start1,
        'end1': max_end1,
        'chr2': link_group[0]['chr2'],
        'start2': min_start2,
        'end2': max_end2,
        'length': max_end1 - min_start1,  # New span length
        'orientation': dominant_orientation,
        'merged_from': original_count,  # Track how many links were merged
        'total_original_length': total_length  # Track original total length
    }
    
    return merged_link

def group_links_by_chromosome_pairs(links):
    """Group links by chromosome pair and calculate statistics."""
    pairs = defaultdict(list)
    
    for link in links:
        pair_key = (link['chr1'], link['chr2'])
        pairs[pair_key].append(link)
    
    pair_stats = {}
    for pair_key, pair_links in pairs.items():
        chr1, chr2 = pair_key
        
        min_start1 = min(link['start1'] for link in pair_links)
        max_end1 = max(link['end1'] for link in pair_links)
        min_start2 = min(link['start2'] for link in pair_links)
        max_end2 = max(link['end2'] for link in pair_links)
        
        orientations = [link['orientation'] for link in pair_links]
        same_count = sum(1 for o in orientations if o == 'same')
        dominant_orientation = 'same' if same_count >= len(orientations)/2 else 'opposite'
        
        total_length = sum(link['length'] for link in pair_links)
        
        pair_stats[pair_key] = {
            'chr1': chr1,
            'chr2': chr2,
            'link_count': len(pair_links),
            'min_start1': min_start1,
            'max_end1': max_end1,
            'min_start2': min_start2,
            'max_end2': max_end2,
            'span1': max_end1 - min_start1,
            'span2': max_end2 - min_start2,
            'total_length': total_length,
            'dominant_orientation': dominant_orientation,
            'orientation_ratio': same_count / len(orientations),
            'links': pair_links  # Keep the actual links for gap analysis
        }
    
    return pair_stats

def detect_coverage_gaps(pair_links, gap_threshold):
    """
    Detect gaps in link coverage and split into continuous segments.
    Returns list of segments, where each segment contains links that are continuous
    within the gap threshold.
    """
    if not pair_links:
        return []
    
    # Sort links by position on chromosome 1
    sorted_links = sorted(pair_links, key=lambda x: x['start1'])
    
    segments = []
    current_segment = [sorted_links[0]]
    
    for i in range(1, len(sorted_links)):
        prev_link = sorted_links[i-1]
        curr_link = sorted_links[i]
        
        # Check gap on chromosome 1
        gap_chr1 = curr_link['start1'] - prev_link['end1']
        
        # Check gap on chromosome 2 (consider orientation)
        if prev_link['orientation'] == curr_link['orientation'] == 'same':
            # Same orientation: gaps should be in same direction
            gap_chr2 = curr_link['start2'] - prev_link['end2']
        else:
            # Mixed or opposite orientation: more complex gap calculation
            # For simplicity, use absolute difference in positions
            gap_chr2 = abs(curr_link['start2'] - prev_link['start2']) - abs(curr_link['end2'] - prev_link['end2'])
        
        # If gap on either chromosome exceeds threshold, start new segment
        if gap_chr1 > gap_threshold or abs(gap_chr2) > gap_threshold:
            segments.append(current_segment)
            current_segment = [curr_link]
        else:
            current_segment.append(curr_link)
    
    # Add the final segment
    segments.append(current_segment)
    
    return segments

def generate_ribbon_links_with_colors(pair_stats, chr_colors, min_links=50, gap_threshold=50000):
    """Generate ribbon-style links using karyotype colors with gap-aware splitting."""
    ribbons = []
    
    sorted_pairs = sorted(pair_stats.items(), key=lambda x: x[1]['link_count'], reverse=True)
    
    print(f"\nTop 10 chromosome pairs by link count:")
    for i, (pair_key, stats) in enumerate(sorted_pairs[:10]):
        chr1_color = chr_colors.get(stats['chr1'], 'black')
        chr2_color = chr_colors.get(stats['chr2'], 'black')
        print(f"{i+1:2d}. {stats['chr1']} ({chr1_color}) ↔ {stats['chr2']} ({chr2_color}): {stats['link_count']} links")
    
    major_pairs = 0
    total_links_represented = 0
    total_segments_generated = 0
    large_gaps_detected = 0
    
    print(f"\nProcessing chromosome pairs with gap threshold: {gap_threshold:,} bp")
    
    for pair_key, stats in sorted_pairs:
        if stats['link_count'] >= min_links:
            # Detect coverage gaps and split into segments
            segments = detect_coverage_gaps(stats['links'], gap_threshold)
            
            # Get the color from the CCMP1545 chromosome (first chromosome)
            ribbon_color = chr_colors.get(stats['chr1'], 'black')
            
            # Generate ribbons for each continuous segment
            segment_ribbons_count = 0
            for i, segment in enumerate(segments):
                if len(segment) > 0:  # Only process non-empty segments
                    # Calculate boundaries for this segment
                    min_start1 = min(link['start1'] for link in segment)
                    max_end1 = max(link['end1'] for link in segment)
                    min_start2 = min(link['start2'] for link in segment)
                    max_end2 = max(link['end2'] for link in segment)
                    
                    # Calculate dominant orientation for this segment
                    orientations = [link['orientation'] for link in segment]
                    same_count = sum(1 for o in orientations if o == 'same')
                    dominant_orientation = 'same' if same_count >= len(orientations)/2 else 'opposite'
                    
                    ribbon = {
                        'chr1': stats['chr1'],
                        'start1': min_start1,
                        'end1': max_end1,
                        'chr2': stats['chr2'],
                        'start2': min_start2,
                        'end2': max_end2,
                        'color': ribbon_color,
                        'link_count': len(segment),
                        'orientation': dominant_orientation,
                        'segment_id': i + 1,
                        'total_segments': len(segments)
                    }
                    
                    ribbons.append(ribbon)
                    segment_ribbons_count += 1
                    total_segments_generated += 1
            
            # Track statistics
            major_pairs += 1
            total_links_represented += stats['link_count']
            
            # Report gaps detected
            gaps_for_this_pair = len(segments) - 1
            if gaps_for_this_pair > 0:
                large_gaps_detected += gaps_for_this_pair
                print(f"  {stats['chr1']} ↔ {stats['chr2']}: Split into {len(segments)} segments ({gaps_for_this_pair} gaps >{gap_threshold:,} bp)")
    
    print(f"\nGap-aware ribbon generation summary:")
    print(f"Processed {major_pairs} chromosome pairs with ≥{min_links} links")
    print(f"Large gaps detected: {large_gaps_detected}")
    print(f"Total ribbon segments generated: {total_segments_generated}")
    print(f"Total links represented: {total_links_represented:,}")
    print(f"Average segments per chromosome pair: {total_segments_generated/max(major_pairs,1):.1f}")
    
    return ribbons

def write_ribbon_links_file(ribbons, output_file):
    """Write ribbon links to circos-compatible format."""
    with open(output_file, 'w') as f:
        for ribbon in ribbons:
            # Format: chr1 start1 end1 chr2 start2 end2 color=colorname
            # Remove z= column and put color= at the end
            f.write(f"{ribbon['chr1']}\t{ribbon['start1']}\t{ribbon['end1']}\t")
            f.write(f"{ribbon['chr2']}\t{ribbon['start2']}\t{ribbon['end2']}\t")
            f.write(f"color={ribbon['color']}\n")
    
    print(f"Wrote {len(ribbons)} ribbon links to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate colored ribbon links for circos')
    parser.add_argument('--karyotype', default='circos-files/CCMP1545_vs_RCC1749_reordered_flipped.kar',
                       help='Karyotype file with color definitions')
    parser.add_argument('--input', default='circos-files/CCMP1545_vs_RCC1749.links.tsv',
                       help='Input links file')
    parser.add_argument('--output', default='circos-files/CCMP1545_vs_RCC1749_ribbons_colored.tsv',
                       help='Output colored ribbon links file')
    parser.add_argument('--min-links', type=int, default=50,
                       help='Minimum number of links required for ribbon generation')
    parser.add_argument('--gap-threshold', type=int, default=50000,
                       help='Maximum gap size (bp) between links before splitting ribbon (default: 50000)')
    parser.add_argument('--min-link-length', type=int, default=1000,
                       help='Minimum link length (bp) to include in ribbon generation (default: 1000)')
    parser.add_argument('--merge-threshold', type=int, default=5000,
                       help='Maximum gap size (bp) between links before merging them (default: 5000)')
    
    args = parser.parse_args()
    
    # Load colors and process links
    chr_colors = parse_karyotype_colors(args.karyotype)
    links = parse_links_file(args.input, args.min_link_length)
    
    # Merge nearby links before grouping
    if args.merge_threshold > 0:
        merged_links = merge_nearby_links(links, args.merge_threshold)
    else:
        merged_links = links
        print("Link merging disabled (merge-threshold = 0)")
    
    pair_stats = group_links_by_chromosome_pairs(merged_links)
    ribbons = generate_ribbon_links_with_colors(pair_stats, chr_colors, args.min_links, args.gap_threshold)
    write_ribbon_links_file(ribbons, args.output)
    
    print(f"\nColored ribbon generation complete!")
    print(f"Original links: {len(links):,}")
    if args.merge_threshold > 0:
        print(f"Merged links: {len(merged_links):,}")
    print(f"Final ribbons: {len(ribbons)}")
    print(f"Total reduction: {len(links)/len(ribbons):.1f}x")

if __name__ == "__main__":
    main()
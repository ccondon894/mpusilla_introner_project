import os
import pysam
import pandas as pd
from collections import defaultdict
import sys
import subprocess
import tempfile
import argparse

# Sample group definitions for Phase 2: Handling Divergence
GROUP1_SAMPLES = {"CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482",
                  "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"}
GROUP2_SAMPLES = {"RCC1749", "RCC3052"}

def load_introner_synteny_context(output_dir, samples):
    """
    Load introner synteny maps with 3-part genomic fingerprints.

    Reads files with format: introner_id, host_gene, upstream_gene, downstream_gene

    Returns: Dictionary mapping sample -> introner_id -> {host_gene, upstream_gene, downstream_gene}
    """
    synteny_context = defaultdict(lambda: defaultdict(dict))

    for sample in samples:
        context_file = os.path.join(output_dir, f"{sample}.introner_synteny_map.tsv")
        if not os.path.exists(context_file):
            print(f"Warning: Synteny map file not found for {sample}: {context_file}")
            continue

        try:
            with open(context_file, 'r') as f:
                header = f.readline().strip().split('\t')
                # Expected: introner_id, host_gene, upstream_gene, downstream_gene

                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 4:
                        continue

                    introner_id = fields[0]
                    host_gene = fields[1] if fields[1] != 'NA' else None
                    upstream_gene = fields[2] if fields[2] != 'NA' else None
                    downstream_gene = fields[3] if fields[3] != 'NA' else None

                    synteny_context[sample][introner_id] = {
                        'host_gene': host_gene,
                        'upstream_gene': upstream_gene,
                        'downstream_gene': downstream_gene
                    }

            print(f"Loaded {len(synteny_context[sample])} introner synteny maps for {sample}")

        except Exception as e:
            print(f"Error loading synteny map for {sample}: {e}")

    return synteny_context


def load_gene_annotations(processed_ann_dir, samples):
    """
    Load gene annotations from BED files for all samples.

    Args:
        processed_ann_dir: Directory containing {sample}.gene.bed files
        samples: List of sample names

    Returns:
        Dictionary mapping sample -> gene_name -> (contig, start, end, strand)
    """
    gene_annotations = defaultdict(dict)

    for sample in samples:
        bed_file = os.path.join(processed_ann_dir, f"{sample}.gene.bed")

        if not os.path.exists(bed_file):
            print(f"Warning: Gene BED file not found for {sample}: {bed_file}")
            continue

        try:
            with open(bed_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 6:
                        continue

                    contig = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    gene_name = fields[3]
                    strand = fields[5]

                    gene_annotations[sample][gene_name] = (contig, start, end, strand)

            print(f"Loaded {len(gene_annotations[sample])} gene annotations for {sample}")

        except Exception as e:
            print(f"Error loading gene annotations for {sample}: {e}")

    return gene_annotations


def load_introner_metadata(bed_dir, samples):
    """
    Load introner metadata from BED files.
    
    Returns: Dictionary mapping sample -> contig -> (start, end) -> metadata
    """
    metadata = defaultdict(lambda: defaultdict(dict))
    
    for sample in samples:
        bed_file = os.path.join(bed_dir, f"{sample}.candidate_loci.filtered.bed")
        if not os.path.exists(bed_file):
            print(f"Warning: BED file not found for {sample}")
            continue
            
        with open(bed_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    continue
                
                # Skip header lines
                if fields[0].startswith('#') or fields[0] == 'chrom':
                    continue
                    
                contig = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                introner_id = fields[3]

                # Parse metadata fields based on BED format
                metadata_dict = {
                    'introner_id': introner_id,
                    'contig': contig,
                    'start': start,
                    'end': end
                }
                
                # Handle different BED formats (with or without structural analysis)
                if len(fields) >= 9:  # Basic format: chrom, start, end, name, family, similarity, orientation, gene_info, splice_info
                    # Parse family (column 4)
                    family_val = fields[4] if fields[4] != 'NA' else ''
                    metadata_dict['family'] = family_val
                    
                    # Parse similarity (column 5)
                    similarity_val = fields[5] if fields[5] != 'NA' else ''
                    metadata_dict['similarity'] = similarity_val
                        
                    # Parse orientation (column 6)
                    orientation_val = fields[6] if fields[6] != 'NA' else ''
                    metadata_dict['orientation'] = orientation_val
                    
                    # Parse gene_info (column 7) and splice_info (column 8)
                    gene_info = fields[7] if fields[7] != 'NA' else ''
                    splice_info = fields[8] if fields[8] != 'NA' else ''
                    
                    metadata_dict['gene'] = gene_info
                    metadata_dict['splice_site'] = splice_info
                    
                    # Debug: Print first few parsed entries for the first sample processed
                    if sample not in metadata and len(metadata_dict) > 3:  # First entry for this sample
                        print(f"DEBUG: Parsed metadata for {sample} {contig}:{start}-{end}")
                        print(f"DEBUG: family='{family_val}', gene='{gene_info}', splice_site='{splice_info}'")
                    
                    # For enhanced format, also parse structural analysis columns
                    if len(fields) >= 14:  # Enhanced format with TIR/TSD columns
                        # Enhanced format: chrom, start, end, name, family, similarity, orientation, 
                        # tir_info, tsd_info, classification, struct_confidence, validation_info, gene_info, splice_info
                        
                        tir_info = fields[7] if fields[7] != 'NA' else ''
                        metadata_dict['tir_info'] = tir_info
                            
                        tsd_info = fields[8] if fields[8] != 'NA' else ''
                        metadata_dict['tsd_info'] = tsd_info
                            
                        classification = fields[9] if fields[9] != 'NA' else ''
                        metadata_dict['classification'] = classification
                        
                        # For enhanced format, gene and splice_site are in different positions
                        gene_info = fields[12] if len(fields) > 12 and fields[12] != 'NA' else ''
                        splice_info = fields[13] if len(fields) > 13 and fields[13] != 'NA' else ''
                        
                        # Override the basic format values with enhanced format positions
                        metadata_dict['gene'] = gene_info
                        metadata_dict['splice_site'] = splice_info
                
                # Fallback: try parsing key=value format for backward compatibility
                else:
                    for field in fields[4:]:
                        if '=' in field:
                            key, value = field.split('=', 1)
                            metadata_dict[key] = value
                
                # Store metadata using coordinates as key
                metadata[sample][contig][(start, end)] = metadata_dict
    
    return metadata

def find_overlapping_introner(metadata, sample, contig, target_start, target_end, min_overlap_fraction=0.6):
    """
    Find an introner in the metadata that overlaps with the target coordinates.

    Args:
        metadata: Introner metadata dictionary
        sample: Sample name
        contig: Contig name
        target_start: Start position of target region
        target_end: End position of target region
        min_overlap_fraction: Minimum fraction of overlap required (relative to the shorter region)

    Returns:
        Tuple of (start, end) coordinates of matching introner if found, None otherwise
    """
    if sample not in metadata or contig not in metadata[sample]:
        return None

    best_match = None
    best_overlap_fraction = 0

    # Check each introner in this contig
    for (start, end), _ in metadata[sample][contig].items():
        # Calculate overlap
        overlap_start = max(start, target_start)
        overlap_end = min(end, target_end)

        if overlap_start <= overlap_end:  # Regions overlap
            overlap_length = overlap_end - overlap_start + 1

            # Calculate overlap as a fraction of the shorter region
            target_length = target_end - target_start + 1
            introner_length = end - start + 1
            shorter_length = min(target_length, introner_length)

            overlap_fraction = overlap_length / shorter_length

            # Update best match if this one has a better overlap
            if overlap_fraction > best_overlap_fraction:
                best_overlap_fraction = overlap_fraction
                best_match = (start, end)

    # Return the best match if it meets the minimum overlap threshold
    if best_overlap_fraction >= min_overlap_fraction:
        return best_match

    return None

def find_all_overlapping_introners(metadata, sample, contig, target_start, target_end, min_overlap_fraction=0.6):
    """
    Find all introners in the metadata that overlap with the target coordinates.

    Args:
        metadata: Introner metadata dictionary
        sample: Sample name
        contig: Contig name
        target_start: Start position of target region
        target_end: End position of target region
        min_overlap_fraction: Minimum fraction of overlap required (relative to the shorter region)

    Returns:
        List of (start, end) coordinates of matching introners, sorted by overlap fraction (best first)
    """
    if sample not in metadata or contig not in metadata[sample]:
        return []

    matches = []

    # Check each introner in this contig
    for (start, end), _ in metadata[sample][contig].items():
        # Calculate overlap
        overlap_start = max(start, target_start)
        overlap_end = min(end, target_end)

        if overlap_start <= overlap_end:  # Regions overlap
            overlap_length = overlap_end - overlap_start + 1

            # Calculate overlap as a fraction of the shorter region
            target_length = target_end - target_start + 1
            introner_length = end - start + 1
            shorter_length = min(target_length, introner_length)

            overlap_fraction = overlap_length / shorter_length

            # Keep match if it meets the minimum overlap threshold
            if overlap_fraction >= min_overlap_fraction:
                matches.append(((start, end), overlap_fraction))

    # Sort by overlap fraction (descending) and return just the coordinates
    matches.sort(key=lambda x: x[1], reverse=True)
    return [coords for coords, _ in matches]

def parse_query_name(query_name, flank_type):
    """Extract information from the query name."""
    # First, split by the pipe character
    parts = query_name.split('|')
    if len(parts) < 2:
        return None, None, {}
    
    # Get coordinates part (before the pipe)
    coords = parts[0]
    
    # Parse the coordinates to extract contig, start, end
    try:
        contig, pos_range = coords.split(':')
        start, end = map(int, pos_range.split('-'))
    except (ValueError, IndexError):
        contig, start, end = None, None, None
    
    # Get the introner ID and metadata (after the pipe)
    introner_info = parts[1].replace(f"_{flank_type}", "")
    
    # Split the introner info into ID and metadata
    introner_parts = introner_info.split()
    introner_id = introner_parts[0]
    
    # Parse metadata if available
    metadata = {}
    for part in introner_parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key] = value
    
    # Create a unique ID by combining introner_id with coordinates
    unique_id = f"{introner_id}_{contig}_{start}_{end}" if contig and start and end else introner_id
    
    return coords, unique_id, metadata

def extract_best_alignments(bam_path, flank_type, export_quality_info=False):
    """
    Extract best alignments for each query from a BAM file.
    
    Args:
        bam_path: Path to the BAM file
        flank_type: Type of flank ('left' or 'right')
        export_quality_info: Whether to export detailed alignment quality info
        
    Returns: Dictionary mapping query_id -> alignment information
    """
    alignments = {}
    all_query_ids = set()  # Track all query IDs seen
    filtered_query_ids = set()  # Track query IDs retained after filtering
    low_quality_alignments = 0  # Track how many alignments fail quality checks
    alignment_quality_info = []  # Detailed quality info for problematic alignments
    
    if not os.path.exists(bam_path):
        return alignments
        
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
                
            coords, query_id, query_metadata = parse_query_name(read.query_name, flank_type)
            if not query_id:
                continue
            
            all_query_ids.add(query_id)
            
            # Check if alignment meets quality criteria
            if not is_good_alignment(read):
                low_quality_alignments += 1
                # If requested, collect detailed quality information for debugging
                if export_quality_info:
                    quality_info = analyze_alignment_quality(read)
                    quality_info['query_id'] = query_id
                    quality_info['flank_type'] = flank_type
                    alignment_quality_info.append(quality_info)
                continue
            
            # Store or update alignment if it's better than previous one
            current_quality = read.mapping_quality
            if query_id not in alignments or current_quality > alignments[query_id]['quality']:
                alignments[query_id] = {
                    'query_coords': coords,
                    'target_contig': bam.get_reference_name(read.reference_id),
                    'target_start': read.reference_start,
                    'target_end': read.reference_end,
                    'quality': current_quality,
                    'query_metadata': query_metadata,
                    'is_reverse': read.is_reverse
                }
                filtered_query_ids.add(query_id)
    
    # Print debug information
    print(f"BAM file: {bam_path}")
    print(f"Total unique query IDs found: {len(all_query_ids)}")
    print(f"Query IDs retained after quality filtering: {len(filtered_query_ids)}")
    print(f"Query IDs lost during filtering: {len(all_query_ids - filtered_query_ids)}")
    print(f"Low quality alignments rejected: {low_quality_alignments} (MAPQ < 30 or matched bases < 70)")
    
    if len(all_query_ids) > 0 and len(filtered_query_ids) > 0:
        # Print a sample of problematic query IDs (up to 5)
        lost_ids = list(all_query_ids - filtered_query_ids)
        if lost_ids:
            print(f"Sample of lost query IDs: {lost_ids[:min(5, len(lost_ids))]}")
    
    # Return different results depending on whether quality info was requested
    if export_quality_info:
        return alignments, alignment_quality_info
    else:
        return alignments


# ==================== Phase 2: Synteny-Guided Rescue Functions ====================

def get_synteny_search_window(query_synteny, gene_annotations, target_sample, window_size=5000):
    """
    Calculate genomic search window using synteny context.

    Args:
        query_synteny: Synteny context dict with host_gene, upstream_gene, downstream_gene
        gene_annotations: Gene annotations dict [sample][gene_name] -> (contig, start, end, strand)
        target_sample: Target sample name
        window_size: Window size to expand around host gene (default 5000bp)

    Returns:
        (target_contig, search_start, search_end) or None if not found
    """
    if not query_synteny or target_sample not in gene_annotations:
        return None

    # Try to find host gene first
    host_gene = query_synteny.get('host_gene')
    if host_gene and host_gene in gene_annotations[target_sample]:
        contig, start, end, strand = gene_annotations[target_sample][host_gene]
        # Expand window around host gene
        search_start = max(0, start - window_size)
        search_end = end + window_size
        return (contig, search_start, search_end)

    # If host gene not found, try upstream gene
    upstream_gene = query_synteny.get('upstream_gene')
    if upstream_gene and upstream_gene in gene_annotations[target_sample]:
        contig, start, end, strand = gene_annotations[target_sample][upstream_gene]
        search_start = max(0, start - window_size)
        search_end = end + window_size * 2  # Larger window since we're searching downstream
        return (contig, search_start, search_end)

    # If upstream not found, try downstream gene
    downstream_gene = query_synteny.get('downstream_gene')
    if downstream_gene and downstream_gene in gene_annotations[target_sample]:
        contig, start, end, strand = gene_annotations[target_sample][downstream_gene]
        search_start = max(0, start - window_size * 2)  # Larger window since we're searching upstream
        search_end = end + window_size
        return (contig, search_start, search_end)

    # No syntenic genes found in target
    return None


def is_good_alignment_relaxed(read, min_mapq=20, min_match_length=50):
    """
    Check if alignment meets relaxed quality criteria for distant samples.

    Uses lower thresholds than standard is_good_alignment():
    - MAPQ >= 20 (vs 30)
    - Matched bases >= 50 (vs 70)

    Args:
        read: pysam AlignedSegment
        min_mapq: Minimum mapping quality (default 20)
        min_match_length: Minimum number of matched bases (default 50)

    Returns:
        True if alignment passes quality filters
    """
    # Check mapping quality
    if read.mapping_quality < min_mapq:
        return False

    # Count matched bases from CIGAR string
    matched_bases = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 0:  # M = alignment match/mismatch
                matched_bases += length

    return matched_bases >= min_match_length


def synteny_guided_rescue_from_bam(query, target, query_id, query_metadata,
                                  synteny_context, gene_annotations, rescue_bam_file):
    """
    Fast rescue using pre-computed BAM file with indexed lookups.

    Uses pysam.fetch() to instantly retrieve alignments in syntenic window.

    Args:
        query: Query sample name
        target: Target sample name
        query_id: Query introner ID
        query_metadata: Metadata for query introner
        synteny_context: Synteny context dictionary
        gene_annotations: Gene annotations dictionary
        rescue_bam_file: Path to pre-computed rescue BAM file

    Returns:
        Updated match dict with new scenario, or None if rescue fails
    """
    # Debug first 3 rescues
    global RESCUE_DEBUG_COUNT
    if 'RESCUE_DEBUG_COUNT' not in globals():
        RESCUE_DEBUG_COUNT = 0

    # Step 1: Get query introner's synteny context
    query_contig = query_metadata.get('contig', '')
    query_start = query_metadata.get('start', 0)
    query_end = query_metadata.get('end', 0)
    query_introner_id = query_metadata.get('introner_id', '')

    if RESCUE_DEBUG_COUNT < 3:
        print(f"      [DEBUG {RESCUE_DEBUG_COUNT}] Rescue for {query_id}")
        print(f"        Metadata: introner_id={query_introner_id}, contig={query_contig}, start={query_start}, end={query_end}")
        RESCUE_DEBUG_COUNT += 1

    query_unique_id = f"{query_introner_id}|{query_contig}|{query_start}|{query_end}"
    query_synteny = synteny_context.get(query, {}).get(query_unique_id)

    if not query_synteny:
        return None

    # Step 2: Find syntenic search window in target genome
    search_window = get_synteny_search_window(query_synteny, gene_annotations, target)

    if not search_window:
        return None

    target_contig, search_start, search_end = search_window

    # Step 3: Extract introner_seq identifier for matching
    # Query IDs have format: introner_seq_XXXX_SAMPLE#X#contig_start_end
    # BAM read names have format: SAMPLE#X#contig:coords|introner_seq_XXXX_left/right
    # We match on the introner_seq_XXXX part
    introner_seq_id = None
    if 'introner_seq_' in query_id:
        # Extract "introner_seq_XXXX" part
        parts = query_id.split('_')
        if len(parts) >= 3:
            introner_seq_id = f"{parts[0]}_{parts[1]}_{parts[2]}"  # "introner_seq_XXXX"

    if not introner_seq_id:
        return None

    # Step 4: Open BAM and fetch alignments in syntenic window
    # This is INSTANT thanks to the .bai index!
    try:
        bamfile = pysam.AlignmentFile(rescue_bam_file, "rb")

        left_aln = None
        right_aln = None

        # pysam.fetch() uses the index to jump directly to this region
        for read in bamfile.fetch(target_contig, search_start, search_end):
            # Check if this is one of our flanks by matching introner_seq_XXXX
            if introner_seq_id in read.query_name:
                if '_left' in read.query_name:
                    # Apply relaxed quality filter
                    if is_good_alignment_relaxed(read):
                        left_aln = read
                elif '_right' in read.query_name:
                    if is_good_alignment_relaxed(read):
                        right_aln = read

            # Optimization: stop if both found
            if left_aln and right_aln:
                break

        bamfile.close()

        # Step 5: Process results
        if left_aln and right_aln:
            # Calculate gap size
            gap_start = left_aln.reference_end
            gap_end = right_aln.reference_start
            gap_size = gap_end - gap_start

            scenario = 1 if gap_size >= 50 else 2

            result = {
                'scenario': scenario,
                'query_contig': query_contig,
                'query_start': query_start,
                'query_end': query_end,
                'target_contig': target_contig,
                'target_start': gap_start,
                'target_end': gap_end,
                'query_metadata': query_metadata,
                'target_metadata': {},
                'match_score': 0.5,
                'rescued': True
            }
            return result
        else:
            return None

    except Exception as e:
        # Log errors for debugging (temporarily)
        import traceback
        print(f"      Rescue failed for {query_id}: {str(e)}")
        # traceback.print_exc()  # Uncomment for detailed debugging
        return None


def build_ortholog_dictionary(samples, bam_dir, bed_dir, genome_dir=None,
                             gene_annotations=None, export_quality_info=False, mode='rescue'):
    """
    Build a comprehensive dictionary of orthologous relationships.

    Args:
        samples: List of sample names
        bam_dir: Directory containing BAM files
        bed_dir: Directory containing BED files
        genome_dir: Directory containing genome FASTA files (for Phase 2 rescue)
        gene_annotations: Gene annotations dictionary (for Phase 2 rescue)
        export_quality_info: Whether to export detailed alignment quality info
        mode: 'initial' for Phase 1 only, 'rescue' for Phase 1 + Phase 2 rescue

    Returns: Dictionary mapping query -> target -> query_id -> alignment info
             And optionally a DataFrame with alignment quality information
    """
    # Load introner metadata from BED files
    print("Loading introner metadata...")
    introner_metadata = load_introner_metadata(bed_dir, samples)

    # Load introner synteny context for disambiguation
    print("Loading introner synteny context...")
    synteny_context = load_introner_synteny_context(bed_dir, samples)
    
    # Initialize result dictionary
    ortholog_dict = defaultdict(lambda: defaultdict(dict))
    
    # Track query IDs for debugging
    total_query_ids_processed = 0
    total_query_ids_kept = 0
    query_ids_by_scenario = {1: 0, 2: 0, 3: 0}  # Track counts by scenario
    
    # Collect alignment quality information if requested
    all_quality_info = []
    
    # Process each query-target pair
    print("Processing alignments...")
    for query in samples:
        for target in samples:
            if query == target:
                continue  # Skip self-comparisons
                
            print(f"Processing {query} -> {target}")
            
            # Get left and right flank alignments
            left_bam = os.path.join(bam_dir, f"{query}_{target}.left_flank.sorted.bam")
            right_bam = os.path.join(bam_dir, f"{query}_{target}.right_flank.sorted.bam")
            if not os.path.exists(left_bam) or not os.path.exists(right_bam):
                print(f"Warning: BAM files not found for {query} -> {target}")
                print(left_bam)
                print(right_bam)
                continue  # Skip this pair if files don't exist
                
            # Extract alignments, potentially with quality info
            if export_quality_info:
                left_alignments, left_quality_info = extract_best_alignments(left_bam, "left", export_quality_info=True)
                right_alignments, right_quality_info = extract_best_alignments(right_bam, "right", export_quality_info=True)
                
                # Collect quality info with query/target identifiers
                for info in left_quality_info + right_quality_info:
                    info['query_sample'] = query
                    info['target_sample'] = target
                    all_quality_info.append(info)
            else:
                left_alignments = extract_best_alignments(left_bam, "left")
                right_alignments = extract_best_alignments(right_bam, "right")
            
            # Find query IDs in either left or right alignments
            all_query_ids = set(left_alignments.keys()) | set(right_alignments.keys())
            print(f"Total unique query IDs for {query} -> {target}: {len(all_query_ids)}")
            
            total_query_ids_processed += len(all_query_ids)
            
            # Process each query ID
            for query_id in all_query_ids:
                left_aln = left_alignments.get(query_id)
                right_aln = right_alignments.get(query_id)
                
                # Get query metadata
                query_contig, coords_range = left_aln['query_coords'].split(':') if left_aln else (right_aln['query_coords'].split(':') if right_aln else (None, None))
                query_start, query_end = map(int, coords_range.split('-')) if coords_range else (None, None)
                query_coords = None
                
                # If missing either left or right alignment, mark as scenario 3
                if not left_aln or not right_aln:
                    ortholog_dict[query][target][query_id] = {
                        'scenario': 3,
                        'query_id': query_id,
                        'query_coords': left_aln['query_coords'] if left_aln else (right_aln['query_coords'] if right_aln else None),
                        'query_contig': query_contig,
                        'query_start': query_start,
                        'query_end': query_end,
                        'target_contig': None,
                        'target_start': None,
                        'target_end': None,
                        'query_metadata': None,
                        'target_metadata': None,
                        'orientation': 'unknown',
                        'left_reverse': None,
                        'right_reverse': None
                    }
                    query_ids_by_scenario[3] += 1
                    continue
                
                # Check if alignments are on same contig
                if left_aln['target_contig'] != right_aln['target_contig']:
                    ortholog_dict[query][target][query_id] = {
                        'scenario': 3,
                        'query_id': query_id,
                        'query_coords': left_aln['query_coords'],
                        'query_contig': query_contig,
                        'query_start': query_start,
                        'query_end': query_end,
                        'target_contig': None,
                        'target_start': None,
                        'target_end': None,
                        'query_metadata': None,
                        'target_metadata': None,
                        'orientation': 'unknown',
                        'left_reverse': None,
                        'right_reverse': None
                    }
                    query_ids_by_scenario[3] += 1
                    continue
                
                # Both flanks aligned to same contig
                target_contig = left_aln['target_contig']
                
                # Get orientation information early for filtering
                left_reverse = left_aln.get('is_reverse', False)
                right_reverse = right_aln.get('is_reverse', False)
                
                # Calculate the gap between flanks (the potential introner region)
                # We need to find which flank comes first in the target sequence
                if left_aln['target_start'] <= right_aln['target_start']:
                    # Left flank comes first
                    gap_start = left_aln['target_end']
                    gap_end = right_aln['target_start']
                else:
                    # Right flank comes first (or they overlap)
                    gap_start = right_aln['target_end']
                    gap_end = left_aln['target_start']
                
                # Calculate the actual gap size between flanks
                gap_size = max(0, gap_end - gap_start)
                
                # Calculate the total span from leftmost start to rightmost end
                target_start = min(left_aln['target_start'], right_aln['target_start'])
                target_end = max(left_aln['target_end'], right_aln['target_end'])
                total_span = target_end - target_start
                
                # Filter out alignments where flanks are too far apart (likely misalignments)
                # Only apply to forward alignments where we expect left < right
                max_reasonable_span = 500
                if (not left_reverse and not right_reverse and 
                    total_span > max_reasonable_span):
                    ortholog_dict[query][target][query_id] = {
                        'scenario': 3,  # Treat as missing data due to misalignment
                        'query_id': query_id,
                        'query_coords': left_aln['query_coords'],
                        'query_contig': query_contig,
                        'query_start': query_start,
                        'query_end': query_end,
                        'target_contig': None,
                        'target_start': None,
                        'target_end': None,
                        'query_metadata': None,
                        'target_metadata': None,
                        'orientation': 'unknown',
                        'left_reverse': None,
                        'right_reverse': None
                    }
                    query_ids_by_scenario[3] += 1
                    continue
                
                # Two criteria for absence:
                # 1. Gap between flanks is too small (< 50bp)
                # 2. Total span is too small (< 200bp) 
                min_gap_size = 50
                min_total_span = 200
                
                if gap_size < min_gap_size or total_span < min_total_span:
                    # Either gap between flanks is too small OR total span is too small - classify as absent (scenario 2)
                    target_coords = None
                    target_metadata = None
                    scenario = 2
                else:
                    # Both gap and total span are large enough - check if target has an introner at this location
                    target_coords = find_overlapping_introner(
                        introner_metadata, target, target_contig, target_start, target_end
                    )
                    target_metadata = introner_metadata[target][target_contig].get(target_coords) if target_coords else None
                    scenario = 1 if target_metadata else 2

                # Find query coordinates in metadata
                query_coords = None
                for (start, end), meta in introner_metadata[query][query_contig].items():
                    if start <= query_start <= end and start <= query_end <= end:
                        query_coords = (start, end)
                        break

                # Get query metadata - prefer from BAM query name, fall back to BED file
                query_name_metadata = left_aln.get('query_metadata', {}) or right_aln.get('query_metadata', {})
                query_bed_metadata = introner_metadata[query][query_contig].get(query_coords) if query_coords else None

                # Combine the metadata, prioritizing the query name metadata
                query_metadata = {}
                if query_bed_metadata:
                    query_metadata.update(query_bed_metadata)
                if query_name_metadata:
                    query_metadata.update(query_name_metadata)

                # Apply synteny-based disambiguation for multiple target candidates
                query_gene = query_metadata.get('gene')
                if scenario == 1 and query_gene:
                    # Check if there are multiple overlapping introners at this location
                    all_target_coords = find_all_overlapping_introners(
                        introner_metadata, target, target_contig, target_start, target_end
                    )

                    if len(all_target_coords) > 1:
                        # Multiple candidates found - use synteny context for disambiguation
                        best_target_coords = None
                        best_synteny_score = -1

                        # Get query synteny context
                        query_unique_id = f"{query_metadata.get('introner_id')}|{query_contig}|{query_start}|{query_end}"
                        query_synteny = synteny_context[query].get(query_unique_id)

                        for candidate_coords in all_target_coords:
                            candidate_metadata = introner_metadata[target][target_contig].get(candidate_coords)
                            if not candidate_metadata:
                                continue

                            # Get candidate synteny context (3-part fingerprint)
                            candidate_introner_id = candidate_metadata.get('introner_id')
                            unique_candidate_id = f"{candidate_introner_id}|{target_contig}|{candidate_coords[0]}|{candidate_coords[1]}"

                            if target in synteny_context and unique_candidate_id in synteny_context[target]:
                                candidate_synteny = synteny_context[target][unique_candidate_id]

                                # Score based on matching synteny components
                                synteny_score = 0

                                # Check host gene match
                                if query_synteny and candidate_synteny.get('host_gene') == query_gene:
                                    synteny_score = 1  # Host gene matches

                                # Update best candidate if this one has better synteny match
                                if synteny_score > best_synteny_score:
                                    best_synteny_score = synteny_score
                                    best_target_coords = candidate_coords
                            else:
                                # No synteny data for this candidate
                                if best_synteny_score < 0:
                                    best_target_coords = candidate_coords
                                    best_synteny_score = 0

                        # Use the best synteny-matched target
                        if best_target_coords and best_target_coords != target_coords:
                            target_coords = best_target_coords
                            target_metadata = introner_metadata[target][target_contig].get(target_coords)
                
                # Determine orientation pattern
                orientation = determine_orientation_pattern(left_reverse, right_reverse)
                
                # Use target coordinates from metadata if found, otherwise use alignment coordinates
                if target_coords and scenario == 1:
                    target_start_final = target_coords[0]
                    target_end_final = target_coords[1]
                else:
                    target_start_final = target_start
                    target_end_final = target_end
                
                ortholog_dict[query][target][query_id] = {
                    'scenario': scenario,
                    'query_id': query_id,
                    'query_coords': left_aln['query_coords'],
                    'query_contig': query_contig,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_contig': target_contig,
                    'target_start': target_start_final,
                    'target_end': target_end_final,
                    'query_metadata': query_metadata,
                    'target_metadata': target_metadata,
                    'orientation': orientation,
                    'left_reverse': left_reverse,
                    'right_reverse': right_reverse
                }
                query_ids_by_scenario[scenario] += 1
                total_query_ids_kept += 1
            
            # Add unaligned introners from metadata
            ortholog_dict = add_unaligned_introners_to_ortholog_dict(ortholog_dict, introner_metadata, query, target)

            # Phase 2: Synteny-guided rescue for distant pairs using pre-computed BAM files
            # Skip Phase 2 if running in 'initial' mode
            is_distant_pair = ((query in GROUP1_SAMPLES and target in GROUP2_SAMPLES) or
                              (query in GROUP2_SAMPLES and target in GROUP1_SAMPLES))

            if mode == 'rescue' and is_distant_pair and gene_annotations:
                # Check if rescue BAM exists
                rescue_bam = os.path.join(bam_dir, "rescue_alignments",
                                         f"{query}_to_{target}.rescue.sorted.bam")

                if os.path.exists(rescue_bam):
                    # Collect all scenario 3 cases for this query-target pair
                    scenario_3_ids = [qid for qid, data in ortholog_dict[query][target].items()
                                    if data.get('scenario') == 3]

                    if scenario_3_ids:
                        print(f"  Attempting synteny-guided rescue for {query}→{target}: "
                              f"{len(scenario_3_ids)} missing loci")

                        rescued_count = 0
                        rescued_to_present = 0
                        rescued_to_absent = 0

                        for query_id in scenario_3_ids:
                            # Parse query_id to get contig and coordinates
                            # Format: introner_seq_XXXX_SAMPLE#X#contig_start_end
                            parts = query_id.split('_')
                            if len(parts) < 4:
                                continue

                            # Extract contig and coordinates
                            remainder = '_'.join(parts[3:])  # "CCMP1545#0#scaffold_18_487378_487748"
                            coords_parts = remainder.rsplit('_', 2)  # Split from right to get last two _
                            if len(coords_parts) != 3:
                                continue

                            contig = coords_parts[0]
                            start = int(coords_parts[1])
                            end = int(coords_parts[2])

                            # Look up metadata correctly
                            query_meta = introner_metadata[query].get(contig, {}).get((start, end), {})

                            # Attempt rescue using fast BAM-based approach
                            rescued_result = synteny_guided_rescue_from_bam(
                                query, target, query_id, query_meta,
                                synteny_context, gene_annotations, rescue_bam
                            )

                            if rescued_result:
                                # Update the ortholog dict with rescued result
                                ortholog_dict[query][target][query_id].update(rescued_result)
                                rescued_count += 1

                                # Track which scenario it was rescued to
                                if rescued_result['scenario'] == 1:
                                    rescued_to_present += 1
                                    query_ids_by_scenario[1] += 1
                                    query_ids_by_scenario[3] -= 1
                                elif rescued_result['scenario'] == 2:
                                    rescued_to_absent += 1
                                    query_ids_by_scenario[2] += 1
                                    query_ids_by_scenario[3] -= 1

                        if rescued_count > 0:
                            print(f"    Successfully rescued {rescued_count}/{len(scenario_3_ids)} loci "
                                  f"(Present: {rescued_to_present}, Absent: {rescued_to_absent})")
                elif len([qid for qid, data in ortholog_dict[query][target].items()
                         if data.get('scenario') == 3]) > 0:
                    print(f"  Warning: Rescue BAM not found for {query}→{target}: {rescue_bam}")

    # Print debug information
    print(f"Total query IDs processed: {total_query_ids_processed}")
    print(f"Total query IDs kept: {total_query_ids_kept}")
    for scenario, count in query_ids_by_scenario.items():
        print(f"Scenario {scenario}: {count} query IDs")
    
    # Return different results depending on whether quality info was requested
    if export_quality_info and all_quality_info:
        quality_df = pd.DataFrame(all_quality_info)
        return ortholog_dict, quality_df
    else:
        return ortholog_dict

def select_best_matches(ortholog_dict):
    """
    For cases with multiple potential targets, select the best matches based on metadata.
    """
    for query in ortholog_dict:
        for target in ortholog_dict[query]:
            for query_id, match in ortholog_dict[query][target].items():
                # Skip cases without both metadata
                if match['scenario'] != 1 or not match['query_metadata'] or not match['target_metadata']:
                    continue
                    
                # Compare gene and introner_id
                score = 0
                
                # Check gene match
                q_gene = match['query_metadata'].get('gene')
                t_gene = match['target_metadata'].get('gene')
                if q_gene and t_gene and q_gene == t_gene:
                    score += 10
                    
                # Check introner_id match
                q_id = match['query_metadata'].get('introner_id')
                t_id = match['target_metadata'].get('introner_id')
                if q_id and t_id and q_id == t_id:
                    score += 5
                    
                # Store score
                match['match_score'] = score
                
    return ortholog_dict

def get_all_introners_from_metadata(metadata, query):
    """
    Extract all introners for a given query sample from metadata.
    
    Args:
        metadata: Dictionary of introner metadata
        query: Query sample name
        
    Returns:
        List of dictionaries with introner information
    """
    introners = []
    
    if query not in metadata:
        return introners
        
    # Iterate through all contigs and introners for this query
    for contig, introner_dict in metadata[query].items():
        for (start, end), meta_dict in introner_dict.items():
            # Create a unique ID for this introner
            introner_id = meta_dict.get('introner_id', 'unknown')
            unique_id = f"{introner_id}_{contig}_{start}_{end}"
            
            # Create introner entry
            introner_info = {
                'query_id': unique_id,
                'query_coords': f"{contig}:{start}-{end}",
                'query_contig': contig,
                'query_start': start,
                'query_end': end,
                'query_metadata': meta_dict
            }
            
            introners.append(introner_info)
            
    return introners

def add_unaligned_introners_to_ortholog_dict(ortholog_dict, introner_metadata, query, target):
    """
    Add introners from metadata that don't have alignments to the ortholog dictionary.
    
    Args:
        ortholog_dict: Dictionary of aligned introners
        introner_metadata: Dictionary of introner metadata
        query: Query sample name
        target: Target sample name
        
    Returns:
        Updated ortholog dictionary
    """
    # Get all introners for this query from metadata
    all_introners = get_all_introners_from_metadata(introner_metadata, query)
    
    # Get set of query_ids already in the ortholog dictionary for this pair
    existing_query_ids = set(ortholog_dict[query][target].keys())
    
    # Count of added introners
    added_count = 0
    
    # Add introners not yet in ortholog_dict
    for introner in all_introners:
        query_id = introner['query_id']
        
        if query_id not in existing_query_ids:
            # Add as scenario 3 (no alignment)
            ortholog_dict[query][target][query_id] = {
                'scenario': 3,
                'query_id': query_id,
                'query_coords': introner['query_coords'],
                'query_contig': introner['query_contig'],
                'query_start': introner['query_start'],
                'query_end': introner['query_end'],
                'target_contig': None,
                'target_start': None,
                'target_end': None,
                'query_metadata': introner['query_metadata'],
                'target_metadata': None,
                'orientation': 'unknown',
                'left_reverse': None,
                'right_reverse': None
            }
            added_count += 1
            
    if added_count > 0:
        print(f"  Added {added_count} unaligned introners from metadata for {query} -> {target}")
        
    return ortholog_dict

def is_good_alignment(read, min_mapq=30, min_match_length=70):
    """
    Check if an alignment is of sufficient quality.
    
    Args:
        read: pysam AlignedSegment object
        min_mapq: Minimum mapping quality score (default: 30)
        min_match_length: Minimum number of matched bases in CIGAR (default: 70)
    
    Returns:
        Boolean indicating if alignment meets quality criteria
    """
    # Check mapping quality
    if read.mapping_quality < min_mapq:
        return False
    
    # Calculate total matched bases from CIGAR
    matched_bases = 0
    for op, length in read.cigartuples:
        # CIGAR operations: 0=M (match/mismatch), 7== (match), 8=X (mismatch)
        if op in [0, 7, 8]:
            matched_bases += length
    
    # Check if matched bases meet minimum threshold
    return matched_bases >= min_match_length

def analyze_alignment_quality(read):
    """
    Analyze an alignment's quality and return diagnostic information.
    
    Args:
        read: pysam AlignedSegment object
    
    Returns:
        Dictionary with alignment quality metrics
    """
    # CIGAR operation types
    cigar_ops = {
        0: 'M',  # match or mismatch
        1: 'I',  # insertion
        2: 'D',  # deletion
        3: 'N',  # skipped region
        4: 'S',  # soft clip
        5: 'H',  # hard clip
        6: 'P',  # padding
        7: '=',  # sequence match
        8: 'X'   # sequence mismatch
    }
    
    # Count each operation type
    op_counts = defaultdict(int)
    total_length = 0
    
    for op, length in read.cigartuples:
        op_char = cigar_ops.get(op, '?')
        op_counts[op_char] += length
        total_length += length if op != 5 else 0  # exclude hard clips from total length
    
    # Calculate aligned length (M/=/X operations)
    aligned_length = op_counts.get('M', 0) + op_counts.get('=', 0) + op_counts.get('X', 0)
    
    # Calculate alignment percentage (aligned / total sequence length)
    align_percent = (aligned_length / total_length) * 100 if total_length > 0 else 0
    
    return {
        'mapq': read.mapping_quality,
        'cigar_str': read.cigarstring,
        'aligned_length': aligned_length,
        'total_length': total_length,
        'align_percent': f"{align_percent:.1f}%",
        'op_counts': dict(op_counts),
        'read_name': read.query_name,
        'reverse_strand': read.is_reverse,
        'ref_name': read.reference_name,
        'ref_start': read.reference_start,
        'ref_end': read.reference_end
    }

def determine_orientation_pattern(left_reverse, right_reverse):
    """
    Determine the overall orientation pattern for an introner based on 
    left and right flank alignment orientations.
    
    Args:
        left_reverse: Boolean indicating if left flank aligned to reverse strand
        right_reverse: Boolean indicating if right flank aligned to reverse strand
    
    Returns:
        String indicating orientation pattern:
        - 'forward': Both flanks align forward
        - 'reverse': Both flanks align reverse (introner on reverse strand)
        - 'mixed': Mixed orientations (potential inversion/rearrangement)
        - 'partial': Only one flank available
        - 'unknown': Insufficient data
    """
    if left_reverse is None and right_reverse is None:
        return 'unknown'
    elif left_reverse is None or right_reverse is None:
        return 'partial'
    elif left_reverse == right_reverse:
        return 'reverse' if left_reverse else 'forward'
    else:
        return 'mixed'

def main():
    # Configuration
    samples = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749",
               "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629",
               "RCC692", "RCC693", "RCC833"]

    # Parse command line arguments with argparse
    parser = argparse.ArgumentParser(
        description='Build ortholog dictionary from BAM alignments',
        usage='%(prog)s <bam_dir> <output_file> [--mode {initial|rescue}] [--export-quality <output_file>]'
    )
    parser.add_argument('bam_dir', help='Directory containing BAM files')
    parser.add_argument('output_file', help='Output TSV file for ortholog results')
    parser.add_argument('--mode', choices=['initial', 'rescue'], default='rescue',
                       help='Execution mode: initial (Phase 1 only) or rescue (Phase 1 + Phase 2)')
    parser.add_argument('--export-quality', dest='quality_output', metavar='FILE',
                       help='Export alignment quality info to specified file')
    parser.add_argument('--bed-dir', dest='bed_dir', default=None,
                       help='Directory containing introner BED files (defaults to bam_dir)')
    parser.add_argument('--genome-dir', dest='genome_dir', default=None,
                       help='Directory containing genome assemblies')
    parser.add_argument('--processed-ann-dir', dest='processed_ann_dir', default=None,
                       help='Directory containing processed annotations')

    args = parser.parse_args()

    bam_dir = args.bam_dir
    bed_dir = args.bed_dir if args.bed_dir else args.bam_dir
    output_file = args.output_file
    mode = args.mode

    # Phase 2: Genome directory and processed annotations directory
    genome_dir = args.genome_dir if args.genome_dir else os.path.join(os.path.dirname(bam_dir), "assemblies")
    processed_ann_dir = args.processed_ann_dir if args.processed_ann_dir else os.path.join(bam_dir, "processed_annotations")

    # Check if we should export alignment quality info
    export_quality_info = args.quality_output is not None
    quality_output_file = args.quality_output if export_quality_info else None

    if export_quality_info:
        print(f"Will export alignment quality info to {quality_output_file}")

    print(f"Running in '{mode}' mode")
    if mode == 'initial':
        print("Phase 1 initial pass: Creating ortholog_results.initial.tsv (no rescue BAMs)")
    else:
        print("Phase 1 rescue pass: Using rescue BAMs to enhance results")

    # Load synteny context for output columns
    print("Loading synteny context for output columns...")
    synteny_context = load_introner_synteny_context(bam_dir, samples)

    # Load gene annotations (needed for Phase 2 rescue mode)
    gene_annotations = None
    if mode == 'rescue':
        print("Loading gene annotations for Phase 2 rescue...")
        gene_annotations = load_gene_annotations(processed_ann_dir, samples)

    # Build the ortholog dictionary
    if export_quality_info:
        ortholog_dict, quality_df = build_ortholog_dictionary(
            samples, bam_dir, bed_dir, genome_dir, gene_annotations,
            export_quality_info=True, mode=mode)
    else:
        ortholog_dict = build_ortholog_dictionary(
            samples, bam_dir, bed_dir, genome_dir, gene_annotations, mode=mode)
    
    # Select best matches based on metadata
    ortholog_dict = select_best_matches(ortholog_dict)
    
    # Count entries before writing
    total_entries = 0
    query_entries_by_sample = {sample: 0 for sample in samples}
    
    for query in ortholog_dict:
        for target in ortholog_dict[query]:
            entries_count = len(ortholog_dict[query][target])
            total_entries += entries_count
            query_entries_by_sample[query] = query_entries_by_sample.get(query, 0) + entries_count
    
    print(f"Total entries to write: {total_entries}")
    print("Entries by query sample:")
    for sample, count in query_entries_by_sample.items():
        print(f"  {sample}: {count}")
    
    # Collect all rows in a list for DataFrame creation
    all_rows = []
    
    for query in ortholog_dict:
        for target in ortholog_dict[query]:
            for query_id, match in ortholog_dict[query][target].items():
                # Extract metadata fields with empty defaults if not available
                query_gene = match['query_metadata'].get('gene', '') if match['query_metadata'] else ''
                target_gene = match['target_metadata'].get('gene', '') if match['target_metadata'] else ''
                
                # Extract family values and ensure they're integers if present
                query_family = match['query_metadata'].get('family', '') if match['query_metadata'] else ''
                if query_family and str(query_family).strip() and query_family != '':
                    try:
                        query_family = int(float(query_family))
                    except (ValueError, TypeError):
                        query_family = ''
                        
                target_family = match['target_metadata'].get('family', '') if match['target_metadata'] else ''
                if target_family and str(target_family).strip() and target_family != '':
                    try:
                        target_family = int(float(target_family))
                    except (ValueError, TypeError):
                        target_family = ''
                
                query_splice_site = match['query_metadata'].get('splice_site', '') if match['query_metadata'] else ''
                target_splice_site = match['target_metadata'].get('splice_site', '') if match['target_metadata'] else ''
                query_introner_id = match['query_metadata'].get('introner_id', '') if match['query_metadata'] else ''
                target_introner_id = match['target_metadata'].get('introner_id', '') if match['target_metadata'] else ''
                match_score = match.get('match_score', 0)

                # Extract orientation information
                orientation = match.get('orientation', 'unknown')
                left_reverse = match.get('left_reverse', None)
                right_reverse = match.get('right_reverse', None)

                # Extract synteny information for query
                query_contig = match.get('query_contig', '')
                query_start_temp = match.get('query_start', '')
                query_end_temp = match.get('query_end', '')
                query_unique_id = f"{query_introner_id}|{query_contig}|{query_start_temp}|{query_end_temp}"
                query_synteny_info = synteny_context.get(query, {}).get(query_unique_id, {})
                query_upstream_gene = query_synteny_info.get('upstream_gene', '')
                query_downstream_gene = query_synteny_info.get('downstream_gene', '')

                # Extract synteny information for target
                target_contig = match.get('target_contig', '')
                target_start_temp = match.get('target_start', '')
                target_end_temp = match.get('target_end', '')
                target_unique_id = f"{target_introner_id}|{target_contig}|{target_start_temp}|{target_end_temp}"
                target_synteny_info = synteny_context.get(target, {}).get(target_unique_id, {})
                target_upstream_gene = target_synteny_info.get('upstream_gene', '')
                target_downstream_gene = target_synteny_info.get('downstream_gene', '')
                
                # Convert boolean values to strings for output
                left_reverse_str = str(left_reverse) if left_reverse is not None else ''
                right_reverse_str = str(right_reverse) if right_reverse is not None else ''
                
                # Format start and end values as integers
                query_start = match['query_start']
                query_end = match['query_end']
                target_start = match['target_start']
                target_end = match['target_end']
                
                # Ensure start is always smaller than end for both query and target
                if query_start is not None and query_end is not None and query_start > query_end:
                    query_start, query_end = query_end, query_start
                    
                if target_start is not None and target_end is not None and target_start > target_end:
                    target_start, target_end = target_end, target_start
                
                # Convert numeric values to integer format for output, preserving None as empty string
                query_start_formatted = int(query_start) if query_start is not None else ''
                query_end_formatted = int(query_end) if query_end is not None else ''
                target_start_formatted = int(target_start) if target_start is not None else ''
                target_end_formatted = int(target_end) if target_end is not None else ''
                
                # Build row dictionary with explicit defaults for all columns
                row = {
                    'query': query,
                    'target': target,
                    'query_id': query_id,
                    'scenario': match['scenario'],
                    'query_contig': match['query_contig'] or '',
                    'query_start': query_start_formatted,
                    'query_end': query_end_formatted,
                    'target_contig': match['target_contig'] or '',
                    'target_start': target_start_formatted,
                    'target_end': target_end_formatted,
                    'query_gene': query_gene,
                    'target_gene': target_gene,
                    'query_family': query_family,
                    'target_family': target_family,
                    'query_splice_site': query_splice_site,
                    'target_splice_site': target_splice_site,
                    'query_introner_id': query_introner_id,
                    'target_introner_id': target_introner_id,
                    'match_score': match_score,
                    'orientation': orientation,
                    'left_reverse': left_reverse_str,
                    'right_reverse': right_reverse_str,
                    'query_upstream_gene': query_upstream_gene,
                    'query_downstream_gene': query_downstream_gene,
                    'target_upstream_gene': target_upstream_gene,
                    'target_downstream_gene': target_downstream_gene
                }
                all_rows.append(row)
    
    # Create DataFrame and write to TSV
    output_df = pd.DataFrame(all_rows)
    output_df.to_csv(output_file, sep='\t', index=False)
    rows_written = len(all_rows)
    
    print(f"Total rows written to output file: {rows_written}")
    print("Done! Results written to ortholog_results.tsv")
    
    # Export alignment quality info if available
    if export_quality_info and 'quality_df' in locals():
        quality_df.to_csv(quality_output_file, sep='\t', index=False)
        print(f"Exported alignment quality information to {quality_output_file}")
    
    return ortholog_dict

if __name__ == "__main__":
    main()
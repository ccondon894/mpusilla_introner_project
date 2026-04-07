#!/usr/bin/env python3
import argparse
import logging
from collections import defaultdict, namedtuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Named tuple to track gene information
GeneInfo = namedtuple('GeneInfo', ['transcript_id', 'gene_id', 'chrom', 'start', 'end', 'identity', 'order'])

def parse_paf_line(line):
    """Parse PAF line to extract transcript ID and alignment stats"""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None

    # PAF format: original_transcript_id is in column 2 (0-indexed: column 1)
    # Matching bases in column 11 (0-indexed: column 10)
    # Total bases in column 12 (0-indexed: column 11)
    original_transcript_id = fields[1]
    try:
        matching_bases = int(fields[10])
        total_bases = int(fields[11])
        percent_identity = (matching_bases / total_bases) * 100 if total_bases > 0 else 0
    except (ValueError, IndexError):
        logger.warning(f"Could not parse alignment stats from PAF line: {line.strip()}")
        return None

    return {
        'transcript_id': original_transcript_id,
        'percent_identity': percent_identity,
        'matching_bases': matching_bases,
        'total_bases': total_bases
    }

def check_overlap(start1, end1, start2, end2):
    """Check if two genomic ranges overlap"""
    return not (end1 < start2 or end2 < start1)

def filter_overlapping_genes(valid_groups):
    """
    Filter overlapping genes, keeping only those with highest identity.
    For ties, keep the first gene encountered.

    Returns: list of (paf_data, gene_group) tuples after overlap filtering
    """
    # Build a comprehensive gene index by contig
    genes_by_contig = defaultdict(list)
    gene_index = {}  # Map from transcript_id -> full gene info

    logger.info("Building gene coordinate index for overlap detection...")

    for order, (paf_data, gene_group) in enumerate(valid_groups):
        if not paf_data:
            continue

        transcript_id = paf_data['transcript_id']

        # Extract chromosome and coordinates from GTF lines
        chrom = None
        exon_starts = []
        exon_ends = []
        gene_id = None

        for line in gene_group:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if not chrom:
                chrom = fields[0]

            # Extract exon coordinates
            if fields[2] in ['exon', 'CDS']:
                exon_starts.append(int(fields[3]))
                exon_ends.append(int(fields[4]))

            # Extract gene_id from attributes
            if not gene_id and fields[2] == 'transcript':
                attrs = parse_gtf_attributes(fields[8])
                gene_id = attrs.get('gene_id')

        if chrom and exon_starts:
            gene_start = min(exon_starts)
            gene_end = max(exon_ends)
            identity = paf_data['percent_identity']

            gene_info = GeneInfo(
                transcript_id=transcript_id,
                gene_id=gene_id,
                chrom=chrom,
                start=gene_start,
                end=gene_end,
                identity=identity,
                order=order
            )

            genes_by_contig[chrom].append(gene_info)
            gene_index[transcript_id] = gene_info

    # Identify overlapping gene pairs and mark lower-identity genes for removal
    genes_to_remove = set()

    for chrom, genes in genes_by_contig.items():
        logger.info(f"Checking {len(genes)} genes on {chrom} for overlaps...")

        for i, gene1 in enumerate(genes):
            if gene1.transcript_id in genes_to_remove:
                continue  # Already marked for removal

            for gene2 in genes[i+1:]:
                if gene2.transcript_id in genes_to_remove:
                    continue

                # Check for overlap
                if check_overlap(gene1.start, gene1.end, gene2.start, gene2.end):
                    # Found an overlap, compare identities
                    if gene1.identity > gene2.identity:
                        genes_to_remove.add(gene2.transcript_id)
                        logger.info(f"  Overlap detected: {gene1.transcript_id} (id={gene1.identity:.1f}%, span={gene1.end-gene1.start}bp) "
                                  f"vs {gene2.transcript_id} (id={gene2.identity:.1f}%, span={gene2.end-gene2.start}bp) "
                                  f"-> Removing {gene2.transcript_id}")
                    elif gene2.identity > gene1.identity:
                        genes_to_remove.add(gene1.transcript_id)
                        logger.info(f"  Overlap detected: {gene1.transcript_id} (id={gene1.identity:.1f}%, span={gene1.end-gene1.start}bp) "
                                  f"vs {gene2.transcript_id} (id={gene2.identity:.1f}%, span={gene2.end-gene2.start}bp) "
                                  f"-> Removing {gene1.transcript_id}")
                        break  # gene1 is marked for removal, skip other comparisons
                    # else: identical identity, keep both (tie-breaking by order/first-come)

    # Filter out marked genes
    filtered_groups = []
    removed_count = 0
    for paf_data, gene_group in valid_groups:
        if paf_data and paf_data['transcript_id'] in genes_to_remove:
            removed_count += 1
            logger.info(f"Removing overlapping gene: {paf_data['transcript_id']} (id={paf_data['percent_identity']:.1f}%)")
        else:
            filtered_groups.append((paf_data, gene_group))

    logger.info(f"Overlap-based filtering removed {removed_count} genes")
    return filtered_groups

def parse_gtf_attributes(attr_str):
    """Parse GTF attribute string into a dictionary"""
    attrs = {}
    # Split by semicolon and strip whitespace
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            try:
                # Split by first space and strip quotes
                key, value = attr.strip().split(' ', 1)
                attrs[key] = value.strip('"')
            except ValueError:
                continue
    return attrs

def format_gtf_attributes(attrs):
    """Format attribute dictionary back to GTF format"""
    formatted = []
    # Ensure gene_id and transcript_id come first if they exist
    priority_attrs = ['gene_id', 'transcript_id']
    for key in priority_attrs:
        if key in attrs:
            formatted.append(f'{key} "{attrs[key]}"')

    # Add remaining attributes
    for key, value in attrs.items():
        if key not in priority_attrs:
            formatted.append(f'{key} "{value}"')

    return '; '.join(formatted) + ';'

def filter_gtf(input_gtf, output_gtf, min_identity=70.0):
    """Filter GTF based on percent identity threshold and update transcript IDs"""
    logger.info("Processing GTF file...")

    with open(input_gtf) as f:
        lines = f.readlines()

    gene_groups = []
    current_paf_data = None
    current_gene_group = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('##PAF'):
            # Process previous group if it exists
            if current_gene_group and current_paf_data:
                gene_groups.append((current_paf_data, current_gene_group))

            # Parse new PAF line and reset group
            current_paf_data = parse_paf_line(line)
            current_gene_group = []
            i += 1
            continue

        # Skip alignment lines (ATN, ATA, AAS, AQA)
        if line.startswith('##AT') or line.startswith('##AQ'):
            i += 1
            continue

        # Skip comment lines and empty lines
        if line.startswith('#') or not line.strip():
            i += 1
            continue

        # Parse GTF line
        fields = line.split('\t')
        if len(fields) != 9:
            i += 1
            continue

        # Add this GTF line to current group
        current_gene_group.append(line)
        i += 1

    # Don't forget the last group!
    if current_gene_group and current_paf_data:
        gene_groups.append((current_paf_data, current_gene_group))

    logger.info(f"Found {len(gene_groups)} total gene groups")

    # Filter groups based on identity threshold
    valid_groups = []
    for paf_data, gene_group in gene_groups:
        if paf_data and paf_data['percent_identity'] >= min_identity:
            valid_groups.append((paf_data, gene_group))
            logger.info(f"Keeping group with {paf_data['percent_identity']:.1f}% identity (transcript: {paf_data['transcript_id']})")
        else:
            logger.info(f"Filtering out group with {paf_data['percent_identity']:.1f}% identity (transcript: {paf_data['transcript_id']})")

    logger.info(f"Found {len(valid_groups)} gene groups above {min_identity}% identity threshold")

    # Apply overlap-based filtering
    logger.info("Applying overlap-based filtering...")
    valid_groups = filter_overlapping_genes(valid_groups)

    logger.info(f"After overlap filtering: {len(valid_groups)} gene groups remaining")

    # Write filtered GTF
    logger.info("Writing filtered GTF...")
    with open(output_gtf, 'w') as out:
        # Write GTF version header
        out.write("##gff-version 3\n")

        for paf_data, gene_group in valid_groups:
            original_transcript_id = paf_data['transcript_id']

            for line in gene_group:
                fields = line.strip().split('\t')
                attrs = parse_gtf_attributes(fields[8])

                # Update transcript_id for transcript, exon, and CDS lines
                if fields[2] in ['transcript', 'exon', 'CDS'] and 'transcript_id' in attrs:
                    attrs['transcript_id'] = original_transcript_id

                # Format and write the updated line
                fields[8] = format_gtf_attributes(attrs)
                out.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Filter miniprot GTF annotations and update transcript IDs')
    parser.add_argument('--input', required=True, help='Input GTF from miniprot')
    parser.add_argument('--output', required=True, help='Output filtered GTF')
    parser.add_argument('--min-identity', type=float, default=70.0,
                       help='Minimum percent identity for feature validation (default: 70.0)')

    args = parser.parse_args()

    try:
        filter_gtf(args.input, args.output, args.min_identity)
        logger.info("Filtering completed successfully")
    except Exception as e:
        logger.error(f"Error during filtering: {str(e)}")
        raise

if __name__ == '__main__':
    main()
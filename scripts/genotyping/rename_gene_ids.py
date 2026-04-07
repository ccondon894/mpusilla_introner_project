import argparse
import sys

def parse_attributes(attribute_string):
    """Parses a GTF/GFF attribute string into a dictionary."""
    attributes = {}
    for part in attribute_string.strip().split(';'):
        if not part:
            continue
        key_value = part.strip().split(' ', 1)
        if len(key_value) == 2:
            key, value = key_value
            attributes[key] = value.strip('"')
    return attributes

def format_attributes(attribute_dict):
    """Formats a dictionary of attributes back into a GTF string."""
    parts = []
    # A common order for GTF attributes
    for key in ['gene_id', 'transcript_id']:
        if key in attribute_dict:
            parts.append(f'{key} "{attribute_dict[key]}"')
    
    for key, value in attribute_dict.items():
        if key not in ['gene_id', 'transcript_id']:
            parts.append(f'{key} "{value}"')
             
    return '; '.join(parts) + ';'

def main():
    parser = argparse.ArgumentParser(
        description="Correctly renames gene_ids in a query GTF and sorts the output for SQANTI3 compatibility."
    )
    parser.add_argument('--reference', required=True, help="Path to the reference GTF file.")
    parser.add_argument('--query', required=True, help="Path to the query (miniprot) GTF file.")
    parser.add_argument('--output', required=True, help="Path for the final, sorted output GTF file.")
    args = parser.parse_args()

    # Step 1: Build the transcript_id -> gene_id map from the reference.
    print(f"Building transcript_id -> gene_id map from: {args.reference}")
    transcript_to_gene_map = {}
    try:
        with open(args.reference, 'r') as ref_file:
            for line in ref_file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                attributes = parse_attributes(parts[8])
                if 'transcript_id' in attributes and 'gene_id' in attributes:
                    transcript_id = attributes['transcript_id']
                    gene_id = attributes['gene_id']
                    if transcript_id not in transcript_to_gene_map:
                        transcript_to_gene_map[transcript_id] = gene_id
    except FileNotFoundError:
        print(f"Error: Reference file not found at {args.reference}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Map built for {len(transcript_to_gene_map)} unique transcripts.")

    # Step 2: Read all features into memory, renaming them along the way.
    print(f"Reading and processing features from: {args.query}")
    header_lines = []
    processed_features = []
    unmapped_transcripts = set()

    try:
        with open(args.query, 'r') as query_file:
            for line in query_file:
                if line.startswith('#'):
                    header_lines.append(line)
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                # Ignore 'gene' lines completely as they interfere with SQANTI3's parser.
                if parts[2] == 'gene':
                    continue

                attributes = parse_attributes(parts[8])
                query_transcript_id = attributes.get('transcript_id')

                if query_transcript_id:
                    correct_gene_id = transcript_to_gene_map.get(query_transcript_id)
                    if correct_gene_id:
                        attributes['gene_id'] = correct_gene_id
                    else:
                        unmapped_transcripts.add(query_transcript_id)
                
                # Store the processed feature data for sorting
                # We need the full line parts + the processed attributes dict
                processed_features.append((parts, attributes))

    except FileNotFoundError:
        print(f"Error: Query file not found at {args.query}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Sort the processed features to ensure correct grouping.
    print("Sorting features to ensure correct transcript grouping...")
    
    # Sort by: 1. Chromosome, 2. Transcript ID, 3. Feature Type (transcript, exon, CDS), 4. Start Coordinate
    # This ensures each transcript's features are grouped together with correct ordering for SQANTI3
    def sort_key(feature_data):
        parts, attrs = feature_data
        chrom = parts[0]
        feature_type = parts[2]
        start_pos = int(parts[3])
        transcript_id = attrs.get('transcript_id', '')

        # Define feature type priority: transcript first, then exon, then CDS
        feature_priority = {
            'transcript': 0,
            'exon': 1,
            'CDS': 2
        }
        type_priority = feature_priority.get(feature_type, 999)

        return (chrom, transcript_id, type_priority, start_pos)

    processed_features.sort(key=sort_key)

    # Step 4: Write the sorted and processed features to the output file.
    print(f"Writing {len(processed_features)} sorted features to: {args.output}")
    with open(args.output, 'w') as out_file:
        for line in header_lines:
            out_file.write(line)
        
        for parts, attributes in processed_features:
            parts[8] = format_attributes(attributes)
            out_file.write('\t'.join(parts) + '\n')

    print("\nProcessing complete.")
    if unmapped_transcripts:
        print(f"Warning: Found {len(unmapped_transcripts)} transcript IDs in the query that were not in the reference map.")

if __name__ == "__main__":
    main()
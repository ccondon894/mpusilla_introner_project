#!/usr/bin/env python3
"""
Create introner-to-gene synteny context mapping with gene neighbors.

This script:
1. Identifies which genes contain/overlap each introner (from bedtools intersect)
2. Parses upstream neighbors (from bedtools closest -id)
3. Parses downstream neighbors (from bedtools closest -iu)
4. Combines into a robust three-part syntenic fingerprint for each introner

The resulting file creates: (Upstream Neighbor) <--- (Host Gene) ---> (Downstream Neighbor)

This is an extremely robust marker for a specific genomic locus that enables
reliable ortholog identification across strains.

Usage:
    python create_introner_context.py <intersect_genes_file> <upstream_file> <downstream_file> <output_file>

Input files:
    - intersect_genes_file: Output from 'bedtools intersect -a introners -b genes -wo'
    - upstream_file: Output from 'bedtools closest -a genes -b genes -id -D ref'
    - downstream_file: Output from 'bedtools closest -a genes -b genes -iu -D ref'

Output format ({sample}.introner_synteny_map.tsv):
    introner_id  host_gene  upstream_gene  downstream_gene
"""

import sys
from collections import defaultdict


def parse_bedtools_intersect_output(intersect_file):
    """
    Parse bedtools intersect output to find genes containing introners.

    bedtools intersect -a {introner.bed} -b {genes.bed} -wo output:
    Columns 1-9: Introner BED
    Columns 10-15: Gene BED
    Column 16: Number of overlapping bp

    Returns dictionary mapping unique_introner_id -> set of containing genes
    """
    introner_genes = defaultdict(set)

    try:
        with open(intersect_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 16:
                    continue

                # Parse introner BED columns (0-8 in 0-indexed)
                introner_contig = fields[0]
                introner_start = int(fields[1])
                introner_end = int(fields[2])
                introner_id = fields[3]

                # Create unique introner identifier
                unique_introner_id = f"{introner_id}|{introner_contig}|{introner_start}|{introner_end}"

                # Parse gene BED columns: chrom(9), start(10), end(11), name(12), score(13), strand(14), overlap(15)
                gene_name = fields[12]

                introner_genes[unique_introner_id].add(gene_name)

    except Exception as e:
        print(f"Error parsing bedtools intersect output: {e}", file=sys.stderr)
        return None

    return introner_genes


def parse_upstream_neighbors(upstream_file):
    """
    Parse bedtools closest output (with -id flag) to find upstream neighbors.

    bedtools closest -id -D ref output format:
    A: gene_chrom, gene_start, gene_end, gene_name, score, strand (columns 0-5)
    B: neighbor_chrom, neighbor_start, neighbor_end, neighbor_name, score, strand (columns 6-11)
    Distance (column 12)

    Returns dictionary mapping gene_name -> upstream_neighbor_gene
    """
    neighbors = {}

    try:
        with open(upstream_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 13:
                    continue

                query_gene = fields[3]
                neighbor_gene = fields[9]

                # Skip self-matches
                if query_gene != neighbor_gene:
                    neighbors[query_gene] = neighbor_gene

    except Exception as e:
        print(f"Error parsing upstream neighbors: {e}", file=sys.stderr)
        return None

    return neighbors


def parse_downstream_neighbors(downstream_file):
    """
    Parse bedtools closest output (with -iu flag) to find downstream neighbors.

    bedtools closest -iu -D ref output format:
    A: gene_chrom, gene_start, gene_end, gene_name, score, strand (columns 0-5)
    B: neighbor_chrom, neighbor_start, neighbor_end, neighbor_name, score, strand (columns 6-11)
    Distance (column 12)

    Returns dictionary mapping gene_name -> downstream_neighbor_gene
    """
    neighbors = {}

    try:
        with open(downstream_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 13:
                    continue

                query_gene = fields[3]
                neighbor_gene = fields[9]

                # Skip self-matches
                if query_gene != neighbor_gene:
                    neighbors[query_gene] = neighbor_gene

    except Exception as e:
        print(f"Error parsing downstream neighbors: {e}", file=sys.stderr)
        return None

    return neighbors


def write_synteny_map_file(introner_genes, upstream_neighbors, downstream_neighbors, output_file):
    """
    Write introner synteny map to TSV file.

    Creates a robust syntenic fingerprint for each introner by combining:
    - Host gene (containing the introner)
    - Upstream neighbor (nearest gene before the host gene)
    - Downstream neighbor (nearest gene after the host gene)

    Output format:
        introner_id  host_gene  upstream_gene  downstream_gene

    For introners in multiple genes, creates one entry per host gene.
    """
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write('introner_id\thost_gene\tupstream_gene\tdownstream_gene\n')

            # Write one line per introner-gene pair
            entries_written = 0
            for unique_introner_id in sorted(introner_genes.keys()):
                host_genes = introner_genes[unique_introner_id]

                # For each host gene, create an entry with synteny context
                for host_gene in sorted(host_genes):
                    # Look up neighbors (default to 'NA' if not found)
                    upstream = upstream_neighbors.get(host_gene, 'NA')
                    downstream = downstream_neighbors.get(host_gene, 'NA')

                    f.write(f"{unique_introner_id}\t{host_gene}\t{upstream}\t{downstream}\n")
                    entries_written += 1

        print(f"Successfully wrote {entries_written} synteny map entries to {output_file}")
        print(f"Coverage: {len(introner_genes)} unique introner locations")
        return True

    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        return False


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <intersect_genes_file> <upstream_file> <downstream_file> <output_file>", file=sys.stderr)
        print(f"", file=sys.stderr)
        print(f"Input files:", file=sys.stderr)
        print(f"  intersect_genes_file: Output from bedtools intersect -a introners -b genes -wo", file=sys.stderr)
        print(f"  upstream_file: Output from bedtools closest -a genes -b genes -id -D ref", file=sys.stderr)
        print(f"  downstream_file: Output from bedtools closest -a genes -b genes -iu -D ref", file=sys.stderr)
        print(f"", file=sys.stderr)
        print(f"Example:", file=sys.stderr)
        print(f"  {sys.argv[0]} intersect.txt upstream.txt downstream.txt CCMP1545.introner_synteny_map.tsv", file=sys.stderr)
        sys.exit(1)

    intersect_file = sys.argv[1]
    upstream_file = sys.argv[2]
    downstream_file = sys.argv[3]
    output_file = sys.argv[4]

    print(f"Step 1: Parsing bedtools intersect output from {intersect_file}...")
    introner_genes = parse_bedtools_intersect_output(intersect_file)

    if introner_genes is None or not introner_genes:
        print("Error: No data parsed from bedtools intersect output", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(introner_genes)} introners with gene overlap information")

    print(f"Step 2: Parsing upstream neighbors from {upstream_file}...")
    upstream_neighbors = parse_upstream_neighbors(upstream_file)

    if upstream_neighbors is None:
        print("Error: No data parsed from upstream neighbors file", file=sys.stderr)
        sys.exit(1)

    print(f"Identified upstream neighbors for {len(upstream_neighbors)} genes")

    print(f"Step 3: Parsing downstream neighbors from {downstream_file}...")
    downstream_neighbors = parse_downstream_neighbors(downstream_file)

    if downstream_neighbors is None:
        print("Error: No data parsed from downstream neighbors file", file=sys.stderr)
        sys.exit(1)

    print(f"Identified downstream neighbors for {len(downstream_neighbors)} genes")

    print(f"Step 4: Creating synteny map for {len(introner_genes)} introners...")
    success = write_synteny_map_file(introner_genes, upstream_neighbors, downstream_neighbors, output_file)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()

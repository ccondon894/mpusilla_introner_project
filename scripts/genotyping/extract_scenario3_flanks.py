#!/usr/bin/env python3
"""
Extract flank sequences for scenario 3 loci that need rescue alignment.

This script:
1. Reads ortholog_results.tsv to identify scenario 3 cases
2. For each distant pair (Group1<->Group2), collects scenario 3 query IDs
3. Extracts corresponding flank sequences from existing FASTA files
4. Writes to rescue_flanks/{query}_to_{target}.scenario3_flanks.fa

Usage:
    python extract_scenario3_flanks.py <ortholog_results.tsv> <flanks_dir> <output_dir>
"""

import sys
import os
from collections import defaultdict
from Bio import SeqIO

# Sample group definitions
GROUP1_SAMPLES = {"CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482",
                  "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"}
GROUP2_SAMPLES = {"RCC1749", "RCC3052"}


def parse_ortholog_results(ortholog_file):
    """
    Parse ortholog results file to find scenario 3 cases for distant pairs.

    Returns: dict[query][target] = set of query_ids with scenario 3
    """
    scenario3_cases = defaultdict(lambda: defaultdict(set))

    with open(ortholog_file, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            query = fields[0]
            target = fields[1]
            query_id = fields[2]
            scenario = fields[3]

            # Check if this is a distant pair
            is_distant = ((query in GROUP1_SAMPLES and target in GROUP2_SAMPLES) or
                         (query in GROUP2_SAMPLES and target in GROUP1_SAMPLES))

            # Collect scenario 3 cases for distant pairs
            if is_distant and scenario == '3':
                scenario3_cases[query][target].add(query_id)

    return scenario3_cases


def extract_flank_sequences(flanks_dir, sample, query_ids):
    """
    Extract flank sequences for specific query IDs from FASTA files.

    Query IDs have format: introner_seq_XXXX_sample_contig_start_end
    FASTA headers have format: contig:flank_start-flank_end|introner_seq_XXXX_left/right

    Strategy: Match on contig and introner_seq number. The coordinates will differ
    (query_id has introner coords, FASTA has flank coords) so we don't match on those.

    Returns: dict[query_id] = {'left': seq, 'right': seq}
    """
    sequences = defaultdict(dict)

    left_flank_file = os.path.join(flanks_dir, f"{sample}.left_flanks.fa")
    right_flank_file = os.path.join(flanks_dir, f"{sample}.right_flanks.fa")

    # Build mapping from (contig, introner_num) -> query_id
    qid_map = {}
    for qid in query_ids:
        # Parse: introner_seq_XXXX_SAMPLE#X#contig_start_end
        # Example: introner_seq_2595_CCMP1545#0#scaffold_13_681659_682026
        if 'introner_seq_' in qid:
            # Extract introner number (e.g., "2595")
            introner_part = qid.split('_')[0:3]  # ['introner', 'seq', '2595']
            introner_num = introner_part[2]

            # Extract contig (everything between sample and last two underscores)
            remainder = '_'.join(qid.split('_')[3:])  # CCMP1545#0#scaffold_13_681659_682026
            contig = '_'.join(remainder.rsplit('_', 2)[0].split('_'))  # Get contig part

            key = (contig, introner_num)
            qid_map[key] = qid

    # Extract left flanks
    if os.path.exists(left_flank_file):
        for record in SeqIO.parse(left_flank_file, "fasta"):
            ref_name = record.id
            if '|introner_seq_' in ref_name:
                coord_part, id_part = ref_name.split('|')
                contig = coord_part.split(':')[0]
                introner_num = id_part.split('_')[2]

                key = (contig, introner_num)
                if key in qid_map:
                    qid = qid_map[key]
                    sequences[qid]['left'] = str(record.seq)
                    sequences[qid]['left_name'] = ref_name

    # Extract right flanks
    if os.path.exists(right_flank_file):
        for record in SeqIO.parse(right_flank_file, "fasta"):
            ref_name = record.id
            if '|introner_seq_' in ref_name:
                coord_part, id_part = ref_name.split('|')
                contig = coord_part.split(':')[0]
                introner_num = id_part.split('_')[2]

                key = (contig, introner_num)
                if key in qid_map:
                    qid = qid_map[key]
                    sequences[qid]['right'] = str(record.seq)
                    sequences[qid]['right_name'] = ref_name

    return sequences


def write_rescue_fasta(output_file, sequences):
    """
    Write rescue flank sequences to FASTA file.

    Args:
        output_file: Path to output FASTA file
        sequences: dict[query_id] = {'left': seq, 'right': seq, 'left_name': name, 'right_name': name}
    """
    with open(output_file, 'w') as f:
        for query_id, seqs in sorted(sequences.items()):
            # Write left flank (skip empty or very short sequences)
            if 'left' in seqs and 'left_name' in seqs and len(seqs['left']) >= 10:
                name = seqs['left_name']
                if not name.endswith('_left'):
                    name = f"{query_id}_left"
                f.write(f">{name}\n{seqs['left']}\n")

            # Write right flank (skip empty or very short sequences)
            if 'right' in seqs and 'right_name' in seqs and len(seqs['right']) >= 10:
                name = seqs['right_name']
                if not name.endswith('_right'):
                    name = f"{query_id}_right"
                f.write(f">{name}\n{seqs['right']}\n")


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <ortholog_results.tsv> <flanks_dir> <output_dir>")
        print(f"")
        print(f"Example:")
        print(f"  {sys.argv[0]} ortholog_results.tsv graph_blast_results graph_blast_results/rescue_flanks/")
        sys.exit(1)

    ortholog_file = sys.argv[1]
    flanks_dir = sys.argv[2]
    output_dir = sys.argv[3]

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"Parsing ortholog results from {ortholog_file}...")
    scenario3_cases = parse_ortholog_results(ortholog_file)

    total_cases = 0
    for query in scenario3_cases:
        for target in scenario3_cases[query]:
            total_cases += len(scenario3_cases[query][target])

    print(f"Found {total_cases} scenario 3 cases across distant pairs")

    # Process each distant pair
    for query in sorted(scenario3_cases.keys()):
        for target in sorted(scenario3_cases[query].keys()):
            query_ids = scenario3_cases[query][target]

            if not query_ids:
                continue

            print(f"Processing {query} → {target}: {len(query_ids)} loci")

            # Extract flank sequences
            sequences = extract_flank_sequences(flanks_dir, query, query_ids)

            # Count complete pairs (both left and right flanks)
            complete_pairs = sum(1 for qid in sequences
                               if 'left' in sequences[qid] and 'right' in sequences[qid])

            print(f"  Found {complete_pairs}/{len(query_ids)} complete flank pairs")

            # Write to rescue FASTA file
            output_file = os.path.join(output_dir, f"{query}_to_{target}.scenario3_flanks.fa")
            write_rescue_fasta(output_file, sequences)

            print(f"  Wrote rescue flanks to {output_file}")

    print("Done!")


if __name__ == '__main__':
    main()

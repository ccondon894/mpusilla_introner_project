#!/usr/bin/env python3
"""
Build MAFFT alignments for shared introner loci.

For each shared locus, gathers consensus sequences from all samples where the
introner is present, creates a combined FASTA, and runs MAFFT alignment.
Does this separately for introner body, left flank, and right flank.

Combined alignments (Group 1 + Group 2 together) are required for accurate
between-group Dxy calculation.
"""

import argparse
import json
import os
import subprocess
import sys
from Bio import SeqIO
from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Build MAFFT alignments for shared introner loci')
    parser.add_argument('--classification-json', required=True)
    parser.add_argument('--consensus-dir', required=True,
                       help='Directory containing per-sample consensus FASTAs')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for alignments')
    parser.add_argument('--regions', nargs='+',
                       default=['introner_body', 'left_flank', 'right_flank'])
    return parser.parse_args()


def load_consensus_sequences(consensus_dir, classification, region):
    """Load consensus sequences for all shared loci from all samples."""
    # Collect all samples that have sequences
    all_samples = set()
    for locus_info in classification['loci'].values():
        all_samples.update(locus_info['group1_present'])
        all_samples.update(locus_info['group2_present'])

    # Load sequences from each sample's consensus FASTA
    # Index by ortholog_id -> sample -> sequence
    locus_sequences = defaultdict(dict)

    for sample in sorted(all_samples):
        consensus_file = os.path.join(consensus_dir, f"{sample}.shared.{region}.consensus.fa")
        if not os.path.exists(consensus_file):
            print(f"  Warning: Missing consensus file: {consensus_file}", file=sys.stderr)
            continue

        for record in SeqIO.parse(consensus_file, 'fasta'):
            # Header format from samtools faidx: contig:start-end
            # We need to match this back to ortholog_id via the BED file
            # The BED 4th column has the ortholog_id, but consensus headers
            # only have the coordinates. We'll map them below.
            header = record.id
            locus_sequences[sample][header] = record

    return locus_sequences, all_samples


def load_bed_mapping(consensus_dir, sample, region):
    """Load BED file to create coordinate -> ortholog_id mapping."""
    bed_file = os.path.join(consensus_dir, f"{sample}.shared.{region}.bed")
    mapping = {}
    if os.path.exists(bed_file):
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                contig, start, end, oid = fields[0], fields[1], fields[2], fields[3]
                coord_key = f"{contig}:{start}-{end}"
                mapping[coord_key] = oid
    return mapping


def build_alignments(consensus_dir, classification, region, output_dir):
    """Build MAFFT alignments for all shared loci for one region type."""
    group1_set = set(classification['summary']['group1_samples'])
    group2_set = set(classification['summary']['group2_samples'])

    all_samples = set()
    for locus_info in classification['loci'].values():
        all_samples.update(locus_info['group1_present'])
        all_samples.update(locus_info['group2_present'])

    # Load BED mappings and consensus sequences per sample
    # ortholog_id -> list of (sample, sequence_str)
    locus_seqs = defaultdict(list)

    for sample in sorted(all_samples):
        bed_mapping = load_bed_mapping(consensus_dir, sample, region)
        consensus_file = os.path.join(consensus_dir, f"{sample}.shared.{region}.consensus.fa")

        if not os.path.exists(consensus_file):
            print(f"  Warning: Missing {consensus_file}", file=sys.stderr)
            continue

        for record in SeqIO.parse(consensus_file, 'fasta'):
            coord_key = record.id
            oid = bed_mapping.get(coord_key)
            if oid is None:
                # Try without the ::ortholog_id suffix that bedtools -name adds
                # Also handle potential format differences
                continue
            locus_seqs[oid].append((sample, str(record.seq)))

    # Build alignment for each locus
    os.makedirs(output_dir, exist_ok=True)
    n_aligned = 0
    n_skipped = 0

    for oid in sorted(classification['loci'].keys()):
        seqs = locus_seqs.get(oid, [])

        if len(seqs) < 2:
            n_skipped += 1
            continue

        # Write unaligned FASTA
        unaligned_path = os.path.join(output_dir, f"{oid}.{region}.fa")
        aligned_path = os.path.join(output_dir, f"{oid}.{region}.mafft.fa")

        with open(unaligned_path, 'w') as f:
            for sample, seq in seqs:
                f.write(f">{sample}\n{seq}\n")

        # Run MAFFT
        try:
            result = subprocess.run(
                ['mafft', '--adjustdirection', '--maxiterate', '1000',
                 '--globalpair', '--quiet', unaligned_path],
                capture_output=True, text=True, timeout=120
            )
            if result.returncode == 0:
                with open(aligned_path, 'w') as f:
                    f.write(result.stdout)
                n_aligned += 1
            else:
                print(f"  MAFFT failed for {oid} {region}: {result.stderr[:200]}",
                      file=sys.stderr)
                n_skipped += 1
        except subprocess.TimeoutExpired:
            print(f"  MAFFT timeout for {oid} {region}", file=sys.stderr)
            n_skipped += 1

        # Clean up unaligned file
        if os.path.exists(aligned_path):
            os.remove(unaligned_path)

    print(f"  {region}: {n_aligned} aligned, {n_skipped} skipped")
    return n_aligned


def main():
    args = parse_arguments()

    print("Loading classification...")
    with open(args.classification_json) as f:
        classification = json.load(f)

    n_loci = len(classification['loci'])
    print(f"Processing {n_loci} shared loci")

    for region in args.regions:
        print(f"\nBuilding {region} alignments...")
        build_alignments(args.consensus_dir, classification, region, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()

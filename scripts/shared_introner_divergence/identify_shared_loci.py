#!/usr/bin/env python3
"""
Identify introner loci shared between Group 1 and Group 2, classify them into
three categories, and generate BED files for introner body and flanking regions.

Categories:
  - all_shared: presence==1 in at least one Group 1 AND one Group 2 sample
  - fixed_shared: presence==1 in ALL Group 1 AND ALL Group 2 samples
  - polymorphic_shared: shared but not fixed (all_shared minus fixed_shared)

Presence encoding:
  1 = present, 2 = absent, 3 = not callable

Genotype matrix coordinates include 100bp flanking on each side.
  Introner body: start+100 to end-100
  Left flank: orientation-aware, 100bp upstream of introner
  Right flank: orientation-aware, 100bp downstream of introner
"""

import argparse
import json
import os
import sys
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Identify shared introner loci and generate BED files')
    parser.add_argument('--genotype-matrix', required=True)
    parser.add_argument('--group1-samples', nargs='+', required=True)
    parser.add_argument('--group2-samples', nargs='+', required=True)
    parser.add_argument('--flank-length', type=int, default=100)
    parser.add_argument('--output-dir', required=True,
                       help='Directory for per-sample BED files')
    parser.add_argument('--classification-json', required=True,
                       help='Output JSON with locus classifications')
    return parser.parse_args()


def classify_shared_loci(df, group1_samples, group2_samples):
    """Classify shared loci into three categories."""
    group1_set = set(group1_samples)
    group2_set = set(group2_samples)

    classification = {}

    for ortholog_id, group in df.groupby('ortholog_id'):
        g1 = group[group['sample'].isin(group1_set)]
        g2 = group[group['sample'].isin(group2_set)]

        g1_present = set(g1[g1['presence'] == 1]['sample'])
        g2_present = set(g2[g2['presence'] == 1]['sample'])

        if not g1_present or not g2_present:
            continue  # Not shared

        # Determine category
        all_g1_present = g1_present == group1_set
        all_g2_present = g2_present == group2_set

        if all_g1_present and all_g2_present:
            category = "fixed_shared"
        else:
            category = "polymorphic_shared"

        classification[ortholog_id] = {
            'category': category,
            'group1_present': sorted(g1_present),
            'group2_present': sorted(g2_present),
            'group1_present_count': len(g1_present),
            'group2_present_count': len(g2_present),
        }

    return classification


def write_bed_files(df, classification, group1_samples, group2_samples,
                    flank_length, output_dir):
    """Write per-sample BED files for introner body, left flank, right flank."""
    all_samples = group1_samples + group2_samples

    # Collect loci per sample
    shared_ids = set(classification.keys())

    for sample in all_samples:
        # Get rows for this sample where introner is present at shared loci
        sample_df = df[(df['sample'] == sample) &
                       (df['presence'] == 1) &
                       (df['ortholog_id'].isin(shared_ids))]

        body_records = []
        left_records = []
        right_records = []

        for _, row in sample_df.iterrows():
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            oid = row['ortholog_id']
            orientation = row['orientation']

            # Introner body: remove 100bp flanking from genotype matrix coords
            body_start = start + flank_length
            body_end = end - flank_length

            if body_end <= body_start:
                continue  # Skip degenerate loci

            body_records.append((contig, body_start, body_end, oid))

            # Flanking regions (orientation-aware, matching rule 10 logic)
            if orientation == 'reverse':
                # Biological left (5') flank is upstream from sequence end
                left_start = end
                left_end = end + flank_length
                # Biological right (3') flank is downstream from sequence start
                right_start = start - flank_length
                right_end = start
            else:
                # Biological left (5') flank is upstream from sequence start
                left_start = start - flank_length
                left_end = start
                # Biological right (3') flank is downstream from sequence end
                right_start = end
                right_end = end + flank_length

            # Skip if flanking coordinates go negative
            if left_start < 0 or right_start < 0:
                continue

            left_records.append((contig, left_start, left_end, oid))
            right_records.append((contig, right_start, right_end, oid))

        # Write BED files
        for region, records in [('introner_body', body_records),
                                 ('left_flank', left_records),
                                 ('right_flank', right_records)]:
            bed_path = os.path.join(output_dir, f"{sample}.shared.{region}.bed")
            with open(bed_path, 'w') as f:
                for contig, s, e, oid in sorted(records, key=lambda x: (x[0], x[1])):
                    f.write(f"{contig}\t{s}\t{e}\t{oid}\n")

            print(f"  {sample} {region}: {len(records)} loci")


def main():
    args = parse_arguments()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading genotype matrix...")
    df = pd.read_csv(args.genotype_matrix, sep='\t')

    print("Classifying shared loci...")
    classification = classify_shared_loci(df, args.group1_samples, args.group2_samples)

    # Add all_shared as a superset category in the JSON
    n_fixed = sum(1 for v in classification.values() if v['category'] == 'fixed_shared')
    n_poly = sum(1 for v in classification.values() if v['category'] == 'polymorphic_shared')
    n_all = len(classification)

    print(f"\nShared loci classification:")
    print(f"  all_shared: {n_all}")
    print(f"  fixed_shared: {n_fixed}")
    print(f"  polymorphic_shared: {n_poly}")

    # Save classification JSON
    output_data = {
        'summary': {
            'all_shared': n_all,
            'fixed_shared': n_fixed,
            'polymorphic_shared': n_poly,
            'group1_samples': args.group1_samples,
            'group2_samples': args.group2_samples,
            'flank_length': args.flank_length,
        },
        'loci': classification
    }

    with open(args.classification_json, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"\nClassification saved to: {args.classification_json}")

    print("\nGenerating BED files...")
    write_bed_files(df, classification, args.group1_samples, args.group2_samples,
                    args.flank_length, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()

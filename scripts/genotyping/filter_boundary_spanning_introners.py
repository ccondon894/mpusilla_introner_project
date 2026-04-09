#!/usr/bin/env python3
"""
Filter introner candidates by gene-boundary overlap.

A biologically valid introner should be either:
  1. Fully inside a single gene (exonic, intronic, or UTR) — plausibly an
     intron-like element inserted into the gene, OR
  2. Fully in an intergenic region — plausibly a transposed introner element
     that lost its host gene or inserted into non-coding space.

Candidates whose body straddles a gene boundary (partially inside, partially
outside) or spans multiple distinct genes are structurally suspicious and
most likely represent false positive BLAST hits where a sequence with
moderate similarity to an introner reference happens to align across a
gene boundary.

This script reads the validated candidate FASTA (output of
validate_introners_blast.py) and the GTF annotation, classifies each
candidate's body by its relationship to gene boundaries, and writes a
filtered FASTA containing only fully-intra-gene or fully-intergenic
candidates. Each kept candidate has a new 'genomic_context' tag added
to its header (intra_gene or intergenic).

Boundary overhangs up to --tolerance bp are treated as extraction noise
and not considered boundary violations.

A log file records every rejected candidate with its failure reason.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict

from Bio import SeqIO


DEFAULT_FLANK_LENGTH = 100
DEFAULT_TOLERANCE = 5


def load_gene_transcripts(gtf_path):
    """Load transcript spans from GTF, grouped by contig.

    A gene's span is the min(start)..max(end) of all its transcripts on a
    single contig. We use 'transcript' features rather than 'gene' because
    miniprot GTFs sometimes lack gene lines.

    Returns a dict: contig -> sorted list of (start_0based, end_0based_excl, gene_id)
    """
    gene_spans = defaultdict(lambda: defaultdict(lambda: {'start': None, 'end': None}))

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue

            contig = fields[0]
            start = int(fields[3])  # 1-based inclusive
            end = int(fields[4])    # 1-based inclusive

            gene_id = None
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break
            if not gene_id:
                continue

            slot = gene_spans[contig][gene_id]
            if slot['start'] is None or start < slot['start']:
                slot['start'] = start
            if slot['end'] is None or end > slot['end']:
                slot['end'] = end

    # Convert to per-contig sorted list of (start_0based, end_0based_excl, gene_id)
    out = {}
    for contig, gene_dict in gene_spans.items():
        entries = []
        for gid, span in gene_dict.items():
            # Convert GTF 1-based inclusive to 0-based half-open
            entries.append((span['start'] - 1, span['end'], gid))
        entries.sort()
        out[contig] = entries
    return out


def find_overlapping_genes(body_start, body_end, gene_list):
    """Find all genes whose span overlaps the given body interval.

    body_start, body_end: 0-based half-open
    gene_list: sorted list of (start, end, gene_id) tuples (0-based half-open)
    """
    overlapping = []
    for gstart, gend, gid in gene_list:
        # gene is past the body — since list is sorted by start, we can break
        if gstart >= body_end:
            break
        if gend > body_start:
            overlapping.append((gstart, gend, gid))
    return overlapping


def classify_body(body_start, body_end, gene_list, tolerance):
    """Classify a body's position relative to genes on its contig.

    Returns a tuple (status, context, detail) where:
      - status: 'keep' or 'reject'
      - context: 'intra_gene' | 'intergenic' | None (when rejected)
      - detail: additional info (gene_id or rejection reason)
    """
    if gene_list is None:
        return ('keep', 'intergenic', 'no_genes_on_contig')

    overlaps = find_overlapping_genes(body_start, body_end, gene_list)

    if not overlaps:
        return ('keep', 'intergenic', 'no_overlap')

    if len(overlaps) > 1:
        # Body spans multiple genes
        gids = ';'.join(g[2] for g in overlaps)
        return ('reject', None, f'spans_multiple_genes:{gids}')

    # Single overlapping gene — check overhang
    gstart, gend, gid = overlaps[0]
    overhang_left = max(0, gstart - body_start)
    overhang_right = max(0, body_end - gend)

    if overhang_left > tolerance or overhang_right > tolerance:
        return ('reject', None,
                f'boundary_overhang:gene={gid};left={overhang_left};right={overhang_right}')

    return ('keep', 'intra_gene', gid)


def main():
    parser = argparse.ArgumentParser(
        description='Filter introner candidates by gene-boundary overlap')
    parser.add_argument('--input', required=True,
                        help='Input validated candidate FASTA '
                             '(output of validate_introners_blast.py)')
    parser.add_argument('--gtf', required=True,
                        help='GTF annotation file')
    parser.add_argument('--output', required=True,
                        help='Output filtered FASTA')
    parser.add_argument('--log', required=True,
                        help='Output log file of rejected candidates')
    parser.add_argument('--flanking-length', type=int, default=DEFAULT_FLANK_LENGTH,
                        help=f'Flanking length to trim for body extraction '
                             f'(default: {DEFAULT_FLANK_LENGTH})')
    parser.add_argument('--tolerance', type=int, default=DEFAULT_TOLERANCE,
                        help=f'Boundary overhang tolerance in bp '
                             f'(default: {DEFAULT_TOLERANCE})')
    args = parser.parse_args()

    # Load gene spans
    print(f"Loading gene transcripts from {args.gtf}...")
    gene_spans = load_gene_transcripts(args.gtf)
    total_genes = sum(len(v) for v in gene_spans.values())
    print(f"  Loaded {total_genes} gene transcript spans on {len(gene_spans)} contigs")

    stats = Counter()
    rejections = []
    kept_records = []

    print(f"Processing candidates from {args.input}...")
    for record in SeqIO.parse(args.input, 'fasta'):
        stats['total'] += 1

        header = record.id + (' ' + record.description.split(' ', 1)[1]
                              if ' ' in record.description else '')
        # Parse the contig and body coordinates from the header
        # Header format: CONTIG:START-END|introner_seq_ID family=... ...
        try:
            header_first = record.id
            region_and_id = header_first.split('|')[0]
            contig, coords = region_and_id.split(':', 1)
            fasta_start_str, fasta_end_str = coords.split('-')
            fasta_start = int(fasta_start_str)
            fasta_end = int(fasta_end_str)
        except (ValueError, IndexError):
            print(f"  Warning: could not parse header: {record.id}", file=sys.stderr)
            stats['parse_error'] += 1
            continue

        # Body = FASTA span minus flanking on each side
        # The FASTA coordinates are the full locus span (body + 2*flanks)
        # so body spans [fasta_start + flank, fasta_end - flank) in 0-based.
        body_start = fasta_start + args.flanking_length
        body_end = fasta_end - args.flanking_length
        if body_end <= body_start:
            stats['body_too_short'] += 1
            rejections.append({
                'introner_id': record.id,
                'reason': 'body_too_short',
                'detail': f'{body_start}..{body_end}',
            })
            continue

        gene_list = gene_spans.get(contig)
        status, context, detail = classify_body(
            body_start, body_end, gene_list, args.tolerance)

        if status == 'reject':
            stats['rejected'] += 1
            reason = detail.split(':', 1)[0]
            stats[f'reject_{reason}'] += 1
            rejections.append({
                'introner_id': record.id,
                'reason': reason,
                'detail': detail,
            })
            continue

        # Keep: add genomic_context tag to the description
        stats['kept'] += 1
        stats[f'kept_{context}'] += 1

        # Rebuild header with genomic_context tag
        # Parse existing description fields
        orig_desc = record.description
        new_desc = f"{orig_desc} genomic_context={context}"
        record.description = new_desc
        # Biopython doesn't always respect description changes; set id and
        # re-output manually to be safe
        kept_records.append(record)

    # Write filtered FASTA
    with open(args.output, 'w') as f:
        for record in kept_records:
            f.write(f">{record.description}\n{str(record.seq)}\n")

    # Write rejection log
    with open(args.log, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['introner_id', 'reason', 'detail'],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rejections)

    # Print summary
    print(f"\n=== Filtering summary ===")
    print(f"  Total candidates:           {stats['total']}")
    print(f"  Kept:                       {stats['kept']}")
    print(f"    intra_gene:               {stats['kept_intra_gene']}")
    print(f"    intergenic:               {stats['kept_intergenic']}")
    print(f"  Rejected:                   {stats['rejected']}")
    for key in sorted(stats):
        if key.startswith('reject_') and stats[key] > 0:
            short = key.replace('reject_', '')
            print(f"    {short:30s}  {stats[key]}")
    if stats['body_too_short']:
        print(f"  Body too short:             {stats['body_too_short']}")
    if stats['parse_error']:
        print(f"  Header parse errors:        {stats['parse_error']}")

    print(f"\nWrote {len(kept_records)} records to {args.output}")
    print(f"Wrote {len(rejections)} rejections to {args.log}")


if __name__ == '__main__':
    main()

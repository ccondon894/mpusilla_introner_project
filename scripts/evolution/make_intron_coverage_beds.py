#!/usr/bin/env python3
"""
Create per-sample BED files of intron loci with 100bp flanking for the coverage caller.

Reads a reference catalog of intron loci and a sample GTF, extracts introns for genes
present in the catalog, adds 100bp flanking on each side, and writes a BED file.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pysam


# ── GTF parsing (copied from non_introner_intron_afs.py) ─────────────────────

def parse_gtf_attribute(attr_str, key="transcript_id"):
    """Extract a value from a GTF attribute string."""
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if field.startswith(key):
            return field.split('"')[1]
    return None


def extract_introns_from_gtf(gtf_path):
    """
    Parse a GTF and extract introns for each gene.

    Returns dict: gene_id -> list of (contig, start, end, intron_index)
    where intron_index is 0-based ordinal in coding order
    (plus strand: genomic order; minus strand: reversed).
    """
    gene_exons = defaultdict(list)
    gene_strand = {}
    gene_contig = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "exon":
                continue

            contig = parts[0]
            start = int(parts[3])  # 1-based inclusive
            end = int(parts[4])    # 1-based inclusive
            strand = parts[6]
            gene_id = parse_gtf_attribute(parts[8])
            if gene_id is None:
                continue

            gene_exons[gene_id].append((start, end))
            gene_strand[gene_id] = strand
            gene_contig[gene_id] = contig

    introns = {}
    for gene_id, exons in gene_exons.items():
        if len(exons) < 2:
            continue

        strand = gene_strand[gene_id]
        contig = gene_contig[gene_id]

        # Sort exons by genomic position
        exons_sorted = sorted(exons, key=lambda x: x[0])

        # Compute intron gaps in genomic order
        genomic_introns = []
        for i in range(len(exons_sorted) - 1):
            intron_start = exons_sorted[i][1] + 1
            intron_end = exons_sorted[i + 1][0] - 1
            genomic_introns.append((intron_start, intron_end))

        # Assign coding-order indices
        if strand == "-":
            # Minus strand: coding order is reversed genomic order
            result = []
            for idx, (istart, iend) in enumerate(reversed(genomic_introns)):
                result.append((contig, istart, iend, idx))
        else:
            result = []
            for idx, (istart, iend) in enumerate(genomic_introns):
                result.append((contig, istart, iend, idx))

        introns[gene_id] = result

    return introns


# ── Main logic ────────────────────────────────────────────────────────────────

def load_reference_catalog(catalog_path):
    """
    Load reference catalog TSV (columns: gene_id, intron_idx, contig, start, end).
    Returns set of gene_ids present in the catalog.
    """
    gene_ids = set()
    with open(catalog_path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gene_ids.add(parts[0])
    return gene_ids


def get_contig_lengths(assembly_path):
    """Get contig lengths from assembly FASTA using pysam."""
    fa = pysam.FastaFile(assembly_path)
    lengths = dict(zip(fa.references, fa.lengths))
    fa.close()
    return lengths


def main():
    parser = argparse.ArgumentParser(
        description="Create per-sample BED files of intron loci with 100bp flanking"
    )
    parser.add_argument("--gtf", required=True, help="Sample GTF file")
    parser.add_argument("--reference_catalog", required=True,
                        help="Reference catalog TSV (gene_id, intron_idx, contig, start, end)")
    parser.add_argument("--assembly", required=True, help="Sample assembly FASTA")
    parser.add_argument("--sample", required=True, help="Sample name (for logging)")
    parser.add_argument("--output", required=True, help="Output BED file")
    args = parser.parse_args()

    # Step 1: Load reference catalog gene_ids
    catalog_genes = load_reference_catalog(args.reference_catalog)
    print(f"[{args.sample}] Loaded {len(catalog_genes)} genes from reference catalog")

    # Step 2: Extract introns from sample GTF (only genes in catalog)
    all_introns = extract_introns_from_gtf(args.gtf)
    introns = {g: v for g, v in all_introns.items() if g in catalog_genes}
    n_introns = sum(len(v) for v in introns.values())
    print(f"[{args.sample}] Extracted {n_introns} introns in {len(introns)} catalog genes")

    # Step 3: Get contig lengths for capping bed_end
    contig_lengths = get_contig_lengths(args.assembly)

    # Step 4: Build BED records with 100bp flanking
    FLANK = 100
    MIN_INTRON_BODY = 10  # skip introns shorter than this (GTF coords)

    bed_records = []
    n_skipped_short = 0
    n_skipped_contig = 0

    for gene_id, intron_list in sorted(introns.items()):
        for contig, intron_start, intron_end, intron_idx in intron_list:
            # intron_start, intron_end are 1-based inclusive (GTF)
            intron_body_len = intron_end - intron_start + 1
            if intron_body_len < MIN_INTRON_BODY:
                n_skipped_short += 1
                continue

            if contig not in contig_lengths:
                n_skipped_contig += 1
                continue

            # Convert to 0-based half-open BED coords with flanking
            bed_start = max(0, intron_start - 1 - FLANK)
            bed_end = min(intron_end + FLANK, contig_lengths[contig])

            intron_id = f"{gene_id}__intron_{intron_idx}"
            bed_records.append((contig, bed_start, bed_end, intron_id, gene_id))

    if n_skipped_short > 0:
        print(f"[{args.sample}] Skipped {n_skipped_short} introns shorter than {MIN_INTRON_BODY}bp")
    if n_skipped_contig > 0:
        print(f"[{args.sample}] Skipped {n_skipped_contig} introns on missing contigs")
    print(f"[{args.sample}] Writing {len(bed_records)} BED records")

    # Step 5: Write BED file
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        f.write("contig\tstart\tend\tintron_id\tgene\tpresence\n")
        for contig, start, end, intron_id, gene_id in bed_records:
            f.write(f"{contig}\t{start}\t{end}\t{intron_id}\t{gene_id}\t1\n")

    print(f"[{args.sample}] Done: {args.output}")


if __name__ == "__main__":
    main()

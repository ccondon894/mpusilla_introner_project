#!/usr/bin/env python3
"""
Build a reference non-introner intron catalog.

Extracts introns from the reference GTF, removes those overlapping introner loci,
and writes a TSV catalog of the remaining (non-introner) introns.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


# ── Intron extraction ─────────────────────────────────────────────────────────

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


# ── Introner filtering ────────────────────────────────────────────────────────

def load_introner_loci(bed_path):
    """
    Load introner loci BED file and return set of (contig, body_start, body_end).
    The BED coordinates include 100bp flanking on each side, so the actual
    introner body is (start + 100, end - 100).
    """
    loci = []
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig = parts[0]
            bed_start = int(parts[1])  # 0-based
            bed_end = int(parts[2])    # 0-based exclusive
            # Strip 100bp flanking on each side to get introner body
            body_start = bed_start + 100
            body_end = bed_end - 100
            loci.append((contig, body_start, body_end))
    return loci


def overlaps(intron_start, intron_end, locus_start, locus_end):
    """Check if two intervals overlap. intron coords are 1-based inclusive,
    locus coords are 0-based half-open from BED."""
    # Convert intron to 0-based half-open for comparison
    return intron_start - 1 < locus_end and locus_start < intron_end


def filter_introner_introns(ref_introns, introner_loci):
    """
    Remove reference introns that overlap with introner loci.
    Returns filtered dict of gene_id -> list of (contig, start, end, intron_index).
    """
    loci_by_contig = defaultdict(list)
    for contig, body_start, body_end in introner_loci:
        loci_by_contig[contig].append((body_start, body_end))

    filtered = {}
    n_removed = 0
    n_kept = 0

    for gene_id, intron_list in ref_introns.items():
        kept = []
        for contig, istart, iend, intron_idx in intron_list:
            is_introner = False
            for lstart, lend in loci_by_contig.get(contig, []):
                if overlaps(istart, iend, lstart, lend):
                    is_introner = True
                    break
            if is_introner:
                n_removed += 1
            else:
                kept.append((contig, istart, iend, intron_idx))
                n_kept += 1
        if kept:
            filtered[gene_id] = kept

    return filtered, n_removed, n_kept


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build reference non-introner intron catalog"
    )
    parser.add_argument("--reference_gtf", required=True,
                        help="Reference GTF file")
    parser.add_argument("--introner_loci", required=True,
                        help="BED file of introner loci (with 100bp flanking)")
    parser.add_argument("--include_mating_region", action="store_true",
                        help="Keep introns in the CCMP1545 mating-type region")
    parser.add_argument("--output", required=True,
                        help="Output TSV catalog")
    args = parser.parse_args()

    # Step 1: Extract introns from reference GTF
    print("Extracting introns from reference GTF...")
    ref_introns = extract_introns_from_gtf(args.reference_gtf)
    n_total = sum(len(v) for v in ref_introns.values())
    print(f"  {len(ref_introns)} genes with introns, {n_total} total introns")

    # Step 2: Load introner loci
    introner_loci = load_introner_loci(args.introner_loci)
    print(f"  Loaded {len(introner_loci)} introner loci")

    # Step 3: Filter out introner introns
    nonintroner, n_removed, n_kept = filter_introner_introns(ref_introns, introner_loci)
    print(f"  Removed {n_removed} introner-overlapping introns")
    print(f"  Retained {n_kept} non-introner introns in {len(nonintroner)} genes")

    # Step 3b: Filter out genes in the mating-type region unless explicitly kept
    MATING_CHROM = "CCMP1545#0#scaffold_2"
    MATING_START = 49808
    MATING_END = 1730591

    if args.include_mating_region:
        n_kept_final = sum(len(v) for v in nonintroner.values())
        print("  Kept mating-type region introns")
        print(f"  Final catalog: {n_kept_final} introns in {len(nonintroner)} genes")
    else:
        n_mating_genes = 0
        n_mating_introns = 0
        filtered = {}
        for gene_id, intron_list in nonintroner.items():
            # Check if any intron in this gene falls within the mating-type region
            in_mating = any(
                contig == MATING_CHROM and istart <= MATING_END and iend >= MATING_START
                for contig, istart, iend, _ in intron_list
            )
            if in_mating:
                n_mating_genes += 1
                n_mating_introns += len(intron_list)
            else:
                filtered[gene_id] = intron_list
        nonintroner = filtered
        n_kept_final = sum(len(v) for v in nonintroner.values())
        print(f"  Removed {n_mating_introns} introns in {n_mating_genes} mating-type region genes")
        print(f"  Final catalog: {n_kept_final} introns in {len(nonintroner)} genes")

    # Step 4: Write output TSV
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        f.write("gene_id\tintron_idx\tcontig\tstart\tend\n")
        for gene_id in sorted(nonintroner):
            for contig, start, end, intron_idx in nonintroner[gene_id]:
                f.write(f"{gene_id}\t{intron_idx}\t{contig}\t{start}\t{end}\n")

    print(f"  Catalog written to: {args.output}")
    print("Done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate Circos histogram tracks from current introner and non-introner intron
genotype matrices.

The introner genotype matrix is long format with per-sample coordinates.
The non-introner intron genotype matrix is wide format with reference
coordinates, so non-reference sample coordinates are recovered from that
sample's GTF using gene_id + intron_index.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


PRESENT_CALLS = {"1", "1.0", "present", "PRESENT", "Present"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--introner-matrix", required=True,
                        help="Long-format genotype_matrix.final.tsv")
    parser.add_argument("--non-introner-matrix", required=True,
                        help="Wide-format intron_genotype_matrix.tsv")
    parser.add_argument("--karyotype", required=True,
                        help="Filtered Circos karyotype file")
    parser.add_argument("--strain1-name", required=True)
    parser.add_argument("--strain2-name", required=True)
    parser.add_argument("--non-introner-reference", required=True,
                        help="Sample represented by ref_start/ref_end columns")
    parser.add_argument("--strain1-gtf", required=True,
                        help="GTF for strain1")
    parser.add_argument("--strain2-gtf", required=True,
                        help="GTF for strain2")
    parser.add_argument("--bin-size", type=int, default=50000)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def circos_contig_name(contig):
    return contig.replace("#", "_")


def load_karyotype(path, strain_names):
    chromosomes = {strain: {} for strain in strain_names}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6 or parts[0] != "chr":
                continue
            chrom = parts[2]
            length = int(parts[5])
            for strain in strain_names:
                if chrom.startswith(f"{strain}_"):
                    chromosomes[strain][chrom] = length
                    break
    return chromosomes


def make_empty_bins(chrom_lengths, bin_size):
    bins = {}
    for chrom, length in chrom_lengths.items():
        n_bins = (length + bin_size - 1) // bin_size
        bins[chrom] = [0] * n_bins
    return bins


def add_feature(bins, bin_size, contig, start, end):
    chrom = circos_contig_name(contig)
    if chrom not in bins:
        return False
    if start <= 0 or end <= 0:
        return False
    if end < start:
        start, end = end, start

    midpoint = (start + end) // 2
    midpoint = max(midpoint, 1)
    idx = (midpoint - 1) // bin_size
    if idx >= len(bins[chrom]):
        idx = len(bins[chrom]) - 1
    bins[chrom][idx] += 1
    return True


def parse_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def is_present(value):
    return str(value).strip() in PRESENT_CALLS


def parse_gtf_attribute(attr_str, key="transcript_id"):
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if field.startswith(key):
            parts = field.split('"')
            if len(parts) >= 2:
                return parts[1]
    return None


def extract_introns_from_gtf(gtf_path):
    gene_exons = defaultdict(list)
    gene_strand = {}
    gene_contig = {}

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue

            contig = parts[0]
            start = parse_int(parts[3])
            end = parse_int(parts[4])
            if start is None or end is None:
                continue
            gene_id = parse_gtf_attribute(parts[8])
            if gene_id is None:
                continue

            gene_exons[gene_id].append((start, end))
            gene_strand[gene_id] = parts[6]
            gene_contig[gene_id] = contig

    introns = {}
    for gene_id, exons in gene_exons.items():
        if len(exons) < 2:
            continue
        exons_sorted = sorted(exons, key=lambda x: x[0])
        genomic_introns = []
        for idx in range(len(exons_sorted) - 1):
            intron_start = exons_sorted[idx][1] + 1
            intron_end = exons_sorted[idx + 1][0] - 1
            if intron_end >= intron_start:
                genomic_introns.append((intron_start, intron_end))

        if gene_strand.get(gene_id) == "-":
            indexed = enumerate(reversed(genomic_introns))
        else:
            indexed = enumerate(genomic_introns)

        contig = gene_contig[gene_id]
        for intron_idx, (start, end) in indexed:
            introns[(gene_id, intron_idx)] = (contig, start, end)

    return introns


def count_introner_matrix(matrix_path, sample, bins, bin_size):
    seen = counted = skipped = 0
    with open(matrix_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "presence", "contig", "start", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"{matrix_path} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            if row["sample"] != sample:
                continue
            seen += 1
            if not is_present(row["presence"]):
                continue
            start = parse_int(row["start"])
            end = parse_int(row["end"])
            if start is None or end is None or not row["contig"]:
                skipped += 1
                continue
            if add_feature(bins, bin_size, row["contig"], start, end):
                counted += 1
            else:
                skipped += 1

    return {"seen": seen, "counted": counted, "skipped": skipped}


def count_non_introner_matrix(matrix_path, sample, reference_sample,
                              sample_introns, bins, bin_size):
    seen = counted = skipped = absent_or_missing = 0
    with open(matrix_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gene_id", "contig", "ref_start", "ref_end",
            "intron_index", sample
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"{matrix_path} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            seen += 1
            if not is_present(row[sample]):
                absent_or_missing += 1
                continue

            if sample == reference_sample:
                contig = row["contig"]
                start = parse_int(row["ref_start"])
                end = parse_int(row["ref_end"])
            else:
                intron_idx = parse_int(row["intron_index"])
                key = (row["gene_id"], intron_idx)
                coords = sample_introns.get(key)
                if coords is None:
                    skipped += 1
                    continue
                contig, start, end = coords

            if start is None or end is None:
                skipped += 1
                continue
            if add_feature(bins, bin_size, contig, start, end):
                counted += 1
            else:
                skipped += 1

    return {
        "seen": seen,
        "counted": counted,
        "skipped": skipped,
        "absent_or_missing": absent_or_missing,
    }


def write_histogram(path, chrom_lengths, bins, bin_size):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        handle.write("# Histogram data for Circos\n")
        handle.write("# chromosome\tstart\tend\tvalue\n")
        for chrom, length in chrom_lengths.items():
            for idx, value in enumerate(bins[chrom]):
                start = idx * bin_size + 1
                end = min((idx + 1) * bin_size, length)
                handle.write(f"{chrom}\t{start}\t{end}\t{value}\n")


def main():
    args = parse_args()
    strains = [args.strain1_name, args.strain2_name]
    chroms = load_karyotype(args.karyotype, strains)
    output_dir = Path(args.output_dir)

    gtf_by_sample = {
        args.strain1_name: args.strain1_gtf,
        args.strain2_name: args.strain2_gtf,
    }
    introns_by_sample = {}
    for sample, gtf in gtf_by_sample.items():
        if sample != args.non_introner_reference:
            introns_by_sample[sample] = extract_introns_from_gtf(gtf)

    for sample in strains:
        if not chroms[sample]:
            raise SystemExit(f"No chromosomes found for {sample} in {args.karyotype}")

        introner_bins = make_empty_bins(chroms[sample], args.bin_size)
        introner_stats = count_introner_matrix(
            args.introner_matrix, sample, introner_bins, args.bin_size
        )
        write_histogram(
            output_dir / f"{sample}.introner_histogram.txt",
            chroms[sample],
            introner_bins,
            args.bin_size,
        )

        non_introner_bins = make_empty_bins(chroms[sample], args.bin_size)
        non_introner_stats = count_non_introner_matrix(
            args.non_introner_matrix,
            sample,
            args.non_introner_reference,
            introns_by_sample.get(sample, {}),
            non_introner_bins,
            args.bin_size,
        )
        write_histogram(
            output_dir / f"{sample}.intron_histogram.txt",
            chroms[sample],
            non_introner_bins,
            args.bin_size,
        )

        print(
            f"{sample}: introners counted={introner_stats['counted']} "
            f"skipped={introner_stats['skipped']}; "
            f"non-introner introns counted={non_introner_stats['counted']} "
            f"skipped={non_introner_stats['skipped']}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()

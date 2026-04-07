#!/usr/bin/env python3
"""
Summarize SNP counts by site class (non-synonymous, synonymous, 4-fold degenerate, non-CDS).

Parses the degenotate BED file to classify CDS positions by degeneracy,
then scans the all-SNPs VCF to count SNPs in each class.
Excludes scaffold_2 (mating type region).

Reports stats for all samples and optionally for a Group 1 (intronerful) subset.
"""

import argparse
import gzip
import re
import sys


GROUP1_SAMPLES = [
    "CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482",
    "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833",
]


def parse_contig_lengths(vcf_file, exclude_scaffold):
    """Extract contig lengths from VCF header, excluding specified scaffold."""
    total = 0
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                break
            m = re.match(r'##contig=<ID=(.+),length=(\d+)>', line)
            if m:
                contig_id = m.group(1)
                length = int(m.group(2))
                if exclude_scaffold not in contig_id:
                    total += length
    return total


def parse_degenotate(bed_file, exclude_scaffold):
    """Parse degenotate BED to build position -> degeneracy lookup.

    Returns dict of (scaffold, 1-based pos) -> degeneracy (0, 2, 3, or 4).
    Positions with '.' degeneracy are skipped (non-coding in BED).
    """
    sites = {}
    opener = gzip.open if bed_file.endswith('.gz') else open
    with opener(bed_file, 'rt') as f:
        for line in f:
            cols = line.split('\t')
            scaffold = cols[0]
            if exclude_scaffold in scaffold:
                continue
            deg = cols[4].strip()
            if deg == '.':
                continue
            pos = int(cols[1]) + 1
            sites[(scaffold, pos)] = int(deg)
    return sites


def is_polymorphic(gt_fields, sample_indices):
    """Check if a site is polymorphic among the given sample indices.

    Parses GT from each sample field. Returns True if more than one
    distinct non-missing allele is observed.
    """
    alleles = set()
    for idx in sample_indices:
        gt = gt_fields[idx].split(':')[0]
        if gt == '.':
            continue
        alleles.add(gt)
    return len(alleles) > 1


def write_summary(out, group_label, site_counts, snp_counts, genome_size):
    """Write summary rows for one sample group."""
    n_0fold, n_syn, n_4fold, n_noncds = site_counts
    s_0fold, s_syn, s_4fold, s_noncds, s_total = snp_counts
    for label, n_sites, n_snps in [
        ("non-synonymous (0-fold)", n_0fold, s_0fold),
        ("synonymous (2+3+4-fold)", n_syn, s_syn),
        ("4-fold degenerate", n_4fold, s_4fold),
        ("non-CDS", n_noncds, s_noncds),
        ("total", genome_size, s_total),
    ]:
        density = n_snps / n_sites if n_sites > 0 else 0
        out.write(f"{group_label}\t{label}\t{n_sites}\t{n_snps}\t{density:.6f}\n")


def main():
    parser = argparse.ArgumentParser(description="Summarize SNP counts by site class.")
    parser.add_argument("--bed", required=True, help="Degenotate degeneracy-all-sites.bed")
    parser.add_argument("--vcf", required=True, help="All-SNPs VCF (mpusilla.snps.with_ref.vcf.gz)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--exclude-scaffold", default="scaffold_2",
                        help="Scaffold to exclude (default: scaffold_2)")
    args = parser.parse_args()

    # Get genome size
    print("Parsing contig lengths from VCF header...", file=sys.stderr)
    genome_size = parse_contig_lengths(args.vcf, args.exclude_scaffold)
    print(f"  Genome size (excl {args.exclude_scaffold}): {genome_size:,}", file=sys.stderr)

    # Parse degeneracy
    print("Parsing degenotate BED file...", file=sys.stderr)
    deg_sites = parse_degenotate(args.bed, args.exclude_scaffold)

    # Count sites by class
    n_0fold = sum(1 for d in deg_sites.values() if d == 0)
    n_2fold = sum(1 for d in deg_sites.values() if d == 2)
    n_3fold = sum(1 for d in deg_sites.values() if d == 3)
    n_4fold = sum(1 for d in deg_sites.values() if d == 4)
    n_syn = n_2fold + n_3fold + n_4fold
    n_cds = len(deg_sites)
    n_noncds = genome_size - n_cds

    print(f"  0-fold (non-synonymous): {n_0fold:,}", file=sys.stderr)
    print(f"  Synonymous (2+3+4-fold): {n_syn:,}", file=sys.stderr)
    print(f"  4-fold degenerate: {n_4fold:,}", file=sys.stderr)
    print(f"  Non-CDS: {n_noncds:,}", file=sys.stderr)

    # Scan VCF — identify sample column indices
    print("Scanning VCF for SNPs...", file=sys.stderr)

    # Counters: [all_samples, group1]
    counts = {
        '0fold': [0, 0], 'syn': [0, 0], '4fold': [0, 0],
        'noncds': [0, 0], 'total': [0, 0],
    }
    group1_indices = None

    opener = gzip.open if args.vcf.endswith('.gz') else open
    with opener(args.vcf, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # Parse header to find sample indices
                header = line.strip().split('\t')
                samples = header[9:]
                all_indices = list(range(len(samples)))
                group1_indices = [i for i, s in enumerate(samples) if s in GROUP1_SAMPLES]
                print(f"  All samples ({len(all_indices)}): {', '.join(samples)}", file=sys.stderr)
                print(f"  Group 1 samples ({len(group1_indices)}): "
                      f"{', '.join(samples[i] for i in group1_indices)}", file=sys.stderr)
                continue

            cols = line.strip().split('\t')
            scaffold = cols[0]
            if args.exclude_scaffold in scaffold:
                continue
            pos = int(cols[1])
            gt_fields = cols[9:]

            # Check polymorphism for each group
            all_poly = is_polymorphic(gt_fields, all_indices)
            g1_poly = is_polymorphic(gt_fields, group1_indices)

            # Classify by site class
            key = (scaffold, pos)
            if key in deg_sites:
                deg = deg_sites[key]
                if deg == 0:
                    cat = '0fold'
                else:
                    cat = 'syn'
            else:
                cat = 'noncds'

            if all_poly:
                counts[cat][0] += 1
                counts['total'][0] += 1
                if cat == 'syn' and deg == 4:
                    counts['4fold'][0] += 1
            if g1_poly:
                counts[cat][1] += 1
                counts['total'][1] += 1
                if cat == 'syn' and deg == 4:
                    counts['4fold'][1] += 1

    for label, idx in [("All samples", 0), ("Group 1", 1)]:
        print(f"  {label} — Total SNPs: {counts['total'][idx]:,}", file=sys.stderr)

    # Write output
    site_counts = (n_0fold, n_syn, n_4fold, n_noncds)
    with open(args.output, 'w') as out:
        out.write("sample_group\tsite_class\ttotal_sites\tsnp_count\tsnp_density\n")
        for group_idx, group_label in [(0, "all_samples"), (1, "group1")]:
            snp_counts = (
                counts['0fold'][group_idx],
                counts['syn'][group_idx],
                counts['4fold'][group_idx],
                counts['noncds'][group_idx],
                counts['total'][group_idx],
            )
            write_summary(out, group_label, site_counts, snp_counts, genome_size)

    print(f"Output written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

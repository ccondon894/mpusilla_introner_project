#!/usr/bin/env python3
"""
Compute codon-level insertion site fingerprints for introner loci.

For each introner with a gene annotation, determines the exact CDS-relative
insertion position using splice site annotations and/or genome sequence analysis.
Produces a (gene_id, codon_number, codon_offset) fingerprint that can be compared
across samples to distinguish ancestral shared introners from independent insertions.

Coordinate conventions:
  - BED coordinates include 100bp flanking regions on each side
  - True introner body = [BED_start + 100, BED_end - 100)
  - splice_info GT@N = GT dinucleotide at position N within the +strand body sequence
  - GTF coordinates are 1-based, closed intervals

For forward-oriented introners (+ strand genes):
  - GT@N on +strand is the functional 5' splice donor
  - Insertion position = BED_start + 100 + N (0-based genomic)

For reverse-oriented introners (- strand genes):
  - GT@N on +strand is NOT the functional donor
  - Script reverse-complements the body and searches for GT on the - strand
  - Insertion position computed from the body's high-coordinate end
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

import pysam
from Bio.Seq import Seq


FLANK_LENGTH = 100
SPLICE_SEARCH_WINDOW = 20

# Insertions within this many codons of an exon end (or start) are reclassified
# as the adjacent intron, since they sit at the splice site boundary. This
# handles cases where miniprot annotated different exon lengths in different
# samples but the introner is at the same biological intron-exon junction.
BOUNDARY_CODON_TOLERANCE = 3

# Number of amino acids on each side of the insertion site to include in
# the flanking context fingerprint. We extract a larger window (5 on each
# side = 10 aa total) so the comparison can handle small shifts in the
# computed insertion position (±1-2 codons from splice site detection
# variability). The comparison uses a sliding 6-mer match within this
# larger window.
FLANKING_AA_N = 5


def parse_gtf_cds(gtf_path):
    """Parse CDS features from GTF, grouped by gene_id.

    Returns dict: gene_id -> list of CDS dicts with contig, start, end, strand, frame.
    """
    cds_by_gene = defaultdict(list)

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue

            contig = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])    # 1-based, inclusive
            strand = fields[6]
            frame = int(fields[7]) if fields[7] != '.' else 0

            gene_id = None
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break

            if gene_id:
                cds_by_gene[gene_id].append({
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'frame': frame
                })

    return cds_by_gene


class GeneContext:
    """Precomputed CDS and protein data for a single gene.

    Builds the concatenated CDS sequence (across all exons in transcript
    order) and the translated amino acid sequence. Provides helpers for
    extracting the flanking amino acid context around an insertion site,
    which is the biologically-meaningful locus identifier that's invariant
    to miniprot exon boundary annotation differences.
    """

    def __init__(self, cds_exons, genome):
        """Build the CDS and protein sequence from the given exon list.

        Args:
            cds_exons: list of dicts with contig, start, end, strand (all on
                the same contig and strand)
            genome: open pysam.FastaFile
        """
        self.contig = cds_exons[0]['contig']
        self.strand = cds_exons[0]['strand']

        if self.strand == '+':
            self.sorted_exons = sorted(cds_exons, key=lambda x: x['start'])
        else:
            self.sorted_exons = sorted(cds_exons, key=lambda x: x['start'],
                                        reverse=True)

        # Build the concatenated CDS sequence in transcript direction.
        # Track cumulative CDS length before each exon so we can locate
        # intron splice sites in the CDS/protein coordinate system.
        cds_parts = []
        self.exon_cum_lengths = [0]  # cum[i] = cum CDS length before exon i
        for exon in self.sorted_exons:
            try:
                # GTF is 1-based inclusive; pysam.fetch is 0-based half-open
                exon_seq = genome.fetch(
                    self.contig, exon['start'] - 1, exon['end']).upper()
            except (ValueError, KeyError):
                self.cds_seq = None
                self.aa_seq = None
                return
            if self.strand == '-':
                exon_seq = reverse_complement(exon_seq)
            cds_parts.append(exon_seq)
            self.exon_cum_lengths.append(
                self.exon_cum_lengths[-1] + len(exon_seq))

        self.cds_seq = ''.join(cds_parts)

        # Translate the CDS. Trim any trailing partial codon since the CDS
        # may not be a multiple of 3 (happens with miniprot annotations that
        # don't align perfectly to the reference protein).
        trimmed = self.cds_seq[:len(self.cds_seq) - (len(self.cds_seq) % 3)]
        try:
            self.aa_seq = str(Seq(trimmed).translate())
        except Exception:
            self.aa_seq = None

    def get_flanking_context(self, location_type, cds_position=None,
                              intron_number=None, n_flank=FLANKING_AA_N):
        """Extract n_flank amino acids on each side of the insertion site.

        Returns a 2*n_flank amino acid string representing the flanking
        context of the insertion. The comparison function (see
        classify_sharing_status.py) uses a sliding-window match on this
        context, so extracting extra flank aa gives robustness to small
        shifts in the computed insertion position (±1-2 codons) caused
        by splice site detection variability across samples.

        For in-CDS insertions (location_type='cds'):
            codon_idx = cds_position // 3 (the codon containing the last
            exonic nucleotide before the insertion)

        For in-intron insertions (location_type='intron'):
            codon_idx = (cum_cds_length_at_end_of_upstream_exon - 1) // 3
            (the codon at the splice site; may be a split codon for phase
            1/2 introns, but it's still translatable from the full CDS)

        Insufficient context at the gene edges is padded with '-'.
        """
        if self.aa_seq is None or not self.aa_seq:
            return None

        if location_type == 'cds':
            if cds_position is None:
                return None
            codon_idx = cds_position // 3
        elif location_type == 'intron':
            if intron_number is None:
                return None
            if intron_number < 1 or intron_number >= len(self.exon_cum_lengths):
                return None
            cum = self.exon_cum_lengths[intron_number]
            if cum == 0:
                return None
            codon_idx = (cum - 1) // 3
        else:
            return None

        # Extract a window of 2*n_flank amino acids centered on the
        # "insertion boundary" (between codon_idx and codon_idx+1).
        # - Upstream n_flank codons end at codon_idx inclusive
        # - Downstream n_flank codons start at codon_idx+1
        up_start = codon_idx - n_flank + 1
        up_end = codon_idx + 1
        dn_start = codon_idx + 1
        dn_end = codon_idx + 1 + n_flank

        # Pad left if needed
        left_pad = '-' * max(0, -up_start)
        real_up_start = max(0, up_start)

        # Pad right if needed
        right_pad = '-' * max(0, dn_end - len(self.aa_seq))
        real_dn_end = min(len(self.aa_seq), dn_end)

        up = self.aa_seq[real_up_start:up_end] if real_up_start < up_end else ''
        dn = self.aa_seq[dn_start:real_dn_end] if dn_start < real_dn_end else ''

        ctx = left_pad + up + dn + right_pad
        # Ensure exact length
        if len(ctx) != 2 * n_flank:
            # Shouldn't happen but handle defensively
            ctx = ctx[:2 * n_flank].ljust(2 * n_flank, '-')
        return ctx


def parse_splice_offset(splice_str):
    """Extract the donor splice site offset from a splice_info string.

    Parses GT@N or GC@N from strings like 'GT@8' or 'GT@3;AG@175'.
    Returns (splice_type, offset) or (None, None).
    """
    if not splice_str or splice_str in ('NA', ''):
        return None, None

    for part in splice_str.split(';'):
        part = part.strip()
        if '@' in part:
            stype, offset = part.split('@', 1)
            if stype in ('GT', 'GC'):
                return stype, int(offset)

    return None, None


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq[::-1].translate(table)


def find_donor_on_minus_strand(body_seq, search_window=SPLICE_SEARCH_WINDOW):
    """Find the 5' splice donor (GT or GC) on the minus strand.

    Reverse-complements the body and searches for GT/GC near the 5' end
    of the minus-strand sequence (= near the 3' end of the + strand body).

    Returns (splice_type, offset_in_revcomp) or (None, None).
    """
    body_rc = reverse_complement(body_seq)
    window = body_rc[:search_window]

    gt = re.search(r'GT', window)
    if gt:
        return 'GT', gt.start()

    gc = re.search(r'GC', window)
    if gc:
        return 'GC', gc.start()

    return None, None


def compute_cds_position(insertion_pos_1based, cds_exons, strand,
                          boundary_codon_tolerance=BOUNDARY_CODON_TOLERANCE):
    """Locate an insertion site within a gene's CDS structure.

    For + strand: CDS exons ordered by ascending genomic position (5' to 3').
    For - strand: CDS exons ordered by descending genomic position (5' to 3').

    Returns a tuple (location_type, cds_position, exon_or_intron_number):
      - ('cds', cds_position, exon_number): inside a CDS exon
      - ('intron', None, intron_number): in a gap between CDS exons OR within
        boundary_codon_tolerance codons of an exon edge. Insertions at the
        last few codons of exon N (or first few codons of exon N+1) are
        treated as "intron N" because they sit at the splice junction. This
        handles miniprot exon length annotation differences for boundary
        insertions.
      - (None, None, None): outside the entire CDS

    Both exon_number and intron_number are 1-indexed in transcript order.
    """
    if strand == '+':
        sorted_exons = sorted(cds_exons, key=lambda x: x['start'])
        cds_offset = 0

        for i, exon in enumerate(sorted_exons):
            exon_len = exon['end'] - exon['start'] + 1

            if insertion_pos_1based < exon['start']:
                if i == 0:
                    return None, None, None
                return 'intron', None, i

            if insertion_pos_1based <= exon['end']:
                # Inside exon i+1 (1-indexed in transcript order)
                # Compute position within exon and check if it's near a boundary
                pos_in_exon = insertion_pos_1based - exon['start']
                codon_in_exon = pos_in_exon // 3
                exon_codon_count = exon_len // 3
                boundary_nt = boundary_codon_tolerance * 3

                # Within boundary_nt of the END of this exon → next intron
                # (only if there IS a next exon)
                if (i + 1 < len(sorted_exons) and
                        (exon['end'] - insertion_pos_1based) < boundary_nt):
                    return 'intron', None, i + 1

                # Within boundary_nt of the START of this exon → previous intron
                # (only if there IS a previous exon)
                if i > 0 and (insertion_pos_1based - exon['start']) < boundary_nt:
                    return 'intron', None, i

                return 'cds', cds_offset + pos_in_exon, i + 1

            cds_offset += exon_len

        return None, None, None

    else:  # '-' strand
        sorted_exons = sorted(cds_exons, key=lambda x: x['start'], reverse=True)
        cds_offset = 0

        for i, exon in enumerate(sorted_exons):
            exon_len = exon['end'] - exon['start'] + 1

            if insertion_pos_1based > exon['end']:
                if i == 0:
                    return None, None, None
                return 'intron', None, i

            if insertion_pos_1based >= exon['start']:
                # Inside exon i+1 in transcript order (- strand: high coords first)
                pos_in_exon = exon['end'] - insertion_pos_1based
                boundary_nt = boundary_codon_tolerance * 3

                # On - strand, the END of an exon (in transcript order) is at
                # the LOW genomic coordinate. So "near the end" means
                # insertion_pos is close to exon['start'].
                if (i + 1 < len(sorted_exons) and
                        (insertion_pos_1based - exon['start']) < boundary_nt):
                    return 'intron', None, i + 1

                # Near the START of the exon (in transcript order) = near the
                # HIGH genomic coordinate.
                if i > 0 and (exon['end'] - insertion_pos_1based) < boundary_nt:
                    return 'intron', None, i

                return 'cds', cds_offset + pos_in_exon, i + 1

            cds_offset += exon_len

        return None, None, None


def main():
    parser = argparse.ArgumentParser(
        description='Compute codon-level insertion fingerprints for introner loci')
    parser.add_argument('--bed', required=True,
                        help='Candidate loci BED file (coords include 100bp flanks)')
    parser.add_argument('--gtf', required=True,
                        help='GTF annotation file with CDS features')
    parser.add_argument('--genome', required=True,
                        help='Indexed genome FASTA (.fa with .fai)')
    parser.add_argument('--sample', required=True,
                        help='Sample name')
    parser.add_argument('--output', required=True,
                        help='Output TSV file')
    args = parser.parse_args()

    # Load data
    cds_by_gene = parse_gtf_cds(args.gtf)
    print(f"Loaded CDS features for {len(cds_by_gene)} genes from GTF")

    genome = pysam.FastaFile(args.genome)

    # Cache of GeneContext objects, keyed by (gene_id, contig). Building
    # these is expensive (requires genome fetch + translation), and many
    # introners share the same gene.
    gene_contexts = {}

    # Track statistics
    stats = defaultdict(int)

    results = []

    with open(args.bed) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            stats['total'] += 1

            contig = row['#chrom']
            bed_start = int(row['start'])  # 0-based, includes flanks
            bed_end = int(row['end'])
            introner_id = row['name']
            family = row['family']
            gene_info = row['gene_info']
            splice_info = row['splice_info']
            orientation = row['orientation']

            # Base result dict (filled in regardless of outcome)
            result = {
                'sample': args.sample,
                'introner_id': introner_id,
                'contig': contig,
                'start': bed_start,
                'end': bed_end,
                'family': family,
                'gene_id': '',
                'strand': '',
                'insertion_pos': '',
                'location_type': '',
                'cds_position': '',
                'codon_number': '',
                'codon_offset': '',
                'exon_number': '',
                'intron_number': '',
                'flanking_aa_context': '',
                'confidence': '',
            }

            # --- Check gene annotation ---
            if not gene_info or gene_info in ('NA', ''):
                stats['no_gene'] += 1
                result['confidence'] = 'no_gene'
                results.append(result)
                continue

            result['gene_id'] = gene_info

            # --- Find CDS features for this gene on this contig ---
            if gene_info not in cds_by_gene:
                stats['gene_not_in_gtf'] += 1
                result['confidence'] = 'gene_not_in_gtf'
                results.append(result)
                continue

            cds_on_contig = [c for c in cds_by_gene[gene_info]
                             if c['contig'] == contig]
            if not cds_on_contig:
                stats['no_cds_on_contig'] += 1
                result['confidence'] = 'no_cds_on_contig'
                results.append(result)
                continue

            strand = cds_on_contig[0]['strand']
            result['strand'] = strand

            # --- Compute insertion position ---
            # True body coordinates (excluding flanks)
            body_start = bed_start + FLANK_LENGTH   # 0-based
            body_end = bed_end - FLANK_LENGTH        # 0-based exclusive
            body_len = body_end - body_start

            if body_len <= 0:
                stats['invalid_body'] += 1
                result['confidence'] = 'invalid_body'
                results.append(result)
                continue

            if strand == '+':
                # Forward: use GT@N from splice_info directly
                splice_type, splice_offset = parse_splice_offset(splice_info)
                if splice_type is None:
                    stats['no_splice_site'] += 1
                    result['confidence'] = 'no_splice_site'
                    results.append(result)
                    continue

                # Last exonic nt before the intron (0-based genomic)
                last_exonic_0 = body_start + splice_offset - 1
                # 1-based for GTF comparison
                insertion_pos_1 = last_exonic_0 + 1

            else:  # strand == '-'
                # Reverse: find GT on the - strand by reverse-complementing the body
                try:
                    body_seq = genome.fetch(contig, body_start, body_end)
                except (ValueError, KeyError):
                    stats['genome_fetch_failed'] += 1
                    result['confidence'] = 'genome_fetch_failed'
                    results.append(result)
                    continue

                rc_type, rc_offset = find_donor_on_minus_strand(body_seq)
                if rc_type is None:
                    stats['no_splice_site'] += 1
                    result['confidence'] = 'no_splice_site'
                    results.append(result)
                    continue

                # Last exonic nt on - strand (0-based genomic)
                # body_rc[P] corresponds to +strand pos (body_end - 1 - P)
                # The GT donor on -strand starts at body_rc[P]
                # Last exonic nt is one position to the RIGHT in genomic coords
                last_exonic_0 = body_end - rc_offset
                insertion_pos_1 = last_exonic_0 + 1  # 1-based

            result['insertion_pos'] = last_exonic_0

            # --- Map to CDS coordinates ---
            location_type, cds_pos, exon_or_intron = compute_cds_position(
                insertion_pos_1, cds_on_contig, strand)

            if location_type is None:
                stats['not_in_cds'] += 1
                result['confidence'] = 'not_in_cds'
                results.append(result)
                continue

            result['location_type'] = location_type

            if location_type == 'cds':
                result['cds_position'] = cds_pos
                result['codon_number'] = cds_pos // 3
                result['codon_offset'] = cds_pos % 3
                result['exon_number'] = exon_or_intron
                stats['cds_success'] += 1
            else:  # 'intron'
                result['intron_number'] = exon_or_intron
                stats['intron_success'] += 1

            # --- Extract flanking amino acid context ---
            gc_key = (gene_info, contig)
            if gc_key not in gene_contexts:
                gene_contexts[gc_key] = GeneContext(cds_on_contig, genome)
            gene_context = gene_contexts[gc_key]

            if location_type == 'cds':
                aa_context = gene_context.get_flanking_context(
                    'cds', cds_position=cds_pos)
            else:  # intron
                aa_context = gene_context.get_flanking_context(
                    'intron', intron_number=exon_or_intron)

            if aa_context:
                result['flanking_aa_context'] = aa_context
                stats['aa_context_success'] += 1
            else:
                stats['aa_context_failed'] += 1

            result['confidence'] = 'high'
            stats['success'] += 1
            results.append(result)

    genome.close()

    # Print summary
    print(f"\nInsertion fingerprint results for {args.sample}:")
    print(f"  Total loci:         {stats['total']}")
    print(f"  Successful:         {stats['success']} "
          f"({100 * stats['success'] / max(1, stats['total']):.1f}%)")
    print(f"    in CDS exon:      {stats['cds_success']}")
    print(f"    in intron:        {stats['intron_success']}")
    print(f"    with aa context:  {stats['aa_context_success']}")
    if stats['aa_context_failed']:
        print(f"    aa context failed:{stats['aa_context_failed']}")
    print(f"  No gene annotation: {stats['no_gene']}")
    print(f"  Gene not in GTF:    {stats['gene_not_in_gtf']}")
    print(f"  No CDS on contig:   {stats['no_cds_on_contig']}")
    print(f"  No splice site:     {stats['no_splice_site']}")
    print(f"  Not in CDS:         {stats['not_in_cds']}")
    if stats['invalid_body']:
        print(f"  Invalid body:       {stats['invalid_body']}")
    if stats['genome_fetch_failed']:
        print(f"  Genome fetch failed:{stats['genome_fetch_failed']}")

    # Write output
    fieldnames = ['sample', 'introner_id', 'contig', 'start', 'end', 'family',
                  'gene_id', 'strand', 'insertion_pos', 'location_type',
                  'cds_position', 'codon_number', 'codon_offset',
                  'exon_number', 'intron_number', 'flanking_aa_context',
                  'confidence']

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote {len(results)} records to {args.output}")


if __name__ == '__main__':
    main()

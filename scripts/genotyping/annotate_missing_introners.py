#!/usr/bin/env python3
"""
Annotate missing gene and family data in genotype matrix.

Strategy:
1. Gene annotation: Use bedtools-style overlap for ALL introners (scenarios 1, 2, 3)
   - Load gene.bed files for each sample
   - Find overlapping genes for each introner coordinate

2. Family annotation: Use consensus sequence matching for ONLY scenario 1 (present) introners
   - Build consensus sequences from reference families using MAFFT alignment + PWM
   - Compare against reference introner family consensus sequences using pairwise alignment
   - Assign family with highest similarity that meets 80-80 rule (similarity ≥ 0.7 AND coverage ≥ 0.8)
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
import tempfile
import re
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import numpy as np
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_gene_bed_files(processed_ann_dir, samples):
    """Load all gene BED files into a dictionary by sample."""
    gene_intervals = defaultdict(lambda: [])

    for sample in samples:
        bed_file = os.path.join(processed_ann_dir, f"{sample}.gene.bed")
        if not os.path.exists(bed_file):
            logger.warning(f"Gene BED file not found for {sample}: {bed_file}")
            continue

        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue

                contig = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                gene_name = fields[3]
                strand = fields[5]

                gene_intervals[sample].append({
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'gene': gene_name,
                    'strand': strand
                })

    logger.info(f"Loaded gene annotations for {len(gene_intervals)} samples")
    return gene_intervals


def find_overlapping_gene(sample, contig, start, end, gene_intervals):
    """Find a gene that overlaps the given introner coordinates."""
    if sample not in gene_intervals:
        return None

    for gene in gene_intervals[sample]:
        # Check if gene's contig matches
        if gene['contig'] != contig:
            continue

        # Check for overlap: gene.start < introner.end AND gene.end > introner.start
        if gene['start'] < end and gene['end'] > start:
            return gene['gene']

    return None


def load_fasta_files(output_dir, samples):
    """Load introner FASTA files for all samples."""
    fasta_dict = defaultdict(dict)

    for sample in samples:
        fasta_file = os.path.join(output_dir, f"{sample}.candidate_loci_plus_flanks.filtered.similarity_checked.fa")
        if not os.path.exists(fasta_file):
            logger.warning(f"FASTA file not found for {sample}: {fasta_file}")
            continue

        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id
                fasta_dict[sample][seq_id] = str(record.seq)

            logger.info(f"Loaded {len(fasta_dict[sample])} sequences for {sample}")
        except Exception as e:
            logger.error(f"Error loading FASTA for {sample}: {e}")

    return fasta_dict


def load_genome_assemblies(genome_dir, samples):
    """Load reference genome FASTA files for all samples (indexed for fast lookup)."""
    genomes = {}

    for sample in samples:
        genome_file = os.path.join(genome_dir, f"{sample}.vg_paths.fa")
        if not os.path.exists(genome_file):
            logger.warning(f"Genome file not found for {sample}: {genome_file}")
            continue

        try:
            # Create a dictionary indexed by contig name for fast lookup
            genome_dict = {}
            for record in SeqIO.parse(genome_file, 'fasta'):
                contig = record.id
                genome_dict[contig] = str(record.seq)

            genomes[sample] = genome_dict
            logger.info(f"Loaded {len(genome_dict)} contigs for {sample}")
        except Exception as e:
            logger.error(f"Error loading genome for {sample}: {e}")

    return genomes


def get_family_from_ortholog_group(ortholog_id, df):
    """Get family assignment from other samples in the same ortholog group.

    Returns the family if found in any other row with the same ortholog_id,
    otherwise returns None.
    """
    group_rows = df[df['ortholog_id'] == ortholog_id]

    # Look for any row with a family assignment in this group
    for idx, row in group_rows.iterrows():
        if pd.notna(row['family']) and row['family'] != '':
            return row['family']

    return None


def extract_sequence_from_genome(sample, contig, start, end, genomes, flanking_length=100):
    """Extract sequence from genome assembly, including flanking regions.

    Returns: (full_seq_with_flanks, success)
    Adjusts coordinates to include flanking regions on both sides.
    """
    if sample not in genomes or contig not in genomes[sample]:
        return None, False

    genome_seq = genomes[sample][contig]
    genome_len = len(genome_seq)

    # Adjust coordinates to include flanking regions
    flank_start = max(0, start - flanking_length)
    flank_end = min(genome_len, end + flanking_length)

    # Extract sequence
    full_seq = genome_seq[flank_start:flank_end]

    # Verify we got a reasonable sequence
    if len(full_seq) < 2 * flanking_length:
        return None, False

    return full_seq, True


# ==================== Consensus Building Functions (from validate_introners.py) ====================

def extract_introner_sequences(ref_dir, family_file):
    """Extract introner sequences (lowercase parts) from reference file"""
    introner_seqs = []

    file_path = os.path.join(ref_dir, family_file)
    for record in SeqIO.parse(file_path, "fasta"):
        seq_str = str(record.seq)

        # Find lowercase section (introner sequence)
        introner_match = re.search(r'[a-z]+', seq_str)
        if introner_match:
            start, end = introner_match.span()
            introner_seq = seq_str[start:end]

            # Create a new SeqRecord for the introner sequence
            introner_rec = SeqRecord(
                Seq(introner_seq),
                id=record.id,
                description="Introner sequence"
            )

            introner_seqs.append(introner_rec)

    return introner_seqs


def run_mafft_alignment(sequences, mafft_path='mafft'):
    """Run MAFFT alignment on a set of sequences"""
    # Create temporary files for input and output
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_in:
        # Write sequences to temporary input file
        SeqIO.write(sequences, temp_in, "fasta")
        temp_in_name = temp_in.name

    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_out:
        temp_out_name = temp_out.name

    # Construct MAFFT command
    cmd = [mafft_path, "--auto", temp_in_name]

    # Run MAFFT
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Write MAFFT output to temporary output file
        with open(temp_out_name, 'w') as f:
            f.write(result.stdout)

        # Parse alignment
        alignment = AlignIO.read(temp_out_name, "fasta")

    except subprocess.CalledProcessError as e:
        logger.error(f"MAFFT alignment failed: {e}")
        logger.error(f"STDERR: {e.stderr}")
        alignment = None
    finally:
        # Clean up temporary files
        os.unlink(temp_in_name)
        os.unlink(temp_out_name)

    return alignment


def create_pwm(alignment):
    """Create position weight matrix from an alignment"""
    # Dictionary to map nucleotides to indices
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3, '-': 4}

    # Initialize the PWM (rows: A, C, G, T, gap; columns: positions)
    pwm = np.zeros((5, alignment.get_alignment_length()))

    # Fill the PWM
    for i in range(alignment.get_alignment_length()):
        for record in alignment:
            nuc = record.seq[i]
            if nuc in nuc_to_idx:
                pwm[nuc_to_idx[nuc], i] += 1

    # Normalize the PWM by column
    col_sums = pwm.sum(axis=0)
    pwm = pwm / col_sums[np.newaxis, :]

    return pwm


def generate_consensus(pwm):
    """Generate consensus sequence from a PWM"""
    # Define mapping from PWM index to nucleotide
    idx_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-'}

    # Generate consensus sequence
    consensus = []
    for i in range(pwm.shape[1]):
        # Get the nucleotide with the highest probability
        max_idx = np.argmax(pwm[:, i])

        # Only add the nucleotide if it's not a gap
        if max_idx != 4 and pwm[max_idx, i] >= 0.3:  # Threshold to avoid ambiguous positions
            consensus.append(idx_to_nuc[max_idx])

    return ''.join(consensus)


def build_consensus_sequences(ref_dir, mafft_path='mafft'):
    """Build consensus sequences for each reference family using MAFFT and PWM"""
    logger.info("Building consensus sequences for reference families...")

    consensus_sequences = {}

    # Find all family FASTA files matching GCA_000151265.1_fam pattern
    fam_files = [f for f in os.listdir(ref_dir) if f.startswith("GCA_000151265.1_fam")]

    if not fam_files:
        logger.warning(f"No family FASTA files found in {ref_dir}")
        return consensus_sequences

    for fam_file in fam_files:
        # Extract family ID from filename
        family_match = re.search(r'fam(\d+)', fam_file)
        if family_match:
            family_id = family_match.group(1)
        else:
            logger.warning(f"Could not extract family ID from {fam_file}, skipping")
            continue

        # Extract introner sequences
        introner_seqs = extract_introner_sequences(ref_dir, fam_file)

        if len(introner_seqs) < 2:
            logger.warning(f"Family {family_id} has only {len(introner_seqs)} sequences, skipping alignment")
            if introner_seqs:
                consensus_sequences[family_id] = str(introner_seqs[0].seq)
            continue

        # Align sequences using MAFFT
        alignment = run_mafft_alignment(introner_seqs, mafft_path)

        if alignment:
            # Create PWM from alignment
            pwm = create_pwm(alignment)

            # Generate consensus sequence
            consensus = generate_consensus(pwm)

            consensus_sequences[family_id] = consensus

            logger.info(f"Generated consensus for family {family_id} (length: {len(consensus)})")
        else:
            logger.error(f"Failed to align sequences for family {family_id}")

    return consensus_sequences


# ==================== Sequence Similarity Functions (from validate_introners.py) ====================

def process_alignment(alignment, seq1, seq2, similarity_cutoff=0.7, coverage_cutoff=0.8):
    """Process alignment to extract similarity and coverage metrics"""
    # Extract aligned sequences
    aligned_seq1, aligned_seq2 = alignment[0], alignment[1]

    # Count matches, mismatches and gaps
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2)
                 if a == b and a != '-' and b != '-')

    # Calculate similarity as percentage of matches relative to aligned length
    aligned_length = len(aligned_seq1)
    non_gap_positions = sum(1 for a, b in zip(aligned_seq1, aligned_seq2)
                           if a != '-' or b != '-')

    similarity = (matches / non_gap_positions) if non_gap_positions > 0 else 0.0

    # Calculate coverage of reference sequence (seq2)
    seq2_aligned_length = sum(1 for b in aligned_seq2 if b != '-')
    coverage = seq2_aligned_length / len(seq2)

    # Check if alignment meets the thresholds
    meets_rule = similarity >= similarity_cutoff and coverage >= coverage_cutoff

    return {
        'similarity': similarity,
        'coverage': coverage,
        'meets_rule': meets_rule
    }


def calculate_similarity_80_80_rule(seq1, seq2, similarity_cutoff=0.7, coverage_cutoff=0.8):
    """
    Calculate similarity between two sequences using pairwise alignment:
    - At least X% sequence identity (default: 70%)
    - Over at least Y% of the reference sequence length (default: 80%)

    Returns:
    - similarity: float - the sequence identity score
    - coverage: float - the proportion of reference sequence covered
    - meets_rule: bool - whether alignment satisfies the given thresholds
    - orientation: str - "forward" or "reverse"
    """
    # Convert to strings if they aren't already
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()

    # Create reverse complement of seq1
    seq1_rc = str(Seq(seq1).reverse_complement())

    # Calculate similarity with original sequence
    forward_alignments = pairwise2.align.globalms(seq1, seq2,
                                                 2, -1,  # match/mismatch
                                                 -2, -0.5,  # gap open/extend
                                                 one_alignment_only=True)

    # Calculate similarity with reverse complement
    reverse_alignments = pairwise2.align.globalms(seq1_rc, seq2,
                                                 2, -1,  # match/mismatch
                                                 -2, -0.5,  # gap open/extend
                                                 one_alignment_only=True)

    # Process forward alignment
    if forward_alignments:
        forward_result = process_alignment(forward_alignments[0], seq1, seq2, similarity_cutoff, coverage_cutoff)
    else:
        forward_result = {'similarity': 0.0, 'coverage': 0.0, 'meets_rule': False}

    # Process reverse alignment
    if reverse_alignments:
        reverse_result = process_alignment(reverse_alignments[0], seq1_rc, seq2, similarity_cutoff, coverage_cutoff)
    else:
        reverse_result = {'similarity': 0.0, 'coverage': 0.0, 'meets_rule': False}

    # Return the result with higher similarity
    if forward_result['similarity'] >= reverse_result['similarity']:
        return (forward_result['similarity'],
                forward_result['coverage'],
                forward_result['meets_rule'],
                "forward")
    else:
        return (reverse_result['similarity'],
                reverse_result['coverage'],
                reverse_result['meets_rule'],
                "reverse")


def find_best_family(query_seq, consensus_sequences, similarity_cutoff=0.7, coverage_cutoff=0.8):
    """Find the best matching reference family consensus for a query sequence.

    Returns: (best_family, best_similarity, best_consensus_length, best_coverage)
    """
    best_family = None
    best_similarity = 0.0
    best_coverage = 0.0
    best_consensus_length = 0

    for family_id, consensus_seq in consensus_sequences.items():
        # Convert consensus to uppercase
        consensus_seq = consensus_seq.upper()
        similarity, coverage, meets_rule, orientation = calculate_similarity_80_80_rule(
            query_seq, consensus_seq, similarity_cutoff, coverage_cutoff)

        # Update best match if this has higher similarity and meets 80-80 rule
        if meets_rule and similarity > best_similarity:
            best_similarity = similarity
            best_coverage = coverage
            best_family = family_id
            best_consensus_length = len(consensus_seq)

    # Only return family if it meets the 80-80 rule
    if best_family is not None:
        return best_family, best_similarity, best_consensus_length, best_coverage

    return None, best_similarity, 0, 0.0


def main():
    parser = argparse.ArgumentParser(
        description='Annotate missing gene and family data in genotype matrix'
    )
    parser.add_argument('--genotype_matrix', required=True,
                       help='Path to genotype_matrix_oriented.tsv')
    parser.add_argument('--processed_ann_dir', required=True,
                       help='Directory containing gene.bed files')
    parser.add_argument('--fasta_dir', required=True,
                       help='Directory containing FASTA files')
    parser.add_argument('--ref_dir', required=True,
                       help='Directory containing reference family FASTA files')
    parser.add_argument('--genome_dir', required=True,
                       help='Directory containing reference genome assemblies (e.g., graph-assemblies)')
    parser.add_argument('--output', required=True,
                       help='Output path for annotated matrix')
    parser.add_argument('--output_log', default=None,
                       help='Optional log file for family annotation results')
    parser.add_argument('--flanking_length', type=int, default=100,
                       help='Length of flanking regions to trim (default: 100)')
    parser.add_argument('--similarity_cutoff', type=float, default=0.7,
                       help='Similarity threshold for family matching (default: 0.7)')
    parser.add_argument('--coverage_cutoff', type=float, default=0.8,
                       help='Coverage threshold for family matching (default: 0.8)')
    parser.add_argument('--mafft_path', default='mafft',
                       help='Path to MAFFT executable (default: mafft)')
    parser.add_argument('--annotation_method', default='ortholog_group',
                       choices=['ortholog_group', 'sequence'],
                       help='Method for family annotation: ortholog_group (use family from other samples in group) or sequence (use consensus matching)')

    args = parser.parse_args()

    # List of samples
    samples = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749",
               "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629",
               "RCC692", "RCC693", "RCC833"]

    # Load data
    logger.info("Loading genotype matrix...")
    df = pd.read_csv(args.genotype_matrix, sep='\t')

    # Drop ortholog groups with unusable within-group status. discordant
    # groups have members at multiple distinct insertion sites within a clade
    # (should have been split but weren't), and uncertain groups have
    # insufficient data to confirm the ortholog relationship. Neither is
    # usable for downstream pi/dxy calculations.
    if 'within_group_status' in df.columns:
        drop_statuses = {'discordant', 'uncertain'}
        before_groups = df['ortholog_id'].nunique()
        before_rows = len(df)
        df = df[~df['within_group_status'].isin(drop_statuses)].copy()
        after_groups = df['ortholog_id'].nunique()
        after_rows = len(df)
        logger.info(
            f"Filtered out {before_groups - after_groups} unusable ortholog groups "
            f"({before_rows - after_rows} rows); {after_groups} groups remain")

    logger.info("Loading gene annotations...")
    gene_intervals = load_gene_bed_files(args.processed_ann_dir, samples)

    logger.info("Loading introner FASTA files...")
    fasta_dict = load_fasta_files(args.fasta_dir, samples)

    logger.info("Loading reference genome assemblies...")
    genomes = load_genome_assemblies(args.genome_dir, samples)

    logger.info("Building consensus sequences for reference families...")
    consensus_sequences = build_consensus_sequences(args.ref_dir, args.mafft_path)

    if not consensus_sequences:
        logger.error("No consensus sequences were built! Check reference directory.")
        sys.exit(1)

    logger.info(f"Built consensus sequences for {len(consensus_sequences)} families")

    # Track statistics
    genes_annotated = 0
    families_annotated = 0

    # Annotate missing genes (for ALL rows)
    logger.info("Annotating missing genes for all rows...")
    for idx, row in df.iterrows():
        if pd.isna(row['gene']) or row['gene'] == '':
            # Try to find overlapping gene
            gene = find_overlapping_gene(
                row['sample'],
                row['contig'],
                row['start'],
                row['end'],
                gene_intervals
            )
            if gene:
                df.at[idx, 'gene'] = gene
                genes_annotated += 1

    # Open log file if specified
    log_file = None
    if args.output_log:
        log_file = open(args.output_log, 'w')
        log_file.write("introner_id\tintroner_length\tbest_family\tsimilarity\tconsensus_length\n")

    # Annotate missing families (ONLY for scenario 1: presence == 1)
    logger.info(f"Annotating missing families using method: {args.annotation_method}...")
    introners_tested = 0

    for idx, row in df.iterrows():
        # Only process scenario 1 (present) introners
        if row['presence'] != 1:
            continue

        # Skip if family is already annotated
        if pd.notna(row['family']) and row['family'] != '':
            continue

        sample = row['sample']
        sequence_id = row['sequence_id']
        ortholog_id = row['ortholog_id']
        best_family = None
        best_similarity = 0.0
        best_consensus_length = 0
        introner_length = 0

        # METHOD 1: Use family from other samples in the same ortholog group
        if args.annotation_method == 'ortholog_group':
            best_family = get_family_from_ortholog_group(ortholog_id, df)
            introners_tested += 1

        # METHOD 2: Use sequence-based consensus matching
        elif args.annotation_method == 'sequence':
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])

            # Extract sequence from genome assembly using genomic coordinates
            full_seq, success = extract_sequence_from_genome(
                sample, contig, start, end, genomes, args.flanking_length
            )

            if not success or full_seq is None:
                continue

            # Extract core introner sequence by trimming flanking regions
            if len(full_seq) <= 2 * args.flanking_length:
                logger.warning(f"Sequence {sequence_id} too short for flanking length {args.flanking_length}, skipping")
                continue

            query_seq = full_seq[args.flanking_length:-args.flanking_length]
            introner_length = len(query_seq)
            introners_tested += 1

            # Find best matching family using core introner sequence against consensus
            best_family, best_similarity, best_consensus_length, best_coverage = find_best_family(
                query_seq,
                consensus_sequences,
                args.similarity_cutoff,
                args.coverage_cutoff
            )

        # Log results (only for sequence method, ortholog method is trivial)
        if log_file and args.annotation_method == 'sequence':
            if best_family is not None:
                log_file.write(f"{sequence_id}\t{introner_length}\t{best_family}\t{best_similarity:.4f}\t{best_consensus_length}\n")
            else:
                log_file.write(f"{sequence_id}\t{introner_length}\tNone\t{best_similarity:.4f}\t0\n")

        if best_family is not None:
            df.at[idx, 'family'] = best_family
            families_annotated += 1

            if families_annotated % 500 == 0:
                logger.info(f"  Annotated {families_annotated} families so far...")

    # Close log file
    if log_file:
        log_file.close()

    # Convert columns back to proper types (avoid .0 for integers and NA for strings)
    # Convert start and end back to int64 (NaN values will be handled as empty strings)
    if 'start' in df.columns:
        df['start'] = df['start'].fillna('').apply(lambda x: str(int(x)) if x != '' else '')
    if 'end' in df.columns:
        df['end'] = df['end'].fillna('').apply(lambda x: str(int(x)) if x != '' else '')
    # Convert family back to string (remove .0 from float assignments)
    if 'family' in df.columns:
        df['family'] = df['family'].fillna('').apply(lambda x: str(int(x)) if x != '' and x != 'nan' else '')

    # Save annotated matrix
    logger.info("Writing annotated matrix...")
    df.to_csv(args.output, sep='\t', index=False, na_rep='')

    # Report statistics
    logger.info("="*60)
    logger.info("ANNOTATION SUMMARY")
    logger.info("="*60)
    logger.info(f"Genes annotated: {genes_annotated}")
    logger.info(f"Families annotated: {families_annotated}")
    logger.info(f"Introners tested for family annotation: {introners_tested}")
    logger.info(f"Consensus sequences used: {len(consensus_sequences)}")
    logger.info(f"Output: {args.output}")
    if args.output_log:
        logger.info(f"Log file: {args.output_log}")
    logger.info(f"Thresholds: similarity={args.similarity_cutoff}, coverage={args.coverage_cutoff}")
    logger.info("="*60)


if __name__ == "__main__":
    main()

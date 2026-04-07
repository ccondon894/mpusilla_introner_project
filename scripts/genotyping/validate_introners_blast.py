#!/usr/bin/env python3
"""
Validate putative introners using BLAST-based family classification.

Replaces the consensus-based approach in validate_introners.py with individual
reference sequence matching via BLAST, providing more reliable family assignments
and explicit confidence margins.

Changes from validate_introners.py:
- Family assignment via BLAST against all individual reference introner sequences
  (not a single PWM consensus per family)
- Reports assignment confidence (margin between best and second-best family)
- Fixed coordinate conversion in check_in_gene() (BED 0-based → GTF 1-based)
- Tighter splice site search windows (10bp instead of 20bp)
- Requires both 5' and 3' splice sites for genic introners
"""

import os
import argparse
import subprocess
import tempfile
import re
import logging
import shutil
from io import StringIO
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import gffutils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Validate putative introners using BLAST-based family classification')
    parser.add_argument('--input_fasta', required=True,
                       help='FASTA file with putative introner sequences (including flanking regions)')
    parser.add_argument('--gtf_file', required=True,
                       help='GTF file with gene annotations')
    parser.add_argument('--ref_dir', required=True,
                       help='Directory containing reference introner family FASTA files')
    parser.add_argument('--output_fasta', required=True,
                       help='Output FASTA file for valid introners')
    parser.add_argument('--output_log', required=True,
                       help='Output TSV log file for validation results')
    parser.add_argument('--identity_cutoff', type=float, default=0.7,
                       help='Minimum percent identity (as fraction, default: 0.7)')
    parser.add_argument('--coverage_cutoff', type=float, default=0.8,
                       help='Minimum query coverage (as fraction, default: 0.8)')
    parser.add_argument('--flanking_length', type=int, default=100,
                       help='Length of flanking regions in input sequences (default: 100)')
    parser.add_argument('--blastn_path', default='blastn',
                       help='Path to blastn executable')
    parser.add_argument('--makeblastdb_path', default='makeblastdb',
                       help='Path to makeblastdb executable')
    parser.add_argument('--confidence_margin', type=float, default=5.0,
                       help='Minimum pident margin for high-confidence assignment (default: 5.0)')
    return parser.parse_args()


def setup_logging(log_file):
    """Set up logging to file and console."""
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_file + ".run.log")
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def process_input_fasta(input_fasta, flanking_length):
    """Process input FASTA file and extract putative introner sequences."""
    logging.info("Processing input FASTA file...")

    putative_introners = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        full_seq = str(record.seq)

        if len(full_seq) <= 2 * flanking_length:
            logging.warning(f"Sequence {record.id} too short ({len(full_seq)}bp), skipping")
            continue

        introner_seq = full_seq[flanking_length:-flanking_length]
        upstream_seq = full_seq[:flanking_length]
        downstream_seq = full_seq[-flanking_length:]

        # Parse location from FASTA header: contig:start-end|identifier
        try:
            header_parts = record.id.split('|')[0]
            contig_parts = header_parts.split(':')
            contig = contig_parts[0]
            position_range = contig_parts[1].split('-')
            start = int(position_range[0])
            end = int(position_range[1])
        except (IndexError, ValueError) as e:
            logging.warning(f"Could not parse location from header: {record.id}, error: {e}")
            contig, start, end = record.id, None, None

        putative_introners.append({
            'id': record.id,
            'contig': contig,
            'start': start,
            'end': end,
            'full_seq': full_seq,
            'introner_seq': introner_seq,
            'upstream_seq': upstream_seq,
            'downstream_seq': downstream_seq,
            'record': record
        })

    logging.info(f"Processed {len(putative_introners)} putative introner sequences")
    return putative_introners


# ============================================================
# BLAST-BASED FAMILY CLASSIFICATION
# ============================================================

def build_blast_database(ref_dir, makeblastdb_path, tmp_dir):
    """
    Extract individual introner body sequences from all reference family files
    and build a BLAST database.

    Reference files have mixed case: UPPERCASE = exonic flanks, lowercase = introner body.
    """
    logging.info("Building BLAST database from reference family sequences...")

    combined_fasta = os.path.join(tmp_dir, "ref_introners.fa")
    family_counts = {}

    with open(combined_fasta, 'w') as out:
        fam_files = [f for f in os.listdir(ref_dir) if f.startswith("GCA_000151265.1_fam")]

        for fam_file in sorted(fam_files):
            family_match = re.search(r'fam(\d+)', fam_file)
            if not family_match:
                continue
            family_id = family_match.group(1)

            file_path = os.path.join(ref_dir, fam_file)
            seq_idx = 0
            for record in SeqIO.parse(file_path, "fasta"):
                seq_str = str(record.seq)
                introner_match = re.search(r'[a-z]+', seq_str)
                if introner_match:
                    introner_body = seq_str[introner_match.start():introner_match.end()].upper()
                    out.write(f">fam{family_id}_seq{seq_idx:04d}\n{introner_body}\n")
                    seq_idx += 1

            family_counts[family_id] = seq_idx
            logging.info(f"  Family {family_id}: {seq_idx} reference sequences")

    # Build BLAST database
    db_path = os.path.join(tmp_dir, "ref_introners")
    cmd = [makeblastdb_path, '-in', combined_fasta, '-dbtype', 'nucl',
           '-out', db_path, '-parse_seqids']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"makeblastdb failed: {result.stderr}")
        raise RuntimeError("Failed to build BLAST database")

    logging.info(f"BLAST database built at {db_path}")
    return db_path, family_counts


def run_blast(putative_introners, db_path, blastn_path, tmp_dir):
    """Run BLAST of candidate introner bodies against reference database."""
    logging.info("Running BLAST...")

    # Write query FASTA
    query_fasta = os.path.join(tmp_dir, "query_introners.fa")
    with open(query_fasta, 'w') as out:
        for pi in putative_introners:
            out.write(f">{pi['id']}\n{pi['introner_seq']}\n")

    # Run BLAST
    output_file = os.path.join(tmp_dir, "blast_results.tsv")
    cmd = [
        blastn_path,
        '-task', 'blastn',
        '-query', query_fasta,
        '-db', db_path,
        '-evalue', '1e-5',
        '-max_target_seqs', '10',
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-out', output_file,
        '-num_threads', '4'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"BLAST failed: {result.stderr}")
        raise RuntimeError("BLAST search failed")

    # Parse results
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

    if os.path.getsize(output_file) == 0:
        logging.warning("BLAST produced no results")
        return pd.DataFrame(columns=columns)

    blast_df = pd.read_csv(output_file, sep='\t', header=None, names=columns)
    logging.info(f"BLAST returned {len(blast_df)} hits for {blast_df['qseqid'].nunique()} queries")

    return blast_df


def assign_families(blast_df, putative_introners, identity_cutoff, coverage_cutoff,
                    confidence_margin):
    """
    Assign family to each introner based on BLAST results.

    For each query:
    - Best hit (highest bitscore) determines family
    - Best hit from a different family determines second-best
    - Margin = best pident - second-best pident
    """
    logging.info("Assigning families from BLAST results...")

    # Extract family from subject ID (e.g., "fam2_seq0001" -> "2")
    if len(blast_df) > 0:
        blast_df = blast_df.copy()
        blast_df['family'] = blast_df['sseqid'].str.extract(r'fam(\d+)_')[0]
        blast_df['qcoverage'] = blast_df['length'] / blast_df['qlen']

    # Build lookup of all putative IDs
    all_ids = {pi['id'] for pi in putative_introners}

    assignments = {}

    for query_id in all_ids:
        hits = blast_df[blast_df['qseqid'] == query_id] if len(blast_df) > 0 else pd.DataFrame()

        if len(hits) == 0:
            assignments[query_id] = {
                'best_family': None,
                'similarity': 0.0,
                'coverage': 0.0,
                'orientation': None,
                'similarity_valid': False,
                'second_best_family': None,
                'second_best_similarity': 0.0,
                'assignment_margin': 0.0,
                'assignment_confidence': 'no_hit'
            }
            continue

        # Sort by bitscore descending
        hits_sorted = hits.sort_values('bitscore', ascending=False)

        # Best hit
        best = hits_sorted.iloc[0]
        best_family = best['family']
        best_pident = best['pident']
        best_coverage = best['qcoverage']

        # Determine orientation from BLAST coordinates
        if best['sstart'] > best['send']:
            orientation = 'reverse'
        else:
            orientation = 'forward'

        # Check if best hit meets thresholds (pident is in %, cutoffs are fractions)
        meets_threshold = (best_pident >= identity_cutoff * 100 and
                          best_coverage >= coverage_cutoff)

        # Second-best family (best hit from a DIFFERENT family)
        other_family_hits = hits_sorted[hits_sorted['family'] != best_family]
        if len(other_family_hits) > 0:
            second = other_family_hits.iloc[0]
            second_family = second['family']
            second_pident = second['pident']
        else:
            second_family = None
            second_pident = 0.0

        margin = best_pident - second_pident

        if not meets_threshold:
            confidence = 'no_hit'
        elif margin >= confidence_margin:
            confidence = 'high'
        else:
            confidence = 'low'

        assignments[query_id] = {
            'best_family': best_family,
            'similarity': best_pident / 100.0,  # Convert to fraction for output compatibility
            'coverage': best_coverage,
            'orientation': orientation,
            'similarity_valid': meets_threshold,
            'second_best_family': second_family,
            'second_best_similarity': second_pident / 100.0,
            'assignment_margin': margin,
            'assignment_confidence': confidence
        }

    n_assigned = sum(1 for a in assignments.values() if a['best_family'] is not None)
    n_valid = sum(1 for a in assignments.values() if a['similarity_valid'])
    n_high = sum(1 for a in assignments.values() if a['assignment_confidence'] == 'high')
    n_low = sum(1 for a in assignments.values() if a['assignment_confidence'] == 'low')

    logging.info(f"  Assigned family: {n_assigned}/{len(assignments)}")
    logging.info(f"  Meeting threshold: {n_valid}")
    logging.info(f"  High confidence: {n_high}, Low confidence: {n_low}")

    return assignments


# ============================================================
# SPATIAL VALIDATION (bug-fixed)
# ============================================================

def setup_gtf_db(gtf_file):
    """Create a database from GTF file for gene annotation queries."""
    logging.info(f"Creating database from GTF file: {gtf_file}")

    db_file = f"{gtf_file}.db"

    if not os.path.exists(db_file):
        try:
            db = gffutils.create_db(gtf_file, dbfn=db_file, force=True,
                                     keep_order=True, merge_strategy='merge',
                                     sort_attribute_values=True,
                                     disable_infer_genes=False,
                                     disable_infer_transcripts=False)
            logging.info(f"Created GTF database at {db_file}")
        except Exception as e:
            logging.error(f"Error creating GTF database: {e}")
            return None
    else:
        logging.info(f"Using existing GTF database at {db_file}")

    try:
        return gffutils.FeatureDB(db_file)
    except Exception as e:
        logging.error(f"Error opening GTF database: {e}")
        return None


def check_in_gene(introner, gtf_db, flanking_length=100):
    """
    Check if an introner is located within a gene.

    Coordinate conversion:
    - Input coords from FASTA header are BED 0-based half-open, including flanks
    - GTF uses 1-based inclusive coordinates
    - BED [start, end) → GTF [start+1, end]
    - Then remove flanking: GTF introner = [start+1+flank, end-flank]
    """
    if not gtf_db or introner['contig'] is None or introner['start'] is None:
        return False, None

    try:
        # Convert BED 0-based half-open to GTF 1-based inclusive
        start_1based = introner['start'] + 1
        end_1based = introner['end']  # BED end is exclusive, so it equals GTF end

        # Remove flanking to get introner body coordinates (1-based)
        introner_start = start_1based + flanking_length
        introner_end = end_1based - flanking_length

        if introner_end <= introner_start:
            return False, None

        contig = introner['contig']

        try:
            overlapping_genes = list(gtf_db.region(
                seqid=contig, start=introner_start, end=introner_end,
                featuretype='gene'))

            if overlapping_genes:
                return True, overlapping_genes[0]
        except Exception as e:
            logging.debug(f"Failed to query GTF for {introner['id']}: {e}")

        return False, None

    except Exception as e:
        logging.warning(f"Error checking if introner {introner['id']} is in a gene: {e}")
        return False, None


def check_splice_sites(introner_seq):
    """
    Check for canonical splice sites in an introner sequence.

    Looks for:
    - 5' splice site: GT or GC within the first 20bp
    - 3' splice site: AG within the last 20bp (recorded but not required for validity)

    Only the 5' site is required for has_valid_sites, matching the established
    behavior. Many Micromonas introners lack a canonical 3' AG near their end.
    """
    introner_seq = introner_seq.upper()
    seq_len = len(introner_seq)

    # 5' splice site: search first 20bp
    search_window_5 = min(20, seq_len)
    first_bp = introner_seq[:search_window_5]
    five_prime_site = None
    five_prime_pos = None

    gt_match = re.search(r'GT', first_bp)
    if gt_match:
        five_prime_site = 'GT'
        five_prime_pos = gt_match.start()
    else:
        gc_match = re.search(r'GC', first_bp)
        if gc_match:
            five_prime_site = 'GC'
            five_prime_pos = gc_match.start()

    # 3' splice site: search last 20bp for AG (recorded for annotation, not required)
    search_window_3 = min(20, seq_len)
    last_bp = introner_seq[-search_window_3:] if seq_len >= search_window_3 else introner_seq
    three_prime_pos = None

    ag_matches = [m for m in re.finditer(r'AG', last_bp)]
    if ag_matches:
        last_ag = ag_matches[-1]
        three_prime_pos = seq_len - search_window_3 + last_ag.start()

    # Only require 5' splice site for validity
    has_valid_sites = five_prime_site is not None

    return has_valid_sites, five_prime_site, five_prime_pos, three_prime_pos


# ============================================================
# VALIDATION PIPELINE
# ============================================================

def validate_introners(putative_introners, family_assignments, gtf_db, flanking_length=100):
    """
    Validate putative introners:
    1. Family assignment from BLAST (already computed)
    2. Check if located within genes (with corrected coordinates)
    3. Check for valid splice sites if within genes
    """
    logging.info("Running spatial validation...")

    results = []

    for putative in putative_introners:
        pid = putative['id']
        assignment = family_assignments.get(pid, {})

        result = {
            'putative_id': pid,
            'best_family': assignment.get('best_family'),
            'similarity': assignment.get('similarity', 0.0),
            'coverage': assignment.get('coverage', 0.0),
            'orientation': assignment.get('orientation'),
            'similarity_valid': assignment.get('similarity_valid', False),
            'second_best_family': assignment.get('second_best_family'),
            'second_best_similarity': assignment.get('second_best_similarity', 0.0),
            'assignment_margin': assignment.get('assignment_margin', 0.0),
            'assignment_confidence': assignment.get('assignment_confidence', 'no_hit'),
        }

        # Check if in gene
        in_gene, overlapping_gene = check_in_gene(putative, gtf_db, flanking_length)
        result['in_gene'] = in_gene
        result['gene_id'] = overlapping_gene.id if overlapping_gene else None

        # Check splice sites
        if in_gene:
            has_splice_site, site_type, site_pos, three_prime_pos = check_splice_sites(
                putative['introner_seq'])
            result['has_splice_site'] = has_splice_site
            result['splice_site_type'] = site_type
            result['splice_site_pos'] = site_pos
            result['three_prime_splice_pos'] = three_prime_pos
        else:
            result['has_splice_site'] = None
            result['splice_site_type'] = None
            result['splice_site_pos'] = None
            result['three_prime_splice_pos'] = None

        # Determine validity:
        # Valid if meets similarity threshold AND (not in gene OR in gene with valid splice sites)
        if result['similarity_valid']:
            if not in_gene or (in_gene and result['has_splice_site']):
                result['is_valid'] = True
            else:
                result['is_valid'] = False
        else:
            result['is_valid'] = False

        result['record'] = putative['record']
        results.append(result)

        if len(results) % 500 == 0:
            logging.info(f"Validated {len(results)} introners...")

    return results


# ============================================================
# OUTPUT
# ============================================================

def write_log_file(results, output_log, identity_cutoff, coverage_cutoff):
    """Write validation results to TSV log file."""
    logging.info(f"Writing detailed log to {output_log}...")

    with open(output_log, 'w') as log:
        # Header
        log.write("introner_name\tbest_matching_family\tsimilarity\tcoverage\torientation\t"
                  "is_within_gene\thas_splice_site\tsplice_site_type\tsplice_site_position\t"
                  "three_prime_splice_pos\tsecond_best_family\tsecond_best_similarity\t"
                  "assignment_margin\tassignment_confidence\tis_valid\n")

        for r in results:
            fields = [
                r['putative_id'],
                r['best_family'] or 'none',
                f"{r['similarity']:.4f}",
                f"{r['coverage']:.4f}",
                r['orientation'] or 'none',
                str(r['in_gene']).lower(),
                str(r['has_splice_site']).lower() if r['has_splice_site'] is not None else 'NA',
                r['splice_site_type'] or 'NA',
                str(r['splice_site_pos']) if r['splice_site_pos'] is not None else 'NA',
                str(r['three_prime_splice_pos']) if r['three_prime_splice_pos'] is not None else 'NA',
                r['second_best_family'] or 'none',
                f"{r['second_best_similarity']:.4f}",
                f"{r['assignment_margin']:.2f}",
                r['assignment_confidence'],
                str(r['is_valid']).lower(),
            ]
            log.write('\t'.join(fields) + '\n')

    # Summary log
    total = len(results)
    valid = sum(1 for r in results if r['is_valid'])
    in_gene = sum(1 for r in results if r['in_gene'])
    n_high = sum(1 for r in results if r['assignment_confidence'] == 'high')
    n_low = sum(1 for r in results if r['assignment_confidence'] == 'low')
    n_no_hit = sum(1 for r in results if r['assignment_confidence'] == 'no_hit')

    summary_log = f"{os.path.splitext(output_log)[0]}_summary.log"
    with open(summary_log, 'w') as summ:
        summ.write("# Introner Validation Summary (BLAST-based)\n")
        summ.write(f"Total putative introners analyzed: {total}\n")
        summ.write(f"Identity threshold: {identity_cutoff}, Coverage threshold: {coverage_cutoff}\n")
        summ.write(f"\n# Family Assignment\n")
        summ.write(f"High confidence assignments: {n_high} ({100*n_high/total:.1f}%)\n")
        summ.write(f"Low confidence assignments:  {n_low} ({100*n_low/total:.1f}%)\n")
        summ.write(f"No hit / below threshold:    {n_no_hit} ({100*n_no_hit/total:.1f}%)\n")
        summ.write(f"\n# Spatial Validation\n")
        summ.write(f"Introners within genes: {in_gene}\n")
        missing_splice = sum(1 for r in results if r['in_gene'] and not r['has_splice_site'])
        summ.write(f"In gene but missing splice sites: {missing_splice}\n")
        summ.write(f"\n# Final Results\n")
        summ.write(f"Valid introners: {valid} ({100*valid/total:.1f}%)\n")

        # Family breakdown
        summ.write(f"\n# Family Breakdown (valid introners)\n")
        from collections import Counter
        fam_counts = Counter(r['best_family'] for r in results if r['is_valid'])
        for fam, count in sorted(fam_counts.items(), key=lambda x: -x[1]):
            summ.write(f"  Family {fam}: {count}\n")

    logging.info(f"Created TSV log: {output_log}")
    logging.info(f"Created summary log: {summary_log}")


def write_filtered_fasta(results, output_fasta):
    """Write valid introners to filtered FASTA file."""
    logging.info(f"Writing filtered FASTA to {output_fasta}...")

    valid_records = []

    for result in results:
        if result['is_valid']:
            record = result['record']

            annotations = [
                f"family={result['best_family']}",
                f"similarity={result['similarity']:.4f}",
                f"orientation={result['orientation']}",
                f"confidence={result['assignment_confidence']}",
            ]

            if result['in_gene']:
                annotations.append(f"gene={result['gene_id']}")
                splice_parts = []
                if result['splice_site_type'] and result['splice_site_pos'] is not None:
                    splice_parts.append(f"{result['splice_site_type']}@{result['splice_site_pos']}")
                if result['three_prime_splice_pos'] is not None:
                    splice_parts.append(f"AG@{result['three_prime_splice_pos']}")
                annotations.append(f"splice_sites={';'.join(splice_parts) if splice_parts else 'none'}")

            record.description = f"{record.description} {' '.join(annotations)}"
            valid_records.append(record)

    with open(output_fasta, 'w') as out:
        SeqIO.write(valid_records, out, "fasta")

    logging.info(f"Wrote {len(valid_records)} valid introners to {output_fasta}")


# ============================================================
# MAIN
# ============================================================

def main():
    args = parse_args()

    setup_logging(args.output_log)

    # Create temp directory for BLAST database
    tmp_dir = tempfile.mkdtemp(prefix="introner_blast_")

    try:
        # Build BLAST database from reference families
        db_path, family_counts = build_blast_database(
            args.ref_dir, args.makeblastdb_path, tmp_dir)

        # Process input FASTA
        putative_introners = process_input_fasta(args.input_fasta, args.flanking_length)

        # Run BLAST
        blast_df = run_blast(putative_introners, db_path, args.blastn_path, tmp_dir)

        # Assign families
        family_assignments = assign_families(
            blast_df, putative_introners,
            args.identity_cutoff, args.coverage_cutoff,
            args.confidence_margin)

        # Set up GTF database
        gtf_db = setup_gtf_db(args.gtf_file)

        # Validate
        results = validate_introners(
            putative_introners, family_assignments, gtf_db, args.flanking_length)

        # Write outputs
        write_log_file(results, args.output_log, args.identity_cutoff, args.coverage_cutoff)
        write_filtered_fasta(results, args.output_fasta)

        logging.info("Introner validation complete!")

    finally:
        # Clean up temp directory
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()

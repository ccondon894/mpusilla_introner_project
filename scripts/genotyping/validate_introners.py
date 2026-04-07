#!/usr/bin/env python3

import os
import argparse
import subprocess
import tempfile
import re
from io import StringIO
import logging
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import gffutils

def parse_args():
    parser = argparse.ArgumentParser(description='Validate putative introners against reference family consensus sequences')
    parser.add_argument('--input_fasta', required=True, help='FASTA file with putative introner sequences (including flanking regions)')
    parser.add_argument('--gtf_file', required=True, help='GTF file with gene annotations')
    parser.add_argument('--ref_dir', required=True, help='Directory containing reference introner family FASTA files')
    parser.add_argument('--output_fasta', required=True, help='Output FASTA file for valid introners')
    parser.add_argument('--output_log', required=True, help='Output log file for validation results')
    parser.add_argument('--similarity_cutoff', type=float, default=0.7, help='Similarity threshold (default: 0.7)')
    parser.add_argument('--coverage_cutoff', type=float, default=0.8, help='Coverage threshold (default: 0.8)')
    parser.add_argument('--mafft_path', default='mafft', help='Path to MAFFT executable')
    parser.add_argument('--flanking_length', type=int, default=100, help='Length of flanking regions in the input sequences (default: 100)')
    parser.add_argument('--disable_structural_analysis', action='store_true', help='Skip structural analysis (no-op if structural analysis is not implemented)')
    parser.add_argument('--disable_micromonas_patterns', action='store_true', help='Skip Micromonas-specific pattern checks (no-op if not implemented)')

    args = parser.parse_args()

    return args

def setup_logging(log_file):
    """Set up logging to file and console"""
    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Create file handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def process_input_fasta(input_fasta, flanking_length):
    """Process input FASTA file and extract putative introner sequences"""
    logging.info("Processing input FASTA file...")
    
    putative_introners = []
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        full_seq = str(record.seq)
        
        # Skip if sequence is too short to have flanking regions
        if len(full_seq) <= 2 * flanking_length:
            logging.warning(f"Warning: Sequence {record.id} is too short ({len(full_seq)}bp), skipping")
            continue
            
        # Extract introner and flanking sequences
        introner_seq = full_seq[flanking_length:-flanking_length]
        upstream_seq = full_seq[:flanking_length]
        downstream_seq = full_seq[-flanking_length:]
        
        # Parse location from FASTA header format: CCMP1545#0#scaffold_3:288947-289326|introner_seq_3514
        try:
            # Split by '|' to get the coordinate part
            header_parts = record.id.split('|')[0]
            
            # Split by ':' to get contig and position range
            contig_parts = header_parts.split(':')
            contig = contig_parts[0]
            
            
            # Split the position range by '-'
            position_range = contig_parts[1].split('-')
            start = int(position_range[0])
            end = int(position_range[1])
            
            logging.debug(f"Parsed location for {record.id}: contig={contig}, start={start}, end={end}")
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
            'record': record  # Keep original record for output
        })
    
    logging.info(f"Processed {len(putative_introners)} putative introner sequences")
    return putative_introners

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

def run_mafft_alignment(sequences, mafft_path):
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
        logging.error(f"MAFFT alignment failed: {e}")
        logging.error(f"STDERR: {e.stderr}")
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

def build_consensus_sequences(ref_dir, mafft_path):
    """Build consensus sequences for each reference family"""
    logging.info("Building consensus sequences for reference families...")
    
    consensus_sequences = {}
    
    # Find all family FASTA files matching GCA_000151265.1_fam pattern
    fam_files = [f for f in os.listdir(ref_dir) if f.startswith("GCA_000151265.1_fam")]
    
    if not fam_files:
        logging.warning(f"No family FASTA files found in {ref_dir}")
        return consensus_sequences
    
    for fam_file in fam_files:
        # Extract family ID from filename
        family_match = re.search(r'fam(\d+)', fam_file)
        if family_match:
            family_id = family_match.group(1)
        else:
            logging.warning(f"Could not extract family ID from {fam_file}, skipping")
            continue
        
        # Extract introner sequences
        introner_seqs = extract_introner_sequences(ref_dir, fam_file)
        
        if len(introner_seqs) < 2:
            logging.warning(f"Family {family_id} has only {len(introner_seqs)} sequences, skipping alignment")
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
            
            logging.info(f"Generated consensus for family {family_id} (length: {len(consensus)})")
        else:
            logging.error(f"Failed to align sequences for family {family_id}")
    
    return consensus_sequences

def calculate_similarity_80_80_rule(seq1, seq2, similarity_cutoff=0.7, coverage_cutoff=0.8):
    """
    Calculate similarity between two sequences:
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

def setup_gtf_db(gtf_file):
    """Create a database from GTF file for gene annotation queries"""
    logging.info(f"Creating database from GTF file: {gtf_file}")
    
    # Create a temporary DB file
    db_file = f"{gtf_file}.db"
    
    # Check if database already exists
    if not os.path.exists(db_file):
        try:
            # Create database from GTF file
            db = gffutils.create_db(gtf_file, 
                                   dbfn=db_file, 
                                   force=True, 
                                   keep_order=True,
                                   merge_strategy='merge',
                                   sort_attribute_values=True,
                                   disable_infer_genes=False,
                                   disable_infer_transcripts=False)
            logging.info(f"Created GTF database at {db_file}")
        except Exception as e:
            logging.error(f"Error creating GTF database: {e}")
            return None
    else:
        logging.info(f"Using existing GTF database at {db_file}")
    
    # Open the database
    try:
        db = gffutils.FeatureDB(db_file)
        return db
    except Exception as e:
        logging.error(f"Error opening GTF database: {e}")
        return None

def check_in_gene(introner, gtf_db, flanking_length=100):
    """Check if an introner is located within a gene"""
    if not gtf_db or introner['contig'] is None or introner['start'] is None or introner['end'] is None:
        return False, None
    
    try:
        # Calculate actual introner coordinates by removing flanking regions
        # The coordinates from the header include flanking regions on each side
        introner_start = introner['start'] + flanking_length
        introner_end = introner['end'] - flanking_length
        
        # The contig name in GTF might not include the scaffold/chromosome prefix
        # Try with and without any prefix
        contig = introner['contig']
        
        # Handle potential differences in contig naming between FASTA and GTF
        contig_variations = [
            contig,
            contig.replace('scaffold_', ''),
            contig.replace('chromosome_', ''),
            contig.replace('chr', ''),
            f"chr{contig}" if not contig.startswith('chr') else contig,
            f"scaffold_{contig}" if not contig.startswith('scaffold_') else contig
        ]
        
        # Try each contig variation
        for contig_var in contig_variations:
            try:
                # Query for genes that overlap with the introner
                logging.debug(f"Searching for genes overlapping {contig_var}:{introner_start}-{introner_end}")
                overlapping_genes = list(gtf_db.region(seqid=contig_var, 
                                                    start=introner_start, 
                                                    end=introner_end, 
                                                    featuretype='gene'))
                
                if overlapping_genes:
                    # Return True and the first overlapping gene
                    logging.debug(f"Found {len(overlapping_genes)} overlapping genes for {introner['id']}")
                    return True, overlapping_genes[0]
            except Exception as e:
                logging.debug(f"Failed to query with contig '{contig_var}': {e}")
                continue
        
        # If we get here, no overlapping genes were found with any contig variation
        return False, None
        
    except Exception as e:
        logging.warning(f"Error checking if introner {introner['id']} is in a gene: {e}")
        return False, None

def check_splice_sites(introner_seq):
    """
    Check if an introner sequence has valid splice sites.
    Looks for 5' splice site (GT/GC at beginning) and 3' splice site (AG at end).
    """
    # Convert to uppercase for consistency
    introner_seq = introner_seq.upper()
    seq_len = len(introner_seq)
    
    # Check 5' splice site (first 20bp)
    first_20bp = introner_seq[:20]
    five_prime_site = None
    five_prime_pos = None
    
    # Check for GT (most common)
    gt_match = re.search(r'GT', first_20bp)
    if gt_match:
        five_prime_site = 'GT'
        five_prime_pos = gt_match.start()
    else:
        # Check for GC (less common)
        gc_match = re.search(r'GC', first_20bp)
        if gc_match:
            five_prime_site = 'GC'
            five_prime_pos = gc_match.start()
    
    # Check 3' splice site (last 20bp, looking for AG)
    last_20bp = introner_seq[-20:] if seq_len >= 20 else introner_seq
    three_prime_pos = None
    
    # Look for AG near the end (working backwards from end)
    ag_matches = [m for m in re.finditer(r'AG', last_20bp)]
    if ag_matches:
        # Take the last AG match (closest to 3' end)
        last_ag = ag_matches[-1]
        three_prime_pos = seq_len - (20 - last_ag.start()) if seq_len >= 20 else last_ag.start()
    
    # Determine if we have valid splice sites
    has_valid_sites = five_prime_site is not None
    
    # Return results - including 3' splice site info for structural analysis
    return has_valid_sites, five_prime_site, five_prime_pos, three_prime_pos

def validate_introners(putative_introners, consensus_sequences, gtf_db, similarity_cutoff=0.7, coverage_cutoff=0.8, flanking_length=100):
    """
    Validate putative introners:
    1. Compare against family consensus sequences
    2. Check if located within genes
    3. Check for valid splice sites if within genes
    """
    logging.info(f"Validating introners (similarity cutoff: {similarity_cutoff}, coverage cutoff: {coverage_cutoff})...")
    
    results = []
    
    for putative in putative_introners:
        putative_id = putative['id']
        putative_seq = putative['introner_seq'].lower()  # Convert to lowercase for comparison
        
        # 1. Compare against consensus sequences
        best_match = {
            'putative_id': putative_id,
            'best_family': None,
            'similarity': 0.0,
            'coverage': 0.0,
            'orientation': None,
            'similarity_valid': False
        }
        
        # Compare against each family consensus
        for family_id, consensus_seq in consensus_sequences.items():
            # Convert consensus to lowercase to match putative sequence
            consensus_seq = consensus_seq.lower()
            similarity, coverage, meets_rule, orientation = calculate_similarity_80_80_rule(
                putative_seq, consensus_seq, similarity_cutoff, coverage_cutoff)
            
            # Update best match if this has higher similarity
            if similarity > best_match['similarity']:
                best_match.update({
                    'best_family': family_id,
                    'similarity': similarity,
                    'coverage': coverage,
                    'orientation': orientation,
                    'similarity_valid': meets_rule
                })
        
        # 2. Check if located within a gene
        in_gene, overlapping_gene = check_in_gene(putative, gtf_db, flanking_length)
        best_match['in_gene'] = in_gene
        best_match['gene_id'] = overlapping_gene.id if overlapping_gene else None
        
        # 3. Check for valid splice sites if within a gene
        if in_gene:
            has_splice_site, site_type, site_pos, three_prime_pos = check_splice_sites(putative_seq)
            best_match['has_splice_site'] = has_splice_site
            best_match['splice_site_type'] = site_type
            best_match['splice_site_pos'] = site_pos
            best_match['three_prime_splice_pos'] = three_prime_pos
        else:
            best_match['has_splice_site'] = None
            best_match['splice_site_type'] = None
            best_match['splice_site_pos'] = None
            best_match['three_prime_splice_pos'] = None
        
        # Structural analysis removed - only family classification is used

        # 5. Determine overall validity
        # Valid if:
        # - Meets similarity threshold AND
        # - Either not in a gene OR (in a gene AND has valid splice site)
        if best_match['similarity_valid']:
            if not in_gene or (in_gene and best_match['has_splice_site']):
                best_match['is_valid'] = True
            else:
                best_match['is_valid'] = False
        else:
            best_match['is_valid'] = False
        
        # Store original record for output
        best_match['record'] = putative['record']
        
        results.append(best_match)
        
        # Print progress for every 100 introners
        if len(results) % 100 == 0:
            logging.info(f"Processed {len(results)} introners...")
    
    return results

def write_log_file(results, output_log, similarity_cutoff, coverage_cutoff):
    """Write validation results to log file"""
    logging.info(f"Writing detailed log to {output_log}...")
    
    # Calculate summary statistics
    total_count = len(results)
    similarity_invalid_count = sum(1 for r in results if not r['similarity_valid'])
    in_gene_count = sum(1 for r in results if r['in_gene'])
    missing_splice_site_count = sum(1 for r in results if r['in_gene'] and not r['has_splice_site'])
    valid_basic_count = sum(1 for r in results if r.get('is_valid_basic', False))
    valid_enhanced_count = sum(1 for r in results if r.get('is_valid_enhanced', False))
    valid_count = sum(1 for r in results if r['is_valid'])
    
    # Structural analysis statistics
    has_tir_count = sum(1 for r in results if r.get('has_tir', False))
    has_tsd_count = sum(1 for r in results if r.get('has_tsd', False))
    has_both_count = sum(1 for r in results if r.get('has_tir', False) and r.get('has_tsd', False))
    
    # Classification statistics
    dna_transposon_count = sum(1 for r in results if r.get('introner_classification') == 'DNA_transposon_like')
    tir_only_count = sum(1 for r in results if r.get('introner_classification') == 'TIR_only')
    tsd_only_count = sum(1 for r in results if r.get('introner_classification') == 'TSD_only')
    no_features_count = sum(1 for r in results if r.get('introner_classification') == 'no_classical_features')
    
    # Create detailed TSV log file
    with open(output_log, 'w') as log:
        # Write header
        log.write("introner_name\tbest_matching_family\tsimilarity\tcoverage\torientation\tis_within_gene\thas_splice_site\tsplice_site_type\tsplice_site_position\thas_tir\ttir_length\ttir_confidence\thas_tsd\ttsd_length\ttsd_sequence\tintroner_classification\tsplice_site_recruitment\tstructure_confidence\tis_valid_basic\tis_valid_enhanced\tis_valid\n")
        
        # Write detailed information for each introner
        for r in results:
            log.write(f"{r['putative_id']}\t")  # introner_name
            log.write(f"{r['best_family'] or 'none'}\t")  # best_matching_family
            log.write(f"{r['similarity']:.4f}\t")  # similarity
            log.write(f"{r['coverage']:.4f}\t")  # coverage
            log.write(f"{r['orientation'] or 'none'}\t")  # orientation
            log.write(f"{str(r['in_gene']).lower()}\t")  # is_within_gene
            log.write(f"{str(r['has_splice_site']).lower() if r['has_splice_site'] is not None else 'NA'}\t")  # has_splice_site
            log.write(f"{r['splice_site_type'] or 'NA'}\t")  # splice_site_type
            log.write(f"{r['splice_site_pos'] if r['splice_site_pos'] is not None else 'NA'}\t")  # splice_site_position
            # Structural analysis columns
            log.write(f"{str(r.get('has_tir', 'NA')).lower()}\t")  # has_tir
            log.write(f"{r.get('tir_length', 'NA')}\t")  # tir_length
            log.write(f"{r.get('tir_confidence', 'NA')}\t")  # tir_confidence
            log.write(f"{str(r.get('has_tsd', 'NA')).lower()}\t")  # has_tsd
            log.write(f"{r.get('tsd_length', 'NA')}\t")  # tsd_length
            log.write(f"{r.get('tsd_sequence', 'none')}\t")  # tsd_sequence
            log.write(f"{r.get('introner_classification', 'NA')}\t")  # introner_classification
            log.write(f"{r.get('splice_site_recruitment', 'NA')}\t")  # splice_site_recruitment
            log.write(f"{r.get('structure_confidence', 'NA')}\t")  # structure_confidence
            log.write(f"{str(r.get('is_valid_basic', 'NA')).lower()}\t")  # is_valid_basic
            log.write(f"{str(r.get('is_valid_enhanced', 'NA')).lower()}\t")  # is_valid_enhanced
            log.write(f"{str(r['is_valid']).lower()}\n")  # is_valid
    
    # Create a summary log file
    summary_log = f"{os.path.splitext(output_log)[0]}_summary.log"
    with open(summary_log, 'w') as summ:
        # Write summary
        summ.write("# Introner Validation Summary\n")
        summ.write(f"Total putative introners analyzed: {total_count}\n")
        summ.write(f"Similarity threshold: {similarity_cutoff}, Coverage threshold: {coverage_cutoff}\n")
        summ.write(f"Introners below similarity threshold: {similarity_invalid_count}\n")
        summ.write(f"Introners within genes: {in_gene_count}\n")
        summ.write(f"Introners within genes missing valid splice sites: {missing_splice_site_count}\n")
        summ.write("\n# Validation Results\n")
        summ.write(f"Valid introners (basic criteria): {valid_basic_count}\n")
        summ.write(f"Valid introners (enhanced criteria): {valid_enhanced_count}\n")
        summ.write(f"Valid introners (legacy): {valid_count}\n")
        summ.write("\n# Structural Feature Distribution\n")
        summ.write(f"Introners with TIRs: {has_tir_count} ({100*has_tir_count/total_count:.1f}%)\n")
        summ.write(f"Introners with TSDs: {has_tsd_count} ({100*has_tsd_count/total_count:.1f}%)\n")
        summ.write(f"Introners with both TIRs and TSDs: {has_both_count} ({100*has_both_count/total_count:.1f}%)\n")
        summ.write("\n# Classification Breakdown\n")
        summ.write(f"DNA transposon-like: {dna_transposon_count} ({100*dna_transposon_count/total_count:.1f}%)\n")
        summ.write(f"TIR only: {tir_only_count} ({100*tir_only_count/total_count:.1f}%)\n")
        summ.write(f"TSD only: {tsd_only_count} ({100*tsd_only_count/total_count:.1f}%)\n")
        summ.write(f"No classical features: {no_features_count} ({100*no_features_count/total_count:.1f}%)\n")
    
    logging.info(f"Created TSV log file: {output_log}")
    logging.info(f"Created summary log file: {summary_log}")

def write_filtered_fasta(results, output_fasta):
    """Write valid introners to filtered FASTA file"""
    logging.info(f"Writing filtered FASTA file to {output_fasta}...")
    
    valid_records = []
    
    for result in results:
        if result['is_valid']:
            # Get the original record
            record = result['record']
            
            # Add validation info to the description
            family_info = f"family={result['best_family']}"
            similarity_info = f"similarity={result['similarity']:.4f}"
            orientation_info = f"orientation={result['orientation']}"
            
            # Add structural analysis info
            tir_info = f"tir={result.get('has_tir', 'NA')}:{result.get('tir_length', 'NA')}:{result.get('tir_confidence', 'NA')}"
            tsd_info = f"tsd={result.get('has_tsd', 'NA')}:{result.get('tsd_length', 'NA')}:{result.get('tsd_sequence', 'none')}"
            classification_info = f"class={result.get('introner_classification', 'NA')}"
            structure_confidence_info = f"struct_conf={result.get('structure_confidence', 'NA')}"
            validation_info = f"valid_basic={result.get('is_valid_basic', 'NA')}:enhanced={result.get('is_valid_enhanced', 'NA')}"
            
            # Add gene info if available
            if result['in_gene']:
                gene_info = f"gene={result['gene_id']}"
                
                # Build splice site info with both 5' and 3' positions
                splice_parts = []
                if result['splice_site_type'] and result['splice_site_pos'] is not None:
                    splice_parts.append(f"{result['splice_site_type']}@{result['splice_site_pos']}")
                if result['three_prime_splice_pos'] is not None:
                    splice_parts.append(f"AG@{result['three_prime_splice_pos']}")
                
                splice_info = f"splice_sites={';'.join(splice_parts) if splice_parts else 'none'}"
                record.description = f"{record.description} {family_info} {similarity_info} {orientation_info} {tir_info} {tsd_info} {classification_info} {structure_confidence_info} {validation_info} {gene_info} {splice_info}"
            else:
                record.description = f"{record.description} {family_info} {similarity_info} {orientation_info} {tir_info} {tsd_info} {classification_info} {structure_confidence_info} {validation_info}"
            
            valid_records.append(record)
    
    # Write valid records to output file
    with open(output_fasta, 'w') as out:
        SeqIO.write(valid_records, out, "fasta")
    
    logging.info(f"Wrote {len(valid_records)} valid introners to {output_fasta}")

def main():
    args = parse_args()
    
    # Set up logging
    setup_logging(args.output_log)
    
    # Build consensus sequences for reference families
    consensus_sequences = build_consensus_sequences(args.ref_dir, args.mafft_path)
    
    # Process input FASTA file
    putative_introners = process_input_fasta(args.input_fasta, args.flanking_length)
    
    # Set up GTF database
    gtf_db = setup_gtf_db(args.gtf_file)
    
    # Validate putative introners
    results = validate_introners(
        putative_introners,
        consensus_sequences,
        gtf_db,
        args.similarity_cutoff,
        args.coverage_cutoff,
        args.flanking_length
    )
    
    # Write log file with validation details
    write_log_file(results, args.output_log, args.similarity_cutoff, args.coverage_cutoff)
    
    # Write filtered FASTA file
    write_filtered_fasta(results, args.output_fasta)
    
    logging.info("Introner validation complete!")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Resolve paralogous protein mappings by selecting the best nucleotide alignment
for each reference protein from multiple miniprot alignments.

This script:
1. Identifies proteins with multiple miniprot alignments
2. Extracts nucleotide sequences for reference and candidate loci
3. Performs nucleotide-level alignment scoring
4. Resolves conflicts to assign one locus per protein
5. Outputs a filtered GTF with redundant alignments removed

Author: Claude Code
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import subprocess
import tempfile
import os

# Data structures
ProteinAlignment = namedtuple('ProteinAlignment', [
    'protein_id', 'scaffold', 'start', 'end', 'strand',
    'alignment_score', 'match_score', 'transcript_id', 'gene_id'
])

GTFRecord = namedtuple('GTFRecord', [
    'seqname', 'source', 'feature', 'start', 'end',
    'score', 'strand', 'frame', 'attributes'
])


class ParalogResolver:
    """Main class for resolving paralogous protein mappings"""

    def __init__(self, ref_gtf, ref_fasta, miniprot_gtf, target_fasta, output_gtf):
        self.ref_gtf = Path(ref_gtf)
        self.ref_fasta = Path(ref_fasta)
        self.miniprot_gtf = Path(miniprot_gtf)
        self.target_fasta = Path(target_fasta)
        self.output_gtf = Path(output_gtf)

        # Data containers
        self.ref_sequences = {}  # transcript_id -> SeqRecord
        self.target_sequences = {}  # scaffold_id -> SeqRecord
        self.protein_alignments = defaultdict(list)  # protein_id -> [ProteinAlignment]
        self.ref_transcript_coords = {}  # transcript_id -> [(start, end, strand, scaffold)]
        self.miniprot_records = []  # All GTF records from miniprot
        self.assigned_loci = set()  # Track assigned genomic loci to prevent double assignment

    def load_sequences(self):
        """Load reference and target FASTA sequences"""
        print("Loading reference sequences...")
        self.ref_sequences = {record.id: record for record in SeqIO.parse(self.ref_fasta, "fasta")}

        print("Loading target sequences...")
        self.target_sequences = {record.id: record for record in SeqIO.parse(self.target_fasta, "fasta")}

        print(f"Loaded {len(self.ref_sequences)} reference sequences")
        print(f"Loaded {len(self.target_sequences)} target sequences")

    def parse_gtf_attributes(self, attr_string):
        """Parse GTF attributes string into dictionary"""
        attrs = {}
        for item in attr_string.strip(';').split(';'):
            if item.strip():
                parts = item.strip().split(' ', 1)
                if len(parts) == 2:
                    key, value = parts
                    attrs[key] = value.strip('"')
                else:
                    # Handle malformed attributes without a value
                    print(f"Warning: Skipping malformed attribute: '{item.strip()}'")
        return attrs

    def build_gtf_attributes(self, attrs_dict):
        """Build GTF attributes string from dictionary"""
        attr_parts = []
        for key, value in attrs_dict.items():
            attr_parts.append(f'{key} "{value}"')
        return '; '.join(attr_parts) + ';'

    def parse_ref_gtf(self):
        """Parse reference GTF to get transcript coordinates"""
        print("Parsing reference GTF...")

        with open(self.ref_gtf) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) != 9:
                    continue

                seqname, source, feature, start, end, score, strand, frame, attributes = fields

                if feature in ['exon', 'CDS']:
                    attrs = self.parse_gtf_attributes(attributes)
                    transcript_id = attrs.get('transcript_id')

                    if transcript_id:
                        if transcript_id not in self.ref_transcript_coords:
                            self.ref_transcript_coords[transcript_id] = []

                        self.ref_transcript_coords[transcript_id].append({
                            'start': int(start),
                            'end': int(end),
                            'strand': strand,
                            'scaffold': seqname,
                            'feature': feature
                        })

        print(f"Parsed coordinates for {len(self.ref_transcript_coords)} reference transcripts")

    def parse_miniprot_gtf(self):
        """Parse miniprot GTF to identify multiple alignments"""
        print("Parsing miniprot GTF...")

        current_protein = None
        current_alignment = {}

        with open(self.miniprot_gtf) as f:
            for line in f:
                line = line.strip()

                # Store all GTF records for later output
                if not line.startswith('#'):
                    fields = line.split('\t')
                    if len(fields) == 9:
                        record = GTFRecord(*fields)
                        self.miniprot_records.append(record)

                # Parse PAF alignment information
                if line.startswith('##PAF'):
                    fields = line.split('\t')
                    protein_id = fields[1]
                    target_scaffold = fields[6]
                    target_start = int(fields[8])  # Target start (0-based)
                    target_end = int(fields[9])    # Target end

                    # Parse alignment score from AS:i: field
                    alignment_score = 0
                    match_score = 0
                    for field in fields[13:]:  # Optional fields start after field 12
                        if field.startswith('AS:i:'):
                            alignment_score = int(field.split(':')[2])
                        elif field.startswith('ms:i:'):
                            match_score = int(field.split(':')[2])

                    # Extract additional info from later lines
                    current_protein = protein_id
                    current_alignment = {
                        'protein_id': protein_id,
                        'scaffold': target_scaffold,
                        'start': target_start,
                        'end': target_end,
                        'alignment_score': alignment_score,
                        'match_score': match_score
                    }

                # Associate with transcript/gene IDs from GTF records
                elif not line.startswith('#') and current_protein:
                    fields = line.split('\t')
                    if len(fields) == 9 and fields[2] == 'transcript':
                        attrs = self.parse_gtf_attributes(fields[8])
                        transcript_id = attrs.get('transcript_id')
                        gene_id = attrs.get('gene_id')
                        strand = fields[6]

                        if transcript_id and current_alignment:
                            alignment = ProteinAlignment(
                                protein_id=current_alignment['protein_id'],
                                scaffold=current_alignment['scaffold'],
                                start=current_alignment['start'],
                                end=current_alignment['end'],
                                strand=strand,
                                alignment_score=current_alignment['alignment_score'],
                                match_score=current_alignment['match_score'],
                                transcript_id=transcript_id,
                                gene_id=gene_id
                            )

                            self.protein_alignments[current_protein].append(alignment)
                            current_alignment = {}

        # Filter to proteins with multiple alignments
        multi_alignment_proteins = {k: v for k, v in self.protein_alignments.items() if len(v) > 1}

        print(f"Found {len(self.protein_alignments)} total protein alignments")
        print(f"Found {len(multi_alignment_proteins)} proteins with multiple alignments")

        # Print some statistics
        alignment_counts = [len(v) for v in multi_alignment_proteins.values()]
        if alignment_counts:
            print(f"Max alignments per protein: {max(alignment_counts)}")
            print(f"Average alignments per multi-alignment protein: {sum(alignment_counts)/len(alignment_counts):.2f}")

        return multi_alignment_proteins

    def extract_reference_sequence(self, protein_id):
        """Extract reference nucleotide sequence for a protein"""
        if protein_id not in self.ref_transcript_coords:
            print(f"Warning: No reference coordinates found for protein {protein_id}")
            return None

        coords = self.ref_transcript_coords[protein_id]
        if not coords:
            return None

        # Group by scaffold and strand
        scaffold = coords[0]['scaffold']
        strand = coords[0]['strand']

        # Sort coordinates by start position
        coords = sorted(coords, key=lambda x: x['start'])

        # Extract sequences for each exon/CDS
        sequences = []
        for coord in coords:
            if coord['feature'] == 'CDS':  # Use CDS regions only
                if scaffold not in self.ref_sequences:
                    print(f"Warning: Scaffold {scaffold} not found in reference FASTA")
                    return None

                seq = self.ref_sequences[scaffold].seq[coord['start']-1:coord['end']]
                sequences.append(seq)

        if not sequences:
            return None

        # Concatenate sequences
        full_sequence = ''.join([str(s) for s in sequences])

        # Reverse complement if on negative strand
        if strand == '-':
            full_sequence = str(Seq(full_sequence).reverse_complement())

        return full_sequence

    def extract_candidate_sequence(self, alignment):
        """Extract candidate loci sequence from target assembly"""
        scaffold = alignment.scaffold
        start = alignment.start
        end = alignment.end
        strand = alignment.strand

        if scaffold not in self.target_sequences:
            print(f"Warning: Target scaffold {scaffold} not found")
            return None

        # Extract sequence
        sequence = self.target_sequences[scaffold].seq[start:end]

        # Reverse complement if on negative strand
        if strand == '-':
            sequence = sequence.reverse_complement()

        return str(sequence)

    def align_sequences(self, ref_seq, candidate_seq, protein_id, candidate_id):
        """Perform nucleotide alignment and return score"""
        if not ref_seq or not candidate_seq:
            return 0

        # Write sequences to temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as ref_file:
            ref_file.write(f">ref_{protein_id}\n{ref_seq}\n")
            ref_path = ref_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as cand_file:
            cand_file.write(f">cand_{candidate_id}\n{candidate_seq}\n")
            cand_path = cand_file.name

        try:
            # Use minimap2 for nucleotide alignment
            cmd = [
                'minimap2', '-a', '-x', 'splice',
                '--secondary=no', '--splice-flank=no',
                cand_path, ref_path
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                print(f"Warning: minimap2 failed for {protein_id} vs {candidate_id}")
                return 0

            # Parse SAM output to get alignment score
            score = 0
            for line in result.stdout.split('\n'):
                if line.startswith('ref_'):
                    fields = line.split('\t')
                    if len(fields) > 10:
                        # Look for alignment score in optional fields
                        for field in fields[11:]:
                            if field.startswith('AS:i:'):
                                score = int(field.split(':')[2])
                                break

            return score

        except subprocess.TimeoutExpired:
            print(f"Warning: minimap2 timeout for {protein_id} vs {candidate_id}")
            return 0
        except Exception as e:
            print(f"Warning: minimap2 error for {protein_id} vs {candidate_id}: {e}")
            return 0
        finally:
            # Clean up temporary files
            try:
                os.unlink(ref_path)
                os.unlink(cand_path)
            except:
                pass

    def resolve_paralogs(self, multi_alignment_proteins):
        """Resolve paralogous mappings using global locus assignment tracking"""
        print("Resolving paralogous mappings with global locus tracking...")

        # Global registry: (scaffold, start, end) -> (protein_id, score, alignment)
        global_locus_assignments = {}

        # Phase 1: Register all single-alignment proteins first
        print("Phase 1: Registering single-alignment proteins...")
        single_count = 0
        for protein_id, alignments in self.protein_alignments.items():
            if len(alignments) == 1:
                alignment = alignments[0]
                locus_key = (alignment.scaffold, alignment.start, alignment.end)

                # Use miniprot alignment score as baseline for single-alignment proteins
                baseline_score = alignment.alignment_score

                # Check for duplicates among single-alignment proteins
                if locus_key in global_locus_assignments:
                    existing_protein, existing_score, _ = global_locus_assignments[locus_key]
                    print(f"  Conflict: {protein_id} vs {existing_protein} at {locus_key}")
                    print(f"    {protein_id} score: {baseline_score}")
                    print(f"    {existing_protein} score: {existing_score}")

                    if baseline_score > existing_score:
                        print(f"    Replacing {existing_protein} with {protein_id}")
                        global_locus_assignments[locus_key] = (protein_id, baseline_score, alignment)
                    else:
                        print(f"    Keeping {existing_protein}")
                else:
                    global_locus_assignments[locus_key] = (protein_id, baseline_score, alignment)
                    single_count += 1

        print(f"Registered {single_count} single-alignment proteins")

        # Phase 2: Process multi-alignment proteins with global conflict checking
        print("Phase 2: Processing multi-alignment proteins...")
        multi_resolved = {}

        for protein_id, alignments in multi_alignment_proteins.items():
            print(f"Processing protein {protein_id} with {len(alignments)} alignments...")

            # Extract reference sequence
            ref_seq = self.extract_reference_sequence(protein_id)
            if not ref_seq:
                print(f"Skipping {protein_id}: no reference sequence")
                continue

            best_score = -1
            best_alignment = None
            candidate_locus = None

            # Score each candidate alignment
            for i, alignment in enumerate(alignments):
                candidate_seq = self.extract_candidate_sequence(alignment)
                if not candidate_seq:
                    continue

                # Perform nucleotide alignment
                nt_score = self.align_sequences(ref_seq, candidate_seq, protein_id, f"{alignment.gene_id}_{i}")
                locus_key = (alignment.scaffold, alignment.start, alignment.end)

                print(f"  Alignment {i+1}: {alignment.scaffold}:{alignment.start}-{alignment.end} nt_score={nt_score}")

                # Check against global registry
                available = True
                if locus_key in global_locus_assignments:
                    existing_protein, existing_score, _ = global_locus_assignments[locus_key]
                    print(f"    Conflict with {existing_protein} (score={existing_score})")

                    if nt_score > existing_score:
                        print(f"    Can displace {existing_protein}")
                        # This locus is available for displacement
                    else:
                        print(f"    Cannot displace {existing_protein}")
                        available = False

                # Only consider this alignment if locus is available or displaceable
                if available and nt_score > best_score:
                    best_score = nt_score
                    best_alignment = alignment
                    candidate_locus = locus_key

            # Assign best available locus to current protein
            if best_alignment and candidate_locus:
                # Handle displacement if necessary
                if candidate_locus in global_locus_assignments:
                    displaced_protein, displaced_score, _ = global_locus_assignments[candidate_locus]
                    print(f"  Displacing {displaced_protein} (score={displaced_score}) with {protein_id} (score={best_score})")

                    # Remove displaced protein from resolved assignments if it was multi-alignment
                    if displaced_protein in multi_resolved:
                        del multi_resolved[displaced_protein]
                else:
                    print(f"  Assigning {protein_id} to available locus {candidate_locus} (score={best_score})")

                # Register the assignment globally
                global_locus_assignments[candidate_locus] = (protein_id, best_score, best_alignment)
                multi_resolved[protein_id] = best_alignment
            else:
                print(f"  No suitable locus found for {protein_id}")

        print(f"Resolved {len(multi_resolved)} multi-alignment proteins")

        # Phase 3: Build final assignments from global registry
        print("Phase 3: Building final assignments...")
        final_assignments = {}

        for (scaffold, start, end), (protein_id, score, alignment) in global_locus_assignments.items():
            final_assignments[protein_id] = alignment

        print(f"Final assignments: {len(final_assignments)} proteins to unique loci")
        print(f"Global loci count: {len(global_locus_assignments)}")

        # Validation: Check for any remaining duplicates
        locus_keys = list(global_locus_assignments.keys())
        if len(locus_keys) != len(set(locus_keys)):
            print("ERROR: Duplicate loci still exist in global registry!")
        else:
            print("✓ Validation passed: All loci are unique")

        return final_assignments

    def write_filtered_gtf(self, resolved_assignments):
        """Write filtered GTF with protein IDs replacing generic transcript IDs"""
        print("Writing filtered GTF with protein IDs...")

        # Create mapping: transcript_id -> protein_id using only resolved assignments
        # The resolved_assignments now contains ALL proteins (single + multi) with unique loci
        transcript_to_protein = {}

        for protein_id, alignment in resolved_assignments.items():
            transcript_to_protein[alignment.transcript_id] = protein_id

        print(f"Created mapping for {len(transcript_to_protein)} transcripts to protein IDs")

        # Write filtered GTF with protein IDs
        records_written = 0
        with open(self.output_gtf, 'w') as out_file:
            current_protein = None
            keep_current = False

            with open(self.miniprot_gtf) as in_file:
                for line in in_file:
                    line = line.strip()

                    if line.startswith('##PAF'):
                        fields = line.split('\t')
                        protein_id = fields[1]
                        current_protein = protein_id

                        # Check if this protein should be kept and coordinates match
                        if protein_id in resolved_assignments:
                            alignment = resolved_assignments[protein_id]
                            keep_current = (line.split('\t')[8] == str(alignment.start) and
                                          line.split('\t')[9] == str(alignment.end))
                        else:
                            keep_current = False

                        if keep_current:
                            out_file.write(line + '\n')

                    elif line.startswith('##'):
                        # Write other comment lines if we're keeping current protein
                        if keep_current:
                            out_file.write(line + '\n')

                    elif line and not line.startswith('#'):
                        # GTF record
                        fields = line.split('\t')
                        if len(fields) == 9:
                            feature_type = fields[2]
                            attrs = self.parse_gtf_attributes(fields[8])
                            transcript_id = attrs.get('transcript_id')

                            # Skip gene features entirely
                            if feature_type == 'gene':
                                continue

                            # Only process features with transcripts we want to keep
                            if transcript_id in transcript_to_protein:
                                protein_id = transcript_to_protein[transcript_id]

                                # Create new attributes with protein ID
                                new_attrs = {'transcript_id': protein_id}
                                new_attr_string = self.build_gtf_attributes(new_attrs)

                                # Reconstruct the GTF line with new attributes
                                fields[8] = new_attr_string
                                modified_line = '\t'.join(fields)

                                out_file.write(modified_line + '\n')
                                records_written += 1

        print(f"Wrote {records_written} GTF records with protein IDs to {self.output_gtf}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Resolve paralogous protein mappings")
    parser.add_argument("--ref-gtf", required=True, help="Reference GTF file")
    parser.add_argument("--ref-fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--miniprot-gtf", required=True, help="Miniprot GTF file with multiple alignments")
    parser.add_argument("--target-fasta", required=True, help="Target assembly FASTA file")
    parser.add_argument("--output-gtf", required=True, help="Output GTF file with resolved alignments")

    args = parser.parse_args()

    # Check input files exist
    for file_path in [args.ref_gtf, args.ref_fasta, args.miniprot_gtf, args.target_fasta]:
        if not Path(file_path).exists():
            print(f"Error: Input file not found: {file_path}")
            sys.exit(1)

    # Initialize resolver
    resolver = ParalogResolver(
        ref_gtf=args.ref_gtf,
        ref_fasta=args.ref_fasta,
        miniprot_gtf=args.miniprot_gtf,
        target_fasta=args.target_fasta,
        output_gtf=args.output_gtf
    )

    # Run the pipeline
    try:
        resolver.load_sequences()
        resolver.parse_ref_gtf()
        multi_alignment_proteins = resolver.parse_miniprot_gtf()

        print(f"\nSuccessfully parsed input files!")
        print(f"Processing {len(multi_alignment_proteins)} proteins with multiple alignments")

        # Resolve paralogous mappings
        resolved_assignments = resolver.resolve_paralogs(multi_alignment_proteins)

        # Write filtered GTF
        resolver.write_filtered_gtf(resolved_assignments)

        print(f"\nParalog resolution complete!")
        print(f"Output written to: {args.output_gtf}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

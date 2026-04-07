#!/usr/bin/env python3
"""
Introner Insertion Pattern Analysis - Phase 1
Data Integration & Basic Patterns (GTF-Based)

Creates gene-level classification matrix from ortholog-level introner data.
Uses GTF annotations for accurate gene identification instead of genotype matrix.
Classifies introners as fixed vs polymorphic within group1 and group2.

Usage:
python scripts/introner_phase1_analysis.py genotype_matrixes/genotype_matrix_updated.tsv gene_annotations/ data/gene2go.json results/phase1_gene_classification.tsv
"""

import pandas as pd
import json
import sys
import os
import re
from collections import defaultdict, Counter

def load_introner_data(genotype_file):
    """Load and clean introner genotype matrix."""
    print(f"Loading introner data from {genotype_file}...")
    
    # Read TSV file, skip duplicate header on line 2
    df = pd.read_csv(genotype_file, sep='\t', header=0)
    df = df[df['ortholog_id'] != 'ortholog_id']  # Remove duplicate header row
    
    print(f"  Loaded {len(df)} records")
    print(f"  Unique ortholog groups: {df['ortholog_id'].nunique()}")
    print(f"  Unique samples: {df['sample'].nunique()}")
    print(f"  Unique genes: {df['gene'].nunique()}")
    
    return df

def load_go_annotations(go_file):
    """Load GO annotations."""
    print(f"Loading GO annotations from {go_file}...")

    with open(go_file, 'r') as f:
        go_data = json.load(f)

    print(f"  Loaded GO annotations for {len(go_data)} genes")
    return go_data

def load_gtf_annotations(gtf_dir):
    """
    Load introner annotations from all sample GTF files and build comprehensive gene-introner database.

    Returns:
        dict: {sample: {gene_id: {ortholog_id: {'introner_id': str, 'family': str,
                                               'contig': str, 'start': int, 'end': int}}}}
    """
    print(f"Loading GTF annotations from {gtf_dir}...")

    gene_introner_db = {}
    gtf_files = [f for f in os.listdir(gtf_dir) if f.endswith('.plus_introners.gtf')]

    for gtf_file in gtf_files:
        sample = gtf_file.replace('.plus_introners.gtf', '')
        gene_introner_db[sample] = defaultdict(dict)

        gtf_path = os.path.join(gtf_dir, gtf_file)
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                if feature_type != 'introner':
                    continue

                contig = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]

                # Parse attributes
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                family_match = re.search(r'family "([^"]+)"', attributes)
                id_match = re.search(r'ID "([^"]+)"', attributes)
                ortholog_id_match = re.search(r'ortholog_id "([^"]+)"', attributes)

                if gene_id_match and family_match and id_match and ortholog_id_match:
                    gene_id = gene_id_match.group(1)
                    family = family_match.group(1)
                    introner_id = id_match.group(1)
                    ortholog_id = ortholog_id_match.group(1)

                    # Clean gene ID for consistency
                    gene_clean = clean_gene_id(gene_id)
                    if gene_clean:
                        # Build nested structure: sample -> gene -> ortholog
                        gene_introner_db[sample][gene_clean][ortholog_id] = {
                            'introner_id': introner_id,
                            'family': family,
                            'contig': contig,
                            'start': start,
                            'end': end
                        }

    # Count total introners and genes
    total_introners = 0
    total_genes = set()
    for sample_data in gene_introner_db.values():
        for gene_id, orthologs in sample_data.items():
            total_introners += len(orthologs)
            total_genes.add(gene_id)

    print(f"  Loaded {total_introners} introner annotations from {len(gtf_files)} GTF files")
    print(f"  Found {len(total_genes)} unique genes with introners")
    print(f"  Samples: {sorted(gene_introner_db.keys())}")

    return gene_introner_db

def classify_samples_into_groups(df):
    """Classify samples into Group1 and Group2 based on CLAUDE.md specification."""
    print("Classifying samples into groups...")

    # Group2 = {RCC1749, RCC3052}, Group1 = all others
    group2_samples = {'RCC1749', 'RCC3052'}
    all_samples = set(df['sample'].unique())
    group1_samples = all_samples - group2_samples

    print(f"  Group1 samples ({len(group1_samples)}): {sorted(group1_samples)}")
    print(f"  Group2 samples ({len(group2_samples)}): {sorted(group2_samples)}")

    return group1_samples, group2_samples

def create_master_gene_ortholog_matrix(gene_introner_db, df, group1_samples, group2_samples):
    """
    Create comprehensive gene-ortholog presence/absence matrix using GTF and genotype data.

    Args:
        gene_introner_db: GTF-derived gene-introner database
        df: Genotype matrix DataFrame
        group1_samples, group2_samples: Sample group classifications

    Returns:
        dict: {gene_id: {'all_orthologs': set, 'sample_presence': {sample: {ortholog: presence}}}}
    """
    print("Creating master gene-ortholog presence/absence matrix...")

    # First, collect all gene-ortholog relationships from GTF
    master_gene_orthologs = defaultdict(lambda: {
        'all_orthologs': set(),
        'sample_presence': defaultdict(dict)
    })

    # Initialize from GTF data - these are the present introners
    for sample, sample_genes in gene_introner_db.items():
        for gene_id, orthologs in sample_genes.items():
            master_gene_orthologs[gene_id]['all_orthologs'].update(orthologs.keys())
            for ortholog_id in orthologs:
                master_gene_orthologs[gene_id]['sample_presence'][sample][ortholog_id] = 1  # Present

    # Now fill in presence/absence from genotype matrix
    genotype_records = 0
    for _, row in df.iterrows():
        ortholog_id = row['ortholog_id']
        sample = row['sample']
        presence = int(row['presence'])

        # Find which gene this ortholog belongs to from GTF data
        target_gene = None
        for gene_id, gene_data in master_gene_orthologs.items():
            if ortholog_id in gene_data['all_orthologs']:
                target_gene = gene_id
                break

        if target_gene:
            # Update presence status from genotype matrix (may override GTF presence)
            master_gene_orthologs[target_gene]['sample_presence'][sample][ortholog_id] = presence
            genotype_records += 1

    # Ensure all samples have entries for all orthologs in each gene
    all_samples = group1_samples | group2_samples
    for gene_id, gene_data in master_gene_orthologs.items():
        for sample in all_samples:
            for ortholog_id in gene_data['all_orthologs']:
                if ortholog_id not in gene_data['sample_presence'][sample]:
                    # If not in genotype matrix and not in GTF for this sample, assume absent
                    gene_data['sample_presence'][sample][ortholog_id] = 2  # Absent

    total_genes = len(master_gene_orthologs)
    total_relationships = sum(len(data['all_orthologs']) for data in master_gene_orthologs.values())

    print(f"  Created matrix for {total_genes} genes with {total_relationships} gene-ortholog relationships")
    print(f"  Integrated {genotype_records} genotype records")

    return dict(master_gene_orthologs)

def clean_gene_id(gene_id):
    """Clean gene ID by removing version suffix (.3.0.228) to match GO file format."""
    if pd.isna(gene_id) or not str(gene_id).strip():
        return None
    gene_id = str(gene_id).strip()
    # Remove version suffix like .3.0.228
    if '.3.0.228' in gene_id:
        gene_id = gene_id.replace('.3.0.228', '')
    return gene_id

def infer_gene_names_for_ortholog_groups(df):
    """
    Infer gene names for records without gene names using other records 
    in the same ortholog group. Absent/missing introners don't have gene names,
    but we can infer them from present introners in the same ortholog group.
    """
    print("  Inferring gene names for absent/missing introners...")
    
    # Create mapping of ortholog_id -> gene_name from records that have gene names
    ortholog_to_gene = {}
    for _, row in df.iterrows():
        if pd.notna(row['gene']) and str(row['gene']).strip():
            ortholog_id = row['ortholog_id']
            gene_name = str(row['gene']).strip()
            if ortholog_id not in ortholog_to_gene:
                ortholog_to_gene[ortholog_id] = gene_name
    
    # Apply inferred gene names to records without gene names
    df_copy = df.copy()
    inferred_count = 0
    for idx, row in df_copy.iterrows():
        if pd.isna(row['gene']) or not str(row['gene']).strip():
            ortholog_id = row['ortholog_id']
            if ortholog_id in ortholog_to_gene:
                df_copy.at[idx, 'gene'] = ortholog_to_gene[ortholog_id]
                inferred_count += 1
    
    print(f"    Inferred gene names for {inferred_count} records")
    return df_copy

def classify_introner_patterns_from_matrix(master_gene_orthologs, group1_samples, group2_samples):
    """
    Classify genes into fixed vs polymorphic introner categories using comprehensive presence/absence matrix.

    Args:
        master_gene_orthologs: Comprehensive gene-ortholog presence/absence matrix
        group1_samples, group2_samples: Sample group classifications

    Returns:
        list: Gene classification results with fixed/polymorphic counts
    """
    print("Classifying introner patterns from comprehensive matrix...")

    results = []

    for gene_id, gene_data in master_gene_orthologs.items():
        group1_fixed = 0
        group1_polymorphic = 0
        group2_fixed = 0
        group2_polymorphic = 0

        # Analyze each ortholog in this gene
        for ortholog_id in gene_data['all_orthologs']:
            # Get presence status for this ortholog across all samples
            group1_presences = []
            group2_presences = []

            for sample in group1_samples:
                if sample in gene_data['sample_presence'] and ortholog_id in gene_data['sample_presence'][sample]:
                    presence = gene_data['sample_presence'][sample][ortholog_id]
                    if presence in [1, 2]:  # Only count present/absent, not missing
                        group1_presences.append(presence)

            for sample in group2_samples:
                if sample in gene_data['sample_presence'] and ortholog_id in gene_data['sample_presence'][sample]:
                    presence = gene_data['sample_presence'][sample][ortholog_id]
                    if presence in [1, 2]:  # Only count present/absent, not missing
                        group2_presences.append(presence)

            # Classify this ortholog for Group1
            if len(group1_presences) > 0:
                if all(p == 1 for p in group1_presences):
                    # Present in all Group1 samples with data
                    if len(group1_presences) == len(group1_samples):
                        group1_fixed += 1  # Fixed (present in all)
                    else:
                        group1_polymorphic += 1  # Polymorphic (present where data exists, missing elsewhere)
                elif all(p == 2 for p in group1_presences):
                    # Absent in all Group1 samples with data - don't count
                    pass
                else:
                    # Mix of present and absent
                    group1_polymorphic += 1

            # Classify this ortholog for Group2
            if len(group2_presences) > 0:
                if all(p == 1 for p in group2_presences):
                    # Present in all Group2 samples with data
                    if len(group2_presences) == len(group2_samples):
                        group2_fixed += 1  # Fixed (present in all)
                    else:
                        group2_polymorphic += 1  # Polymorphic (present where data exists, missing elsewhere)
                elif all(p == 2 for p in group2_presences):
                    # Absent in all Group2 samples with data - don't count
                    pass
                else:
                    # Mix of present and absent
                    group2_polymorphic += 1

        total_orthologs = len(gene_data['all_orthologs'])

        results.append({
            'gene_id': gene_id,
            'group1_fixed': group1_fixed,
            'group1_polymorphic': group1_polymorphic,
            'group2_fixed': group2_fixed,
            'group2_polymorphic': group2_polymorphic,
            'total_orthologs': total_orthologs
        })

    print(f"  Classified {len(results)} genes")
    return results

def classify_introner_patterns(gene_data, group1_samples, group2_samples):
    """
    Classify genes into fixed vs polymorphic introner categories.
    
    NOTE: Ortholog groups with only missing data (presence=3) are excluded 
    from all counts including total_orthologs. Only ortholog groups with 
    at least some presence/absence data (presence=1,2) are counted.
    """
    print("Classifying introner patterns...")
    
    results = []
    
    for gene, data in gene_data.items():
        # Count orthologs with presence=1 (introner present) by group
        group1_present_orthologs = set()
        group1_total_orthologs = set()
        group2_present_orthologs = set()
        group2_total_orthologs = set()
        
        # Group1 analysis
        group1_by_ortholog = defaultdict(list)
        for record in data['group1_records']:
            ortholog_id = record['ortholog_id']
            presence = record['presence']
            group1_by_ortholog[ortholog_id].append(presence)
            # Only count orthologs with actual presence/absence data (not missing data)
            if presence in [1, 2]:
                group1_total_orthologs.add(ortholog_id)
            if presence == 1:
                group1_present_orthologs.add(ortholog_id)
        
        # Group2 analysis
        group2_by_ortholog = defaultdict(list)
        for record in data['group2_records']:
            ortholog_id = record['ortholog_id']
            presence = record['presence']
            group2_by_ortholog[ortholog_id].append(presence)
            # Only count orthologs with actual presence/absence data (not missing data)
            if presence in [1, 2]:
                group2_total_orthologs.add(ortholog_id)
            if presence == 1:
                group2_present_orthologs.add(ortholog_id)
        
        # Classify fixed vs polymorphic for each group
        group1_fixed = 0
        group1_polymorphic = 0
        group2_fixed = 0
        group2_polymorphic = 0
        
        # Group1 classification
        for ortholog_id in group1_total_orthologs:
            presences = group1_by_ortholog[ortholog_id]
            if len(presences) > 0:
                # Fixed: present in all samples (all presence = 1)
                # Polymorphic: present in some samples (mix of 1 and 2)
                valid_presences = [p for p in presences if p in [1, 2]]
                if len(valid_presences) > 0 and all(p == 1 for p in valid_presences):
                    group1_fixed += 1
                elif 1 in valid_presences and 2 in valid_presences:
                    group1_polymorphic += 1
        
        # Group2 classification
        for ortholog_id in group2_total_orthologs:
            presences = group2_by_ortholog[ortholog_id]
            if len(presences) > 0:
                valid_presences = [p for p in presences if p in [1, 2]]
                if len(valid_presences) > 0 and all(p == 1 for p in valid_presences):
                    group2_fixed += 1
                elif 1 in valid_presences and 2 in valid_presences:
                    group2_polymorphic += 1
        
        # Calculate total orthologs (excluding missing-data-only groups)
        all_classified_orthologs = group1_total_orthologs | group2_total_orthologs
        
        results.append({
            'gene_id': gene,
            'group1_fixed': group1_fixed,
            'group1_polymorphic': group1_polymorphic,
            'group2_fixed': group2_fixed,
            'group2_polymorphic': group2_polymorphic,
            'total_orthologs': len(all_classified_orthologs)
        })
    
    return results

def add_go_annotations(results, go_data):
    """Add GO annotations to gene classification results."""
    print("Adding GO annotations...")
    
    go_matched = 0
    for result in results:
        gene_id = result['gene_id']
        if gene_id in go_data:
            result['go_terms'] = go_data[gene_id]
            go_matched += 1
        else:
            result['go_terms'] = []
    
    print(f"  Matched GO annotations for {go_matched}/{len(results)} genes")
    return results

def generate_statistics(results):
    """Generate descriptive statistics."""
    print("\n=== PHASE 1 DESCRIPTIVE STATISTICS ===")
    
    total_genes = len(results)
    print(f"Total genes analyzed: {total_genes}")
    
    # Genes with any introners
    genes_with_introners = [r for r in results if 
                           r['group1_fixed'] > 0 or r['group1_polymorphic'] > 0 or 
                           r['group2_fixed'] > 0 or r['group2_polymorphic'] > 0]
    
    print(f"Genes with introners: {len(genes_with_introners)} ({len(genes_with_introners)/total_genes*100:.1f}%)")
    
    # Group-specific statistics
    group1_genes = [r for r in results if r['group1_fixed'] > 0 or r['group1_polymorphic'] > 0]
    group2_genes = [r for r in results if r['group2_fixed'] > 0 or r['group2_polymorphic'] > 0]
    both_groups = [r for r in results if 
                   (r['group1_fixed'] > 0 or r['group1_polymorphic'] > 0) and
                   (r['group2_fixed'] > 0 or r['group2_polymorphic'] > 0)]
    
    print(f"Genes with Group1 introners: {len(group1_genes)}")
    print(f"Genes with Group2 introners: {len(group2_genes)}")
    print(f"Genes with both Group1 and Group2 introners: {len(both_groups)}")
    
    # Fixed vs polymorphic counts
    total_group1_fixed = sum(r['group1_fixed'] for r in results)
    total_group1_polymorphic = sum(r['group1_polymorphic'] for r in results)
    total_group2_fixed = sum(r['group2_fixed'] for r in results)
    total_group2_polymorphic = sum(r['group2_polymorphic'] for r in results)
    
    print(f"\nIntroner ortholog counts:")
    print(f"  Group1 fixed: {total_group1_fixed}")
    print(f"  Group1 polymorphic: {total_group1_polymorphic}")
    print(f"  Group2 fixed: {total_group2_fixed}")
    print(f"  Group2 polymorphic: {total_group2_polymorphic}")
    
    # GO annotation coverage
    genes_with_go = [r for r in results if len(r['go_terms']) > 0]
    print(f"\nGO annotation coverage: {len(genes_with_go)}/{total_genes} ({len(genes_with_go)/total_genes*100:.1f}%)")

def save_results(results, output_file):
    """Save gene classification matrix to file."""
    print(f"Saving results to {output_file}...")
    
    # Convert to DataFrame for easy saving
    df_results = pd.DataFrame(results)
    
    # Convert GO terms list to string for TSV output
    df_results['go_terms_str'] = df_results['go_terms'].apply(
        lambda x: ';'.join(x) if x else ''
    )
    
    # Reorder columns
    columns = [
        'gene_id', 'group1_fixed', 'group1_polymorphic', 
        'group2_fixed', 'group2_polymorphic', 'total_orthologs', 'go_terms_str'
    ]
    df_results = df_results[columns]
    
    # Save as TSV
    df_results.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved {len(df_results)} gene records")

def main():
    if len(sys.argv) != 5:
        print("Usage: python introner_phase1_analysis.py <genotype_matrix.tsv> <gtf_dir> <gene2go.json> <output.tsv>")
        print("\nExample:")
        print("python introner_phase1_analysis.py \\")
        print("  genotype_matrixes/genotype_matrix_updated.tsv \\")
        print("  gene_annotations/ \\")
        print("  data/gene2go.json \\")
        print("  results/phase1_gene_classification.tsv")
        sys.exit(1)

    genotype_file = sys.argv[1]
    gtf_dir = sys.argv[2]
    go_file = sys.argv[3]
    output_file = sys.argv[4]

    print("=== INTRONER INSERTION PATTERN ANALYSIS - PHASE 1 (GTF-BASED) ===")
    print(f"Introner genotype matrix: {genotype_file}")
    print(f"GTF annotations directory: {gtf_dir}")
    print(f"GO annotations: {go_file}")
    print(f"Output file: {output_file}")
    print()

    try:
        # Load data
        df = load_introner_data(genotype_file)
        gene_introner_db = load_gtf_annotations(gtf_dir)
        go_data = load_go_annotations(go_file)

        # Classify samples into groups
        group1_samples, group2_samples = classify_samples_into_groups(df)

        # Create comprehensive gene-ortholog presence/absence matrix
        master_gene_orthologs = create_master_gene_ortholog_matrix(gene_introner_db, df, group1_samples, group2_samples)

        # Classify introner patterns using the comprehensive matrix
        results = classify_introner_patterns_from_matrix(master_gene_orthologs, group1_samples, group2_samples)

        # Add GO annotations
        results = add_go_annotations(results, go_data)

        # Generate statistics
        generate_statistics(results)

        # Save results
        save_results(results, output_file)

        print("\n✓ Phase 1 analysis completed successfully!")

    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
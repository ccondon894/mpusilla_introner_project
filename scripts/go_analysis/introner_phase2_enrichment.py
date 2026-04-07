#!/usr/bin/env python3
"""
Introner Insertion Pattern Analysis - Phase 2
Functional Enrichment Analysis

Performs GO enrichment analysis for different introner categories:
- Group1-Fixed, Group1-Polymorphic, Group2-Fixed, Group2-Polymorphic

Uses GO Slim categories for better statistical power and Fisher's exact test.
python scripts/introner_phase2_enrichment.py results/phase1_gene_classification.tsv results/phase2_enrichment_results.tsv
"""

import argparse
import pandas as pd
import json
import sys
import numpy as np
from collections import defaultdict, Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def load_phase1_results(phase1_file, go_json):
    """
    Load Phase 1 gene classification results and merge with full GO annotation data.

    CRITICAL FIX: Use all GO-annotated genes as background, not just introner-matrix genes.
    """
    print(f"Loading Phase 1 results from {phase1_file}...")

    # Load Phase 1 introner classifications
    phase1_df = pd.read_csv(phase1_file, sep='\t')
    print(f"  Loaded {len(phase1_df)} genes from introner analysis")

    # Parse GO terms back from string format
    def parse_go_terms(x):
        if pd.isna(x) or not str(x).strip():
            return []
        return str(x).strip().split(';')

    phase1_df['go_terms'] = phase1_df['go_terms_str'].apply(parse_go_terms)

    # Load FULL GO annotation data (this is the proper background!)
    go_file = go_json
    print(f"  Loading full GO annotation data from {go_file}...")
    
    try:
        with open(go_file, 'r') as f:
            full_go_data = json.load(f)
        print(f"  Loaded GO annotations for {len(full_go_data)} total genes")
    except FileNotFoundError:
        print(f"  ERROR: Could not find {go_file}")
        print(f"  Falling back to Phase 1 genes only (incorrect background)")
        df_with_go = phase1_df[phase1_df['go_terms'].apply(len) > 0]
        print(f"  {len(df_with_go)} genes have GO annotations")
        return df_with_go
    
    # Create comprehensive dataset with all GO-annotated genes
    all_genes_data = []
    
    for gene_id, go_terms in full_go_data.items():
        # Check if this gene appears in Phase 1 introner data
        phase1_match = phase1_df[phase1_df['gene_id'] == gene_id]
        
        if len(phase1_match) > 0:
            # Gene has introner data - use Phase 1 classifications
            row = phase1_match.iloc[0].to_dict()
            row['go_terms'] = go_terms  # Use GO data from gene2go.json
        else:
            # Gene not in introner data - classify as introner-free
            row = {
                'gene_id': gene_id,
                'group1_fixed': 0,
                'group1_polymorphic': 0,
                'group2_fixed': 0,
                'group2_polymorphic': 0,
                'total_orthologs': 0,
                'go_terms_str': ';'.join(go_terms),
                'go_terms': go_terms
            }
        
        all_genes_data.append(row)
    
    # Create comprehensive dataframe
    df_complete = pd.DataFrame(all_genes_data)
    
    # Add introner classification flags
    df_complete['has_group1'] = (df_complete['group1_fixed'] > 0) | (df_complete['group1_polymorphic'] > 0)
    df_complete['has_group2'] = (df_complete['group2_fixed'] > 0) | (df_complete['group2_polymorphic'] > 0)
    df_complete['has_any_introner'] = df_complete['has_group1'] | df_complete['has_group2']
    
    # Filter genes with GO annotations
    df_with_go = df_complete[df_complete['go_terms'].apply(len) > 0]
    
    print(f"  CORRECTED BACKGROUND: {len(df_with_go)} total genes with GO annotations")
    print(f"    Genes with introners: {df_with_go['has_any_introner'].sum()}")
    print(f"    Genes without introners: {(~df_with_go['has_any_introner']).sum()}")
    print(f"    Group1 genes: {df_with_go['has_group1'].sum()}")
    print(f"    Group2 genes: {df_with_go['has_group2'].sum()}")
    
    return df_with_go

def create_go_slim_mapping():
    """
    Create GO Slim mapping for broad functional categories.
    Maps specific GO terms to broader categories for better statistical power.
    """
    print("Creating GO Slim mapping...")
    
    # Define broad GO Slim categories based on common functional themes
    go_slim_mapping = {
        # Biological Process categories
        'GO:0008152': 'Metabolic Process',  # metabolic process
        'GO:0009058': 'Biosynthetic Process',  # biosynthetic process
        'GO:0009056': 'Catabolic Process',  # catabolic process
        'GO:0006810': 'Transport',  # transport
        'GO:0051234': 'Establishment of Localization',  # establishment of localization
        'GO:0006412': 'Translation',  # translation
        'GO:0006351': 'Transcription DNA-templated',  # transcription, DNA-templated
        'GO:0006355': 'Regulation of Transcription',  # regulation of transcription
        'GO:0006396': 'RNA Processing',  # RNA processing
        'GO:0008219': 'Cell Death',  # cell death
        'GO:0007049': 'Cell Cycle',  # cell cycle
        'GO:0000278': 'Mitotic Cell Cycle',  # mitotic cell cycle
        'GO:0006464': 'Cellular Protein Modification',  # cellular protein modification process
        'GO:0016043': 'Cellular Component Organization',  # cellular component organization
        'GO:0050896': 'Response to Stimulus',  # response to stimulus
        'GO:0009607': 'Response to Biotic Stimulus',  # response to biotic stimulus
        'GO:0009628': 'Response to Abiotic Stimulus',  # response to abiotic stimulus
        'GO:0009719': 'Response to Endogenous Stimulus',  # response to endogenous stimulus
        'GO:0015979': 'Photosynthesis',  # photosynthesis
        'GO:0009765': 'Photosynthesis Light Reactions',  # photosynthesis, light reactions
        'GO:0019684': 'Photosynthesis Light Harvesting',  # photosynthesis, light harvesting
        
        # Molecular Function categories  
        'GO:0003824': 'Catalytic Activity',  # catalytic activity
        'GO:0016787': 'Hydrolase Activity',  # hydrolase activity
        'GO:0016740': 'Transferase Activity',  # transferase activity
        'GO:0016491': 'Oxidoreductase Activity',  # oxidoreductase activity
        'GO:0016874': 'Ligase Activity',  # ligase activity
        'GO:0005515': 'Protein Binding',  # protein binding
        'GO:0003677': 'DNA Binding',  # DNA binding
        'GO:0003723': 'RNA Binding',  # RNA binding
        'GO:0022857': 'Transmembrane Transporter Activity',  # transmembrane transporter activity
        'GO:0005524': 'ATP Binding',  # ATP binding
        'GO:0000166': 'Nucleotide Binding',  # nucleotide binding
        'GO:0046872': 'Metal Ion Binding',  # metal ion binding
        
        # Cellular Component categories
        'GO:0005634': 'Nucleus',  # nucleus
        'GO:0005737': 'Cytoplasm',  # cytoplasm
        'GO:0005739': 'Mitochondrion',  # mitochondrion
        'GO:0009507': 'Chloroplast',  # chloroplast
        'GO:0005783': 'Endoplasmic Reticulum',  # endoplasmic reticulum
        'GO:0005794': 'Golgi Apparatus',  # Golgi apparatus
        'GO:0005886': 'Plasma Membrane',  # plasma membrane
        'GO:0016020': 'Membrane',  # membrane
        'GO:0005622': 'Intracellular',  # intracellular
        'GO:0005840': 'Ribosome',  # ribosome
        'GO:0005856': 'Cytoskeleton',  # cytoskeleton
        'GO:0005773': 'Vacuole',  # vacuole
    }
    
    print(f"  Defined {len(go_slim_mapping)} GO Slim categories")
    return go_slim_mapping

def map_genes_to_go_slim(df, go_slim_mapping):
    """Map genes to GO Slim categories."""
    print("Mapping genes to GO Slim categories...")
    
    gene_to_slim = defaultdict(set)
    slim_to_genes = defaultdict(set)
    
    for _, row in df.iterrows():
        gene_id = row['gene_id']
        go_terms = row['go_terms']
        
        for go_term in go_terms:
            if go_term in go_slim_mapping:
                slim_category = go_slim_mapping[go_term]
                gene_to_slim[gene_id].add(slim_category)
                slim_to_genes[slim_category].add(gene_id)
    
    # Convert to regular dicts with lists
    gene_to_slim = {gene: list(categories) for gene, categories in gene_to_slim.items()}
    slim_to_genes = {category: list(genes) for category, genes in slim_to_genes.items()}
    
    print(f"  Mapped {len(gene_to_slim)} genes to GO Slim categories")
    print(f"  Found {len(slim_to_genes)} GO Slim categories with genes")
    
    # Print category sizes
    for category, genes in sorted(slim_to_genes.items(), key=lambda x: len(x[1]), reverse=True):
        print(f"    {category}: {len(genes)} genes")
    
    return gene_to_slim, slim_to_genes

def create_introner_gene_sets(df):
    """Create overlapping gene sets for simplified introner categories."""
    print("Creating simplified overlapping introner gene sets...")

    # Create simplified overlapping categories
    gene_sets = {}

    # all-introner: any gene containing at least one introner from either group
    gene_sets['all-introners'] = set(df[df['has_any_introner']]['gene_id'])

    # group1-introner: any gene containing at least one Group1 introner (fixed or polymorphic)
    mask = df['has_group1']
    gene_sets['group1'] = set(df[mask]['gene_id'])

    # group2-introner: any gene containing at least one Group2 introner (fixed or polymorphic)
    mask = df['has_group2']
    gene_sets['group2'] = set(df[mask]['gene_id'])

    background_genes = set(df['gene_id'])

    print(f"  Background genes: {len(background_genes)}")
    for category, genes in gene_sets.items():
        print(f"  {category}: {len(genes)} genes")

    # Validation: report overlaps (these are expected and allowed now)
    print("  Overlap analysis (expected for overlapping categories):")
    group1_and_group2 = gene_sets['group1'] & gene_sets['group2']
    print(f"    Genes with both Group1 and Group2 introners: {len(group1_and_group2)}")

    # Verify that all-introner is the union
    union = gene_sets['group1'] | gene_sets['group2']
    if gene_sets['all-introners'] == union:
        print(f"    Validation: all-introner correctly equals union of group1 and group2")
    else:
        print(f"    WARNING: all-introner does not equal union of groups!")
        print(f"      all-introner: {len(gene_sets['all-introners'])}")
        print(f"      union: {len(union)}")

    introner_free_genes = len(df[~df['has_any_introner']])
    print(f"  Introner-free genes: {introner_free_genes} genes")

    return gene_sets, background_genes

def perform_enrichment_analysis(gene_sets, background_genes, slim_to_genes, min_category_size=5):
    """
    Perform Fisher's exact test for GO enrichment analysis.
    """
    print("Performing enrichment analysis...")
    
    results = []
    
    for introner_type, introner_genes in gene_sets.items():
        print(f"  Analyzing {introner_type} ({len(introner_genes)} genes)...")
        
        for go_category, category_genes in slim_to_genes.items():
            category_genes_set = set(category_genes)
            
            # Skip categories with too few genes
            if len(category_genes_set) < min_category_size:
                continue
            
            # Create 2x2 contingency table
            # [introner_with_go, introner_without_go, non_introner_with_go, non_introner_without_go]
            introner_with_go = len(introner_genes & category_genes_set)
            introner_without_go = len(introner_genes - category_genes_set)
            non_introner_with_go = len(category_genes_set - introner_genes)
            non_introner_without_go = len(background_genes - introner_genes - category_genes_set)
            
            # Skip if no overlap
            if introner_with_go == 0:
                continue
            
            # Fisher's exact test
            contingency_table = [[introner_with_go, introner_without_go],
                               [non_introner_with_go, non_introner_without_go]]
            
            odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
            
            # Calculate fold enrichment
            expected = (len(introner_genes) * len(category_genes_set)) / len(background_genes)
            fold_enrichment = introner_with_go / expected if expected > 0 else float('inf')
            
            results.append({
                'introner_type': introner_type,
                'go_category': go_category,
                'genes_with_introner_and_go': introner_with_go,
                'genes_with_introner_total': len(introner_genes),
                'genes_with_go_total': len(category_genes_set),
                'background_total': len(background_genes),
                'odds_ratio': odds_ratio,
                'fold_enrichment': fold_enrichment,
                'p_value': p_value,
                'contingency_table': str(contingency_table)
            })
    
    print(f"  Completed {len(results)} enrichment tests")
    return results

def apply_multiple_testing_correction(results, fdr_threshold=0.05):
    """Apply FDR correction for multiple testing."""
    print("Applying multiple testing correction...")

    if not results:
        return results

    df_results = pd.DataFrame(results)

    # Apply FDR correction within each introner type
    for introner_type in df_results['introner_type'].unique():
        mask = df_results['introner_type'] == introner_type
        p_values = df_results.loc[mask, 'p_value'].values

        if len(p_values) > 0:
            _, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')
            df_results.loc[mask, 'fdr'] = fdr_values

    # Add significance flags
    df_results['significant'] = df_results['fdr'] < fdr_threshold
    df_results['enrichment_type'] = df_results.apply(
        lambda row: 'Enriched' if row['odds_ratio'] > 1 and row['significant']
                   else 'Depleted' if row['odds_ratio'] < 1 and row['significant']
                   else 'Not Significant', axis=1
    )

    return df_results

def generate_summary_statistics(df_results):
    """Generate summary statistics for enrichment results."""
    print("\n=== PHASE 2 ENRICHMENT ANALYSIS SUMMARY ===")
    
    for introner_type in df_results['introner_type'].unique():
        type_results = df_results[df_results['introner_type'] == introner_type]
        
        print(f"\n{introner_type}:")
        print(f"  Total tests: {len(type_results)}")
        print(f"  Significant (FDR < 0.05): {type_results['significant'].sum()}")
        print(f"  Enriched: {(type_results['enrichment_type'] == 'Enriched').sum()}")
        print(f"  Depleted: {(type_results['enrichment_type'] == 'Depleted').sum()}")
        
        # Show top enriched categories
        enriched = type_results[type_results['enrichment_type'] == 'Enriched'].sort_values('fdr')
        if len(enriched) > 0:
            print(f"  Top enriched categories:")
            for _, row in enriched.head(3).iterrows():
                print(f"    {row['go_category']}: OR={row['odds_ratio']:.2f}, FDR={row['fdr']:.3f}")

def save_results(df_results, output_file):
    """Save enrichment results to file."""
    print(f"Saving results to {output_file}...")
    
    # Sort by significance and effect size
    df_sorted = df_results.sort_values(['introner_type', 'fdr', 'fold_enrichment'], 
                                      ascending=[True, True, False])
    
    # Select columns for output
    output_columns = [
        'introner_type', 'go_category', 'genes_with_introner_and_go', 
        'genes_with_introner_total', 'genes_with_go_total', 'background_total',
        'odds_ratio', 'fold_enrichment', 'p_value', 'fdr', 'enrichment_type', 'significant'
    ]
    
    df_output = df_sorted[output_columns]
    df_output.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    
    print(f"  Saved {len(df_output)} enrichment results")

def main():
    parser = argparse.ArgumentParser(
        description='Introner Phase 2 - Functional Enrichment Analysis'
    )
    parser.add_argument('phase1_file', help='Phase 1 gene classification TSV')
    parser.add_argument('output_file', help='Output enrichment results TSV')
    parser.add_argument('--go_json', required=True, help='Path to gene2go JSON file')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                        help='FDR significance threshold (default: 0.05)')
    parser.add_argument('--min_genes', type=int, default=2,
                        help='Minimum genes in GO category to test (default: 2)')
    args = parser.parse_args()

    print("=== INTRONER INSERTION PATTERN ANALYSIS - PHASE 2 ===")
    print(f"Phase 1 results: {args.phase1_file}")
    print(f"GO JSON: {args.go_json}")
    print(f"Output file: {args.output_file}")
    print(f"FDR threshold: {args.fdr_threshold}")
    print(f"Min genes: {args.min_genes}")
    print()

    try:
        # Load Phase 1 results
        df = load_phase1_results(args.phase1_file, args.go_json)

        # Create GO Slim mapping
        go_slim_mapping = create_go_slim_mapping()

        # Map genes to GO Slim categories
        gene_to_slim, slim_to_genes = map_genes_to_go_slim(df, go_slim_mapping)

        # Create introner gene sets
        gene_sets, background_genes = create_introner_gene_sets(df)

        # Perform enrichment analysis
        results = perform_enrichment_analysis(gene_sets, background_genes, slim_to_genes,
                                              min_category_size=args.min_genes)

        # Apply multiple testing correction
        df_results = apply_multiple_testing_correction(results,
                                                       fdr_threshold=args.fdr_threshold)

        # Generate summary statistics
        generate_summary_statistics(df_results)

        # Save results
        save_results(df_results, args.output_file)

        print("\n✓ Phase 2 enrichment analysis completed successfully!")

    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
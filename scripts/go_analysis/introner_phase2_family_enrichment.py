#!/usr/bin/env python3
"""
Introner Insertion Pattern Analysis - Phase 2 (Family-Based)
Functional Enrichment Analysis by Introner Family Type

Performs GO enrichment analysis stratified by introner family (sequence composition type):
- Family 1, Family 2, Family 3, Family 4, Family 14, Family 15
- Excludes Family -1 (missing data placeholder)
- Uses overlapping categories: genes can appear in multiple families

Usage:
python scripts/introner_phase2_family_enrichment.py genotype_matrixes/genotype_matrix_updated_12112025.tsv results/phase1_gene_classification_12112025.tsv results/phase2_family_enrichment_results.tsv
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

def load_genotype_matrix(genotype_file):
    """
    Load genotype matrix and extract family information.
    NOTE: Only families for genes in Phase 1 (GTF-validated) will be used.
    """
    print(f"Loading genotype matrix from {genotype_file}...")

    # Read TSV file, skip duplicate header on line 2
    df = pd.read_csv(genotype_file, sep='\t', header=0)
    df = df[df['ortholog_id'] != 'ortholog_id']  # Remove duplicate header row

    print(f"  Loaded {len(df)} records")
    print(f"  Unique ortholog groups: {df['ortholog_id'].nunique()}")

    # Get unique families (excluding -1)
    families = set(df['family'].unique())
    families.discard(-1)
    families = sorted([f for f in families if isinstance(f, (int, float))])

    print(f"  Found {len(families)} introner families: {families}")

    return df, families

def load_phase1_results(phase1_file, genotype_df, families, go_json):
    """
    Load Phase 1 gene classification results and create family-based gene sets.

    CRITICAL: Only uses genes that Phase 1 identified as having introners (GTF-validated).
    Family information is looked up from genotype matrix, but ONLY for Phase 1 genes.
    """
    print(f"Loading Phase 1 results from {phase1_file}...")

    # Load Phase 1 introner classifications (these are GTF-validated genes)
    phase1_df = pd.read_csv(phase1_file, sep='\t')
    print(f"  Loaded {len(phase1_df)} genes from introner analysis (GTF-validated)")

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
        return df_with_go, {}

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

    # Filter genes with GO annotations
    df_with_go = df_complete[df_complete['go_terms'].apply(len) > 0]

    print(f"  CORRECTED BACKGROUND: {len(df_with_go)} total genes with GO annotations")

    # Build family gene sets from genotype matrix
    # CRITICAL: Only assign families to genes that Phase 1 actually identified
    phase1_genes_set = set(phase1_df['gene_id'])

    family_genes = defaultdict(set)

    for _, row in genotype_df.iterrows():
        gene = row['gene']
        family = row['family']

        # Skip missing data and empty genes
        if pd.isna(family) or family == -1:
            continue
        if pd.isna(gene) or not str(gene).strip():
            continue

        # Clean gene ID
        gene_clean = str(gene).strip()
        if '.3.0.228' in gene_clean:
            gene_clean = gene_clean.replace('.3.0.228', '')

        # CRITICAL: Only include if this gene was detected by Phase 1
        if gene_clean in phase1_genes_set:
            family_genes[int(family)].add(gene_clean)

    print(f"\n  Family gene sets (PHASE 1-VALIDATED ONLY):")
    for family in sorted(family_genes.keys()):
        print(f"    Family {family}: {len(family_genes[family])} genes")

    return df_with_go, family_genes

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

def perform_enrichment_analysis(family_genes, background_genes, slim_to_genes, min_category_size=5):
    """
    Perform Fisher's exact test for GO enrichment analysis by family.
    CRITICAL: Filter family gene sets to only genes with GO annotations (in background_genes).
    """
    print("Performing enrichment analysis by family...")

    results = []

    for family, introner_genes in family_genes.items():
        # CRITICAL: Only include genes that are in the background (have GO annotations)
        introner_genes_set = set(introner_genes) & background_genes
        print(f"  Analyzing Family {family} ({len(introner_genes_set)} genes with GO)...")

        for go_category, category_genes in slim_to_genes.items():
            category_genes_set = set(category_genes)

            # Skip categories with too few genes
            if len(category_genes_set) < min_category_size:
                continue

            # Create 2x2 contingency table
            # [family_with_go, family_without_go, non_family_with_go, non_family_without_go]
            family_with_go = len(introner_genes_set & category_genes_set)
            family_without_go = len(introner_genes_set - category_genes_set)
            non_family_with_go = len(category_genes_set - introner_genes_set)
            non_family_without_go = len(background_genes - introner_genes_set - category_genes_set)

            # Skip if no overlap
            if family_with_go == 0:
                continue

            # Fisher's exact test
            contingency_table = [[family_with_go, family_without_go],
                               [non_family_with_go, non_family_without_go]]

            odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')

            # Calculate fold enrichment
            expected = (len(introner_genes_set) * len(category_genes_set)) / len(background_genes)
            fold_enrichment = family_with_go / expected if expected > 0 else float('inf')

            results.append({
                'family': int(family),
                'go_category': go_category,
                'genes_with_family_and_go': family_with_go,
                'genes_with_family_total': len(introner_genes_set),
                'genes_with_go_total': len(category_genes_set),
                'background_total': len(background_genes),
                'odds_ratio': odds_ratio,
                'fold_enrichment': fold_enrichment,
                'p_value': p_value,
                'contingency_table': str(contingency_table)
            })

    print(f"  Completed {len(results)} enrichment tests")
    return results

def apply_multiple_testing_correction(results):
    """Apply FDR correction for multiple testing."""
    print("Applying multiple testing correction...")

    if not results:
        return results

    df_results = pd.DataFrame(results)

    # Apply FDR correction within each family
    for family in df_results['family'].unique():
        mask = df_results['family'] == family
        p_values = df_results.loc[mask, 'p_value'].values

        if len(p_values) > 0:
            _, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')
            df_results.loc[mask, 'fdr'] = fdr_values

    # Add significance flags
    df_results['significant'] = df_results['fdr'] < 0.05
    df_results['enrichment_type'] = df_results.apply(
        lambda row: 'Enriched' if row['odds_ratio'] > 1 and row['significant']
                   else 'Depleted' if row['odds_ratio'] < 1 and row['significant']
                   else 'Not Significant', axis=1
    )

    return df_results

def generate_summary_statistics(df_results):
    """Generate summary statistics for enrichment results."""
    print("\n=== PHASE 2 (FAMILY-BASED) ENRICHMENT ANALYSIS SUMMARY ===")

    for family in sorted(df_results['family'].unique()):
        family_results = df_results[df_results['family'] == family]

        print(f"\nFamily {int(family)}:")
        print(f"  Total tests: {len(family_results)}")
        print(f"  Significant (FDR < 0.05): {family_results['significant'].sum()}")
        print(f"  Enriched: {(family_results['enrichment_type'] == 'Enriched').sum()}")
        print(f"  Depleted: {(family_results['enrichment_type'] == 'Depleted').sum()}")

        # Show top enriched categories
        enriched = family_results[family_results['enrichment_type'] == 'Enriched'].sort_values('fdr')
        if len(enriched) > 0:
            print(f"  Top enriched categories:")
            for _, row in enriched.head(3).iterrows():
                print(f"    {row['go_category']}: OR={row['odds_ratio']:.2f}, FDR={row['fdr']:.3f}")

def save_results(df_results, output_file):
    """Save enrichment results to file."""
    print(f"Saving results to {output_file}...")

    # Sort by family and significance
    df_sorted = df_results.sort_values(['family', 'fdr', 'fold_enrichment'],
                                      ascending=[True, True, False])

    # Select columns for output
    output_columns = [
        'family', 'go_category', 'genes_with_family_and_go',
        'genes_with_family_total', 'genes_with_go_total', 'background_total',
        'odds_ratio', 'fold_enrichment', 'p_value', 'fdr', 'enrichment_type', 'significant'
    ]

    df_output = df_sorted[output_columns]
    df_output.to_csv(output_file, sep='\t', index=False, float_format='%.6f')

    print(f"  Saved {len(df_output)} enrichment results")

def main():
    parser = argparse.ArgumentParser(
        description='Introner Phase 2 - Family-Based Functional Enrichment Analysis'
    )
    parser.add_argument('genotype_file', help='Genotype matrix TSV')
    parser.add_argument('phase1_file', help='Phase 1 gene classification TSV')
    parser.add_argument('output_file', help='Output enrichment results TSV')
    parser.add_argument('--go_json', required=True, help='Path to gene2go JSON file')
    args = parser.parse_args()

    print("=== INTRONER INSERTION PATTERN ANALYSIS - PHASE 2 (FAMILY-BASED) ===")
    print(f"Genotype matrix: {args.genotype_file}")
    print(f"Phase 1 results: {args.phase1_file}")
    print(f"GO JSON: {args.go_json}")
    print(f"Output file: {args.output_file}")
    print()

    try:
        # Load data
        genotype_df, families = load_genotype_matrix(args.genotype_file)
        df, family_genes = load_phase1_results(args.phase1_file, genotype_df, families, args.go_json)

        # Create GO Slim mapping
        go_slim_mapping = create_go_slim_mapping()

        # Map genes to GO Slim categories
        gene_to_slim, slim_to_genes = map_genes_to_go_slim(df, go_slim_mapping)

        # Perform enrichment analysis
        background_genes = set(df['gene_id'])
        results = perform_enrichment_analysis(family_genes, background_genes, slim_to_genes)

        # Apply multiple testing correction
        df_results = apply_multiple_testing_correction(results)

        # Generate summary statistics
        generate_summary_statistics(df_results)

        # Save results
        save_results(df_results, args.output_file)

        print("\n✓ Phase 2 (family-based) enrichment analysis completed successfully!")

    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

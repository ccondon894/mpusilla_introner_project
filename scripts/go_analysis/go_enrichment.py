import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
import os

def bonferroni_correction(p_values):
    """
    Apply Bonferroni correction to a list of p-values.
    
    Parameters:
    -----------
    p_values : list or array
        List of uncorrected p-values
        
    Returns:
    --------
    list : Bonferroni corrected p-values (capped at 1.0)
    """
    p_values = np.array(p_values)
    n_tests = len(p_values)
    corrected = p_values * n_tests
    # Cap at 1.0 since p-values cannot exceed 1
    return np.minimum(corrected, 1.0)

def fdr_correction(p_values):
    """
    Apply Benjamini-Hochberg FDR correction to a list of p-values.
    
    Parameters:
    -----------
    p_values : list or array
        List of uncorrected p-values
        
    Returns:
    --------
    list : FDR corrected p-values
    """
    p_values = np.array(p_values)
    n_tests = len(p_values)
    
    # Sort p-values and keep track of original indices
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Calculate FDR corrections
    corrected = np.zeros_like(sorted_p_values)
    
    # Work backwards from largest p-value
    for i in range(n_tests - 1, -1, -1):
        rank = i + 1  # rank starts from 1
        bh_value = sorted_p_values[i] * n_tests / rank
        
        if i == n_tests - 1:
            corrected[i] = min(bh_value, 1.0)
        else:
            corrected[i] = min(bh_value, corrected[i + 1])
    
    # Restore original order
    original_order = np.argsort(sorted_indices)
    return corrected[original_order]

def perform_go_enrichment(go_terms_json, genes_of_interest_tsv, output_prefix='go_enrichment', obo_file=None):
    """
    Perform GO enrichment analysis on genes of interest using expanded GO terms.

    Parameters:
    -----------
    go_terms_json : str
        Path to JSON file containing gene-to-GO mappings
    genes_of_interest_tsv : str
        Path to TSV file containing genes of interest with relationship column
    output_prefix : str
        Prefix for output files
    obo_file : str
        Path to GO OBO ontology file
    """
    # Step 1: Read the expanded GO terms
    print("Loading GO annotations...")
    with open(go_terms_json, 'r') as f:
        gene_to_go = json.load(f)
    
    # Ensure GO terms are properly formatted (with GO: prefix)
    formatted_gene_to_go = {}
    for gene, terms in gene_to_go.items():
        formatted_gene_to_go[gene] = []
        for term in terms:
            if not term.startswith('GO:'):
                term = f"GO:{term}"
            formatted_gene_to_go[gene].append(term)
    
    # Step 2: Read and filter the genes of interest
    print("Processing genes of interest...")
    genes_df = pd.read_csv(genes_of_interest_tsv, sep='\t')
    filtered_genes = genes_df[genes_df['relationship'].isin(['overlapping'])]
    genes_of_interest = [gene.split(".3.0.228")[0] for gene in filtered_genes['gene_id']]
    genes_of_interest = set(genes_of_interest)
    
    print(f"Total genes in input file: {len(genes_df)}")
    print(f"Genes filtered by relationship: {len(filtered_genes)}")
    print(f"Unique genes of interest: {len(genes_of_interest)}")
    
    # Step 3: Prepare data for goatools
    # Create population (all genes with GO terms)
    population = set(formatted_gene_to_go.keys())
    
    # Find study genes (genes of interest) that have GO annotations
    study_genes = genes_of_interest.intersection(population)
    
    print(f"Total genes with GO annotations: {len(population)}")
    print(f"Genes of interest with GO annotations: {len(study_genes)}")
    
    if len(study_genes) == 0:
        print("No genes of interest have GO annotations. Analysis cannot proceed.")
        return None
    
    # Step 4: Create gene-to-GO association dictionary
    associations = {}
    for gene, terms in formatted_gene_to_go.items():
        if gene in population:
            associations[gene] = set(terms)
    
    # Step 5: Run GO enrichment analysis
    print("Running GO enrichment analysis...")
    obodag = GODag(obo_file)
    
    # Create the GO enrichment study object (without multiple testing corrections)
    g = GOEnrichmentStudy(
        pop=population,
        assoc=associations,
        obo_dag=obodag,
        propagate_counts=True,  # propagate counts up the GO hierarchy
        alpha=0.05  # significance level
    )
    
    # Run the enrichment analysis
    results = g.run_study(study_genes)
    
    # Step 7: Process and save results
    print("Processing results...")
    # Convert results to dataframe for easier manipulation
    rows = []
    p_values = []
    
    for r in results:
        # Get namespace from GODag using the GO term
        namespace = "Unknown"
        if r.GO in obodag:
            namespace = obodag[r.GO].namespace
            
        rows.append({
            'GO_term': r.GO,
            'term_name': r.name,
            'namespace': namespace,
            'p_uncorrected': r.p_uncorrected,
            'ratio_in_study': f"{r.study_count}/{r.study_n}",
            'ratio_in_pop': f"{r.pop_count}/{r.pop_n}",
            'depth': r.depth,
            'study_items': ','.join(r.study_items),
            'enrichment': 'enriched' if r.enrichment == 'e' else 'depleted'  # 'e' for enriched, 'd' for depleted
        })
        p_values.append(r.p_uncorrected)
    
    df_results = pd.DataFrame(rows)
    
    # Apply custom multiple testing corrections
    print("Applying multiple testing corrections...")
    df_results['p_bonferroni'] = bonferroni_correction(p_values)
    df_results['p_fdr'] = fdr_correction(p_values)
    
    # Sort by p-value
    df_results = df_results.sort_values('p_uncorrected')
    
    # Save to CSV
    output_file = output_prefix if output_prefix.endswith('.tsv') else f"{output_prefix}_results.tsv"
    df_results.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    print(f"\nGO Enrichment Analysis Summary:")
    print(f"--------------------------------")
    print(f"Total GO terms tested: {len(df_results)}")
    print(f"Significant GO terms (p<0.05, uncorrected): {len(df_results[df_results['p_uncorrected'] < 0.05])}")
    print(f"Significant GO terms (p<0.05, Bonferroni): {len(df_results[df_results['p_bonferroni'] < 0.05])}")
    print(f"Significant GO terms (p<0.05, FDR): {len(df_results[df_results['p_fdr'] < 0.05])}")
    print(f"\nResults saved to: {output_file}")
    
    # Create visualizations - derive prefix from output file path
    viz_prefix = output_file.rsplit('.tsv', 1)[0] if output_file.endswith('.tsv') else output_file
    visualize_go_results(df_results, viz_prefix)
    
    return df_results

def visualize_go_results(df_results, output_prefix, top_n=20):
    """
    Create visualizations for GO enrichment results.
    
    Parameters:
    -----------
    df_results : pandas.DataFrame
        DataFrame containing GO enrichment results
    output_prefix : str
        Prefix for output files
    top_n : int
        Number of top GO terms to visualize
    """
    # Filter for significant results and take top N
    sig_results = df_results[df_results['p_fdr'] < 0.05]
    
    if len(sig_results) == 0:
        print("No significant GO terms to visualize.")
        return
    
    # Sort and take top N
    sig_results = sig_results.sort_values('p_fdr').head(top_n)
    
    # Split by namespace for better visualization
    namespaces = sig_results['namespace'].unique()
    
    for namespace in namespaces:
        ns_results = sig_results[sig_results['namespace'] == namespace]
        if len(ns_results) == 0:
            continue
        
        # Create barplot for this namespace
        plt.figure(figsize=(10, max(6, len(ns_results) * 0.3)))
        
        # Calculate -log10(p-value) for better visualization
        ns_results['neg_log_p'] = -np.log10(ns_results['p_fdr'])
        
        # Extract study ratios for calculating enrichment fold
        ns_results['study_count'] = ns_results['ratio_in_study'].apply(lambda x: int(x.split('/')[0]))
        ns_results['study_total'] = ns_results['ratio_in_study'].apply(lambda x: int(x.split('/')[1]))
        ns_results['pop_count'] = ns_results['ratio_in_pop'].apply(lambda x: int(x.split('/')[0]))
        ns_results['pop_total'] = ns_results['ratio_in_pop'].apply(lambda x: int(x.split('/')[1]))
        
        # Calculate fold enrichment
        ns_results['fold_enrichment'] = (ns_results['study_count'] / ns_results['study_total']) / \
                                        (ns_results['pop_count'] / ns_results['pop_total'])
        
        # Sort for better visualization
        ns_results = ns_results.sort_values('neg_log_p')
        
        # Create barplot
        sns.barplot(
            data=ns_results,
            y='term_name',
            x='neg_log_p',
            color=sns.color_palette("Set2")[list(namespaces).index(namespace)],
            alpha=0.8
        )
        
        plt.xlabel('-log10(FDR p-value)')
        plt.ylabel('GO Term')
        plt.title(f'Top Enriched GO Terms - {namespace}')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{namespace}_barplot.png", dpi=300)
        plt.close()
        
        # Create bubble plot for fold enrichment vs significance
        plt.figure(figsize=(12, 8))
        
        # Create scatter plot with size proportional to gene count
        plt.scatter(
            x=ns_results['fold_enrichment'],
            y=ns_results['neg_log_p'],
            s=ns_results['study_count'] * 10,
            alpha=0.7,
            color=sns.color_palette("Set2")[list(namespaces).index(namespace)]
        )
        
        # Add labels
        for _, row in ns_results.iterrows():
            plt.annotate(
                row['term_name'],
                xy=(row['fold_enrichment'], row['neg_log_p']),
                xytext=(5, 0),
                textcoords='offset points',
                ha='left',
                va='center',
                fontsize=8
            )
        
        plt.xlabel('Fold Enrichment')
        plt.ylabel('-log10(FDR p-value)')
        plt.title(f'GO Term Enrichment - {namespace}')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{namespace}_bubble.png", dpi=300)
        plt.close()
    
    print(f"Visualizations saved with prefix: {output_prefix}")

# Example usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Perform GO enrichment analysis on genes of interest.')
    parser.add_argument('--go_json', required=True, help='Path to JSON file with gene-to-GO mappings')
    parser.add_argument('--genes', required=True, help='Path to TSV file with genes of interest')
    parser.add_argument('--go_obo', required=True, help='Path to GO OBO ontology file')
    parser.add_argument('--output', default='go_enrichment', help='Prefix for output files')

    args = parser.parse_args()

    results = perform_go_enrichment(args.go_json, args.genes, args.output, args.go_obo)
    
    if results is not None:
        print("\nTop 10 enriched GO terms:")
        print(results[['GO_term', 'term_name', 'namespace', 'p_uncorrected', 'p_fdr', 'ratio_in_study', 'enrichment']].head(10))
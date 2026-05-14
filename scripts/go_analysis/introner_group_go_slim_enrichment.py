#!/usr/bin/env python3
"""
GO Slim enrichment for introner groups derived from genotype_matrix.final.tsv.

This companion analysis uses the same introner gene sets as the full GO test,
but collapses propagated GO annotations onto a small set of broad GO Slim
categories. That reduces the multiple-testing burden and improves power for
general functional patterns.
"""

import argparse
import math
import sys
from collections import defaultdict

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from introner_group_go_enrichment import (
    INTRONER_GROUPS,
    build_background,
    classify_loci,
    load_gene2go,
    parse_csv_arg,
    parse_obo,
    propagate_gene2go,
)


GO_SLIM_TERMS = {
    # Biological process
    "GO:0008152": "Metabolic Process",
    "GO:0009058": "Biosynthetic Process",
    "GO:0009056": "Catabolic Process",
    "GO:0006810": "Transport",
    "GO:0051234": "Establishment of Localization",
    "GO:0006412": "Translation",
    "GO:0006351": "Transcription DNA-templated",
    "GO:0006355": "Regulation of Transcription",
    "GO:0006396": "RNA Processing",
    "GO:0008219": "Cell Death",
    "GO:0007049": "Cell Cycle",
    "GO:0000278": "Mitotic Cell Cycle",
    "GO:0006464": "Cellular Protein Modification",
    "GO:0016043": "Cellular Component Organization",
    "GO:0050896": "Response to Stimulus",
    "GO:0009607": "Response to Biotic Stimulus",
    "GO:0009628": "Response to Abiotic Stimulus",
    "GO:0009719": "Response to Endogenous Stimulus",
    "GO:0015979": "Photosynthesis",
    "GO:0009765": "Photosynthesis Light Reactions",
    "GO:0019684": "Photosynthesis Light Harvesting",
    # Molecular function
    "GO:0003824": "Catalytic Activity",
    "GO:0016787": "Hydrolase Activity",
    "GO:0016740": "Transferase Activity",
    "GO:0016491": "Oxidoreductase Activity",
    "GO:0016874": "Ligase Activity",
    "GO:0005515": "Protein Binding",
    "GO:0003677": "DNA Binding",
    "GO:0003723": "RNA Binding",
    "GO:0022857": "Transmembrane Transporter Activity",
    "GO:0005524": "ATP Binding",
    "GO:0000166": "Nucleotide Binding",
    "GO:0046872": "Metal Ion Binding",
    # Cellular component
    "GO:0005634": "Nucleus",
    "GO:0005737": "Cytoplasm",
    "GO:0005739": "Mitochondrion",
    "GO:0009507": "Chloroplast",
    "GO:0005783": "Endoplasmic Reticulum",
    "GO:0005794": "Golgi Apparatus",
    "GO:0005886": "Plasma Membrane",
    "GO:0016020": "Membrane",
    "GO:0005622": "Intracellular",
    "GO:0005840": "Ribosome",
    "GO:0005856": "Cytoskeleton",
    "GO:0005773": "Vacuole",
}


def build_slim_to_genes(background_genes, gene2go, terms, min_genes):
    slim_to_genes = defaultdict(set)
    for gene in background_genes:
        for go_id in gene2go[gene]:
            if go_id in GO_SLIM_TERMS and go_id in terms:
                slim_to_genes[go_id].add(gene)
    return {
        go_id: genes
        for go_id, genes in slim_to_genes.items()
        if len(genes) >= min_genes
    }


def run_slim_enrichment(gene_set_df, background_genes, slim_to_genes, terms, fdr_threshold):
    rows = []
    for group_name in INTRONER_GROUPS:
        if gene_set_df.empty:
            study_genes = set()
        else:
            study_genes = set(gene_set_df.loc[gene_set_df["introner_group"] == group_name, "gene_id"])
        study_genes &= background_genes
        nonstudy_genes = background_genes - study_genes

        for go_id, category_genes in slim_to_genes.items():
            study_with_go = len(study_genes & category_genes)
            study_without_go = len(study_genes - category_genes)
            background_with_go = len(nonstudy_genes & category_genes)
            background_without_go = len(nonstudy_genes - category_genes)
            table = [[study_with_go, study_without_go], [background_with_go, background_without_go]]

            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
            expected = len(study_genes) * len(category_genes) / len(background_genes) if background_genes else 0
            fold_enrichment = study_with_go / expected if expected else math.nan
            enrichment_type = (
                "enriched" if odds_ratio > 1
                else "depleted" if odds_ratio < 1
                else "neutral"
            )

            rows.append({
                "introner_group": group_name,
                "GO_slim_term": go_id,
                "go_slim_category": GO_SLIM_TERMS[go_id],
                "namespace": terms[go_id].get("namespace", ""),
                "study_with_go": study_with_go,
                "study_total": len(study_genes),
                "background_with_go": len(category_genes),
                "background_total": len(background_genes),
                "odds_ratio": odds_ratio,
                "fold_enrichment": fold_enrichment,
                "p_value": p_value,
                "enrichment_type": enrichment_type,
                "study_genes": ";".join(sorted(study_genes & category_genes)),
            })

    results = pd.DataFrame(rows)
    if results.empty:
        return results

    _, results["fdr_bh_global"], _, _ = multipletests(results["p_value"], method="fdr_bh")
    results["fdr_bh_by_group"] = math.nan
    for group_name, idx in results.groupby("introner_group").groups.items():
        _, corrected, _, _ = multipletests(results.loc[idx, "p_value"], method="fdr_bh")
        results.loc[idx, "fdr_bh_by_group"] = corrected
    results["significant_global"] = (
        (results["fdr_bh_global"] < fdr_threshold)
        & (results["enrichment_type"] != "neutral")
    )
    results["significant_by_group"] = (
        (results["fdr_bh_by_group"] < fdr_threshold)
        & (results["enrichment_type"] != "neutral")
    )
    results["significant"] = results["significant_global"]
    return results.sort_values(
        ["introner_group", "fdr_bh_global", "p_value", "fold_enrichment"],
        ascending=[True, True, True, False],
    )


def write_summary(path, genotype_file, go_jsons, group1_samples, group2_samples,
                  locus_df, gene_set_df, background_genes, slim_to_genes,
                  results, fdr_threshold, background_scope):
    with open(path, "w") as handle:
        handle.write("Introner Group GO Slim Enrichment Summary\n")
        handle.write("=" * 45 + "\n\n")
        handle.write(f"Genotype matrix: {genotype_file}\n")
        handle.write(f"GO JSON sources: {', '.join(go_jsons)}\n")
        handle.write(f"Group 1 samples: {', '.join(group1_samples)}\n")
        handle.write(f"Group 2 samples: {', '.join(group2_samples)}\n")
        handle.write("Fisher alternative: two-sided (enrichment and depletion)\n")
        handle.write(f"Background scope: {background_scope}\n")
        handle.write("Multiple testing: Benjamini-Hochberg across all group-category tests; q-value is the BH-adjusted p-value/FDR estimate.\n")
        handle.write(f"Primary significance flag: global q < {fdr_threshold}.\n\n")

        handle.write(f"Ortholog groups classified: {len(locus_df)}\n")
        handle.write(f"Background genes with GO terms: {len(background_genes)}\n")
        handle.write(f"GO Slim categories tested: {len(slim_to_genes)}\n\n")

        for group_name in INTRONER_GROUPS:
            loci = locus_df["introner_groups"].fillna("").str.contains(group_name, regex=False).sum()
            if gene_set_df.empty:
                group_genes = set()
            else:
                group_genes = set(gene_set_df.loc[gene_set_df["introner_group"] == group_name, "gene_id"])
            annotated_genes = group_genes & background_genes
            sig = 0 if results.empty else (
                (results["introner_group"] == group_name) & (results["significant"])
            ).sum()
            enriched = 0 if results.empty else (
                (results["introner_group"] == group_name)
                & (results["significant"])
                & (results["enrichment_type"] == "enriched")
            ).sum()
            depleted = 0 if results.empty else (
                (results["introner_group"] == group_name)
                & (results["significant"])
                & (results["enrichment_type"] == "depleted")
            ).sum()
            handle.write(
                f"{group_name}: {loci} loci, {len(group_genes)} genes, "
                f"{len(annotated_genes)} GO-annotated tested genes, "
                f"{sig} significant GO Slim categories at global q<{fdr_threshold} "
                f"({enriched} enriched, {depleted} depleted)\n"
            )


def main():
    parser = argparse.ArgumentParser(description="GO Slim enrichment for final-matrix introner groups.")
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--go-json", required=True, nargs="+")
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--group1-samples", required=True)
    parser.add_argument("--group2-samples", required=True)
    parser.add_argument("--enrichment-output", required=True)
    parser.add_argument("--significant-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--background-scope", default="all_annotated",
                        choices=["all_annotated", "genotyped"],
                        help="Use all GO-annotated genes or only genes represented in the genotype matrix as background.")
    parser.add_argument("--exclude-within-status", default="low_identity")
    args = parser.parse_args()

    group1_samples = parse_csv_arg(args.group1_samples)
    group2_samples = parse_csv_arg(args.group2_samples)
    excluded_within_status = parse_csv_arg(args.exclude_within_status)

    terms = parse_obo(args.go_obo)
    gene2go_raw = load_gene2go(args.go_json)
    gene2go = propagate_gene2go(gene2go_raw, terms)
    locus_df, gene_set_df = classify_loci(
        args.genotype_matrix,
        group1_samples,
        group2_samples,
        excluded_within_status,
    )
    background_genes = build_background(locus_df, gene2go, args.background_scope)
    if not background_genes:
        raise ValueError("No genotyped genes have GO annotations after propagation.")

    slim_to_genes = build_slim_to_genes(background_genes, gene2go, terms, args.min_genes)
    if not slim_to_genes:
        raise ValueError("No GO Slim categories met the minimum gene threshold.")

    results = run_slim_enrichment(
        gene_set_df,
        background_genes,
        slim_to_genes,
        terms,
        args.fdr_threshold,
    )
    results.to_csv(args.enrichment_output, sep="\t", index=False, float_format="%.6g")
    significant = results[results["significant"]] if not results.empty else results
    significant.to_csv(args.significant_output, sep="\t", index=False, float_format="%.6g")
    write_summary(
        args.summary_output,
        args.genotype_matrix,
        args.go_json,
        group1_samples,
        group2_samples,
        locus_df,
        gene_set_df,
        background_genes,
        slim_to_genes,
        results,
        args.fdr_threshold,
        args.background_scope,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

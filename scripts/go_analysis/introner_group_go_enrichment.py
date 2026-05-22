#!/usr/bin/env python3
"""
GO enrichment for introner groups derived directly from genotype_matrix.final.tsv.

The analysis classifies introner ortholog groups into biologically interpretable
sets, maps those sets to genes, propagates gene GO annotations to ontology
ancestors, and tests one-sided GO enrichment against the genotyped gene universe.
"""

import argparse
import json
import math
import sys
from collections import defaultdict
from functools import lru_cache

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


PRESENT = 1
ABSENT = 2
MISSING = 3

ANCESTRAL_STATUS = {"ancestral", "likely_ancestral", "ancestral_low_identity"}
INDEPENDENT_STATUS = {"independent", "likely_independent"}
BASE_INTRONER_GROUPS = [
    "all_introners",
    "ancestral",
    "independent_insertion",
    "polymorphic_within_group1",
    "consistent_group1",
    "consistent_group2",
]
INTRONER_GROUPS = BASE_INTRONER_GROUPS


def clean_gene_id(gene_id):
    if pd.isna(gene_id):
        return None
    gene_id = str(gene_id).strip()
    if not gene_id or gene_id.lower() == "nan":
        return None
    return gene_id.replace(".3.0.228", "")


def parse_presence(value):
    if pd.isna(value):
        return None
    try:
        presence = int(float(value))
    except (TypeError, ValueError):
        return None
    if presence in {PRESENT, ABSENT, MISSING}:
        return presence
    return None


def parse_family(value):
    if pd.isna(value):
        return None
    try:
        family = int(float(value))
    except (TypeError, ValueError):
        return None
    if family < 0:
        return None
    return family


def clean_go_id(go_id):
    go_id = str(go_id).strip()
    if not go_id:
        return None
    if go_id.startswith("GO:"):
        return go_id
    if go_id.isdigit():
        return f"GO:{go_id.zfill(7)}"
    return go_id


def parse_csv_arg(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_obo(obo_file):
    terms = {}
    current = None

    def commit(term):
        if term and term.get("id") and not term.get("is_obsolete"):
            terms[term["id"]] = term

    with open(obo_file) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                commit(current)
                current = {"parents": set()}
                continue
            if line.startswith("["):
                commit(current)
                current = None
                continue
            if current is None or not line:
                continue
            if line.startswith("id: "):
                current["id"] = line.split("id: ", 1)[1]
            elif line.startswith("name: "):
                current["name"] = line.split("name: ", 1)[1]
            elif line.startswith("namespace: "):
                current["namespace"] = line.split("namespace: ", 1)[1]
            elif line.startswith("is_a: "):
                current["parents"].add(line.split("is_a: ", 1)[1].split()[0])
            elif line.startswith("relationship: part_of "):
                current["parents"].add(line.split("relationship: part_of ", 1)[1].split()[0])
            elif line == "is_obsolete: true":
                current["is_obsolete"] = True
    commit(current)
    return terms


def load_gene2go(paths):
    gene2go = defaultdict(set)
    for path in paths:
        with open(path) as handle:
            data = json.load(handle)
        for gene, terms in data.items():
            gene_id = clean_gene_id(gene)
            if not gene_id:
                continue
            for term in terms:
                go_id = clean_go_id(term)
                if go_id:
                    gene2go[gene_id].add(go_id)
    return dict(gene2go)


def propagate_gene2go(gene2go, terms):
    @lru_cache(maxsize=None)
    def ancestors(go_id):
        if go_id not in terms:
            return frozenset()
        found = set()
        for parent in terms[go_id]["parents"]:
            if parent in terms:
                found.add(parent)
                found.update(ancestors(parent))
        return frozenset(found)

    propagated = {}
    for gene, go_ids in gene2go.items():
        expanded = set()
        for go_id in go_ids:
            if go_id not in terms:
                continue
            expanded.add(go_id)
            expanded.update(ancestors(go_id))
        if expanded:
            propagated[gene] = expanded
    return propagated


def first_nonempty(values):
    for value in values:
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return ""


def validate_samples(df, group1_samples, group2_samples):
    matrix_samples = set(df["sample"].dropna().astype(str))
    requested = set(group1_samples) | set(group2_samples)
    missing = sorted(requested - matrix_samples)
    if missing:
        raise ValueError(
            "Configured samples are absent from genotype matrix: "
            + ", ".join(missing)
            + ". Matrix samples: "
            + ", ".join(sorted(matrix_samples))
        )
    overlap = set(group1_samples) & set(group2_samples)
    if overlap:
        raise ValueError("Samples cannot be in both groups: " + ", ".join(sorted(overlap)))


def family_group_label(family):
    return f"family_{family}"


def ordered_group_names(gene_set_df=None):
    family_groups = []
    if gene_set_df is not None and not gene_set_df.empty:
        observed = set(gene_set_df["introner_group"].dropna().astype(str))
        family_groups = sorted(
            (group for group in observed if group.startswith("family_")),
            key=lambda group: int(group.split("_", 1)[1]),
        )
    return BASE_INTRONER_GROUPS + family_groups


def classify_loci(genotype_file, group1_samples, group2_samples, excluded_within_status):
    df = pd.read_csv(genotype_file, sep="\t")
    required = {
        "ortholog_id", "sample", "gene", "presence",
        "within_group_status", "cross_group_status",
    }
    missing_columns = sorted(required - set(df.columns))
    if missing_columns:
        raise ValueError("Missing required columns: " + ", ".join(missing_columns))

    validate_samples(df, group1_samples, group2_samples)
    group1_samples = set(group1_samples)
    group2_samples = set(group2_samples)
    excluded_within_status = set(excluded_within_status)

    locus_rows = []
    group_to_genes = defaultdict(lambda: defaultdict(set))
    has_family = "family" in df.columns

    for ortholog_id, locus in df.groupby("ortholog_id", sort=False):
        genes = sorted(
            {str(gene) for gene in locus["gene"].map(clean_gene_id) if gene}
        )
        within_status = first_nonempty(locus["within_group_status"])
        cross_status = first_nonempty(locus["cross_group_status"])

        group1_all_calls = [
            parse_presence(p)
            for p in locus.loc[locus["sample"].isin(group1_samples), "presence"]
        ]
        group2_all_calls = [
            parse_presence(p)
            for p in locus.loc[locus["sample"].isin(group2_samples), "presence"]
        ]
        group1_calls = [p for p in group1_all_calls if p in {PRESENT, ABSENT}]
        group2_calls = [p for p in group2_all_calls if p in {PRESENT, ABSENT}]

        group1_present = PRESENT in group1_calls
        group2_present = PRESENT in group2_calls
        group1_polymorphic = PRESENT in group1_calls and ABSENT in group1_calls
        group1_consistent_present = (
            len(group1_all_calls) == len(group1_samples)
            and all(p == PRESENT for p in group1_all_calls)
        )
        group2_consistent_present = (
            len(group2_all_calls) == len(group2_samples)
            and all(p == PRESENT for p in group2_all_calls)
        )
        use_locus = within_status not in excluded_within_status

        labels = []
        if use_locus and (group1_present or group2_present):
            labels.append("all_introners")
        if use_locus and cross_status in ANCESTRAL_STATUS and group1_present and group2_present:
            labels.append("ancestral")
        if use_locus and cross_status in INDEPENDENT_STATUS and group1_present and group2_present:
            labels.append("independent_insertion")
        if use_locus and group1_polymorphic:
            labels.append("polymorphic_within_group1")
        if use_locus and group1_consistent_present:
            labels.append("consistent_group1")
        if use_locus and group2_consistent_present:
            labels.append("consistent_group2")

        label_to_genes = {label: set(genes) for label in labels}
        if use_locus and has_family:
            present_rows = locus.loc[locus["presence"].map(parse_presence) == PRESENT]
            for family, family_rows in present_rows.groupby(
                present_rows["family"].map(parse_family), dropna=True
            ):
                if family is None:
                    continue
                family_label = family_group_label(int(family))
                family_genes = {
                    gene for gene in family_rows["gene"].map(clean_gene_id) if gene
                }
                if not family_genes:
                    family_genes = set(genes)
                if family_genes:
                    labels.append(family_label)
                    label_to_genes[family_label] = family_genes

        locus_rows.append({
            "ortholog_id": ortholog_id,
            "gene_ids": ";".join(genes),
            "within_group_status": within_status,
            "cross_group_status": cross_status,
            "group1_nonmissing": len(group1_calls),
            "group1_present": sum(1 for p in group1_calls if p == PRESENT),
            "group1_absent": sum(1 for p in group1_calls if p == ABSENT),
            "group1_missing": sum(1 for p in group1_all_calls if p == MISSING or p is None),
            "group1_consistent_present": group1_consistent_present,
            "group2_nonmissing": len(group2_calls),
            "group2_present": sum(1 for p in group2_calls if p == PRESENT),
            "group2_absent": sum(1 for p in group2_calls if p == ABSENT),
            "group2_missing": sum(1 for p in group2_all_calls if p == MISSING or p is None),
            "group2_consistent_present": group2_consistent_present,
            "introner_groups": ";".join(labels),
            "excluded_by_within_status": not use_locus,
        })

        for label in labels:
            for gene in label_to_genes.get(label, genes):
                group_to_genes[label][gene].add(ortholog_id)

    gene_set_rows = []
    for label in ordered_group_names(pd.DataFrame(
        {"introner_group": list(group_to_genes.keys())}
    )):
        for gene, orthologs in sorted(
            group_to_genes[label].items(),
            key=lambda item: str(item[0]),
        ):
            if gene is None or pd.isna(gene):
                continue
            gene_set_rows.append({
                "introner_group": label,
                "gene_id": gene,
                "n_loci": len(orthologs),
                "ortholog_ids": ";".join(sorted(orthologs)),
            })

    return pd.DataFrame(locus_rows), pd.DataFrame(gene_set_rows)


def build_background(locus_df, gene2go, scope="all_annotated"):
    if scope == "all_annotated":
        return {gene for gene, go_terms in gene2go.items() if go_terms}
    if scope != "genotyped":
        raise ValueError(f"Unsupported background scope: {scope}")
    genes_in_matrix = set()
    for value in locus_df["gene_ids"]:
        if pd.notna(value) and str(value).strip():
            genes_in_matrix.update(str(value).split(";"))
    return {gene for gene in genes_in_matrix if gene in gene2go and gene2go[gene]}


def run_enrichment(gene_set_df, background_genes, gene2go, terms, min_genes, fdr_threshold):
    go_to_genes = defaultdict(set)
    for gene in background_genes:
        for go_id in gene2go[gene]:
            go_to_genes[go_id].add(gene)

    test_terms = {
        go_id: genes
        for go_id, genes in go_to_genes.items()
        if len(genes) >= min_genes and go_id in terms
    }

    rows = []
    for group_name in ordered_group_names(gene_set_df):
        if gene_set_df.empty:
            study_genes = set()
        else:
            study_genes = set(gene_set_df.loc[gene_set_df["introner_group"] == group_name, "gene_id"])
        study_genes &= background_genes
        nonstudy_genes = background_genes - study_genes

        for go_id, annotated_genes in test_terms.items():
            study_with_go = len(study_genes & annotated_genes)
            study_without_go = len(study_genes - annotated_genes)
            background_with_go = len(nonstudy_genes & annotated_genes)
            background_without_go = len(nonstudy_genes - annotated_genes)
            table = [[study_with_go, study_without_go], [background_with_go, background_without_go]]

            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
            expected = len(study_genes) * len(annotated_genes) / len(background_genes) if background_genes else 0
            fold_enrichment = study_with_go / expected if expected else math.nan
            enrichment_type = (
                "enriched" if odds_ratio > 1
                else "depleted" if odds_ratio < 1
                else "neutral"
            )

            rows.append({
                "introner_group": group_name,
                "GO_term": go_id,
                "term_name": terms[go_id].get("name", ""),
                "namespace": terms[go_id].get("namespace", ""),
                "study_with_go": study_with_go,
                "study_total": len(study_genes),
                "background_with_go": len(annotated_genes),
                "background_total": len(background_genes),
                "odds_ratio": odds_ratio,
                "fold_enrichment": fold_enrichment,
                "p_value": p_value,
                "enrichment_type": enrichment_type,
                "study_genes": ";".join(sorted(study_genes & annotated_genes)),
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
                  locus_df, gene_set_df, background_genes, results, fdr_threshold,
                  background_scope):
    with open(path, "w") as handle:
        handle.write("Introner Group GO Enrichment Summary\n")
        handle.write("=" * 40 + "\n\n")
        handle.write(f"Genotype matrix: {genotype_file}\n")
        handle.write(f"GO JSON sources: {', '.join(go_jsons)}\n")
        handle.write(f"Group 1 samples: {', '.join(group1_samples)}\n")
        handle.write(f"Group 2 samples: {', '.join(group2_samples)}\n")
        handle.write("Presence coding: 1=present, 2=absent, 3=missing/excluded\n")
        handle.write("Polymorphic calls ignore missing data; consistent groups require every configured sample to be present.\n")
        handle.write("Fisher alternative: two-sided (enrichment and depletion)\n")
        handle.write(f"Background scope: {background_scope}\n")
        handle.write("Multiple testing: Benjamini-Hochberg across all group-term tests; q-value is the BH-adjusted p-value/FDR estimate.\n")
        handle.write(f"Primary significance flag: global q < {fdr_threshold}.\n\n")

        handle.write(f"Ortholog groups classified: {len(locus_df)}\n")
        handle.write(f"Background genes with GO terms: {len(background_genes)}\n\n")

        for group_name in ordered_group_names(gene_set_df):
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
                f"{sig} significant GO terms at global q<{fdr_threshold} "
                f"({enriched} enriched, {depleted} depleted)\n"
            )


def main():
    parser = argparse.ArgumentParser(description="GO enrichment for final-matrix introner groups.")
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--go-json", required=True, nargs="+",
                        help="One or more gene-to-GO JSON files; terms are unioned by gene.")
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--group1-samples", required=True,
                        help="Comma-separated Group 1 sample IDs.")
    parser.add_argument("--group2-samples", required=True,
                        help="Comma-separated Group 2 sample IDs.")
    parser.add_argument("--locus-output", required=True)
    parser.add_argument("--gene-set-output", required=True)
    parser.add_argument("--enrichment-output", required=True)
    parser.add_argument("--significant-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--background-scope", default="all_annotated",
                        choices=["all_annotated", "genotyped"],
                        help="Use all GO-annotated genes or only genes represented in the genotype matrix as background.")
    parser.add_argument("--exclude-within-status", default="low_identity",
                        help="Comma-separated within_group_status values to exclude.")
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

    results = run_enrichment(
        gene_set_df,
        background_genes,
        gene2go,
        terms,
        args.min_genes,
        args.fdr_threshold,
    )

    locus_df.to_csv(args.locus_output, sep="\t", index=False)
    gene_set_df.to_csv(args.gene_set_output, sep="\t", index=False)
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

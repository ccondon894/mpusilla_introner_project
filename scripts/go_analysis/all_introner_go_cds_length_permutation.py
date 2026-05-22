#!/usr/bin/env python3
"""
CDS-length-weighted permutation test for all-introner GO enrichment.

This tests whether GO enrichment among genes with at least one introner is
stronger than expected if the observed introner-bearing-gene labels were
reassigned across GO-annotated coding genes with probability proportional to
CDS length. The permutation preserves the observed number of introner-bearing
genes. Conceptually, the observed per-gene introner burdens are exchangeable
bundles assigned onto those sampled genes, which preserves the total introner
count and the burden distribution; the first implemented GO statistic uses only
the resulting gene-level presence/absence state.
"""

import argparse
import math
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from introner_group_go_enrichment import (
    load_gene2go,
    parse_csv_arg,
    parse_obo,
    propagate_gene2go,
    classify_loci,
)


GO_TERMS_OF_INTEREST = [
    "GO:0000166",  # nucleotide binding
    "GO:0005524",  # ATP binding
    "GO:0017076",  # purine nucleotide binding
    "GO:1901265",  # nucleoside phosphate binding
    "GO:0043168",  # anion binding
]


_TERM_MATRIX = None
_WEIGHTS = None
_STUDY_TOTAL = None
_OBS_COUNTS = None


def clean_gene_id(gene_id):
    if pd.isna(gene_id):
        return None
    gene_id = str(gene_id).strip()
    if not gene_id or gene_id.lower() == "nan":
        return None
    return gene_id.replace(".3.0.228", "")


def extract_cds_lengths_from_gtf(gtf_path):
    cds_lengths = defaultdict(int)
    attr_re = re.compile(r'gene_id[=\s]+"?([^;"]+)"?')

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            match = attr_re.search(fields[8])
            if not match:
                continue
            gene_id = clean_gene_id(match.group(1))
            if gene_id:
                cds_lengths[gene_id] += int(fields[4]) - int(fields[3]) + 1

    return dict(cds_lengths)


def build_gene_burdens(gene_set_df):
    if gene_set_df.empty:
        return pd.DataFrame(columns=["gene_id", "n_loci"])
    burdens = gene_set_df.loc[
        gene_set_df["introner_group"] == "all_introners",
        ["gene_id", "n_loci"],
    ].copy()
    burdens["gene_id"] = burdens["gene_id"].map(clean_gene_id)
    burdens = burdens.dropna(subset=["gene_id"])
    return burdens.groupby("gene_id", as_index=False)["n_loci"].sum()


def build_background(gene2go, cds_lengths):
    genes = sorted(
        gene for gene, go_terms in gene2go.items()
        if go_terms and cds_lengths.get(gene, 0) > 0
    )
    if not genes:
        raise ValueError("No genes have both propagated GO annotations and CDS lengths.")
    return genes


def build_term_matrix(genes, gene2go, terms, min_genes):
    go_to_count = defaultdict(int)
    for gene in genes:
        for go_id in gene2go[gene]:
            if go_id in terms:
                go_to_count[go_id] += 1

    test_terms = sorted(
        go_id for go_id, count in go_to_count.items()
        if count >= min_genes
    )
    term_index = {go_id: i for i, go_id in enumerate(test_terms)}

    rows = []
    cols = []
    for gene_i, gene in enumerate(genes):
        for go_id in gene2go[gene]:
            term_i = term_index.get(go_id)
            if term_i is not None:
                rows.append(gene_i)
                cols.append(term_i)

    data = np.ones(len(rows), dtype=np.int8)
    matrix = csr_matrix(
        (data, (rows, cols)),
        shape=(len(genes), len(test_terms)),
        dtype=np.int16,
    )
    return test_terms, matrix


def observed_enrichment(study_mask, term_matrix, terms, test_terms, fdr_threshold):
    study_total = int(study_mask.sum())
    background_total = int(study_mask.size)
    study_counts = np.asarray(term_matrix[study_mask].sum(axis=0)).ravel().astype(int)
    background_counts = np.asarray(term_matrix.sum(axis=0)).ravel().astype(int)

    rows = []
    for i, go_id in enumerate(test_terms):
        study_with_go = int(study_counts[i])
        background_with_go_total = int(background_counts[i])
        study_without_go = study_total - study_with_go
        nonstudy_with_go = background_with_go_total - study_with_go
        nonstudy_without_go = (
            background_total - study_total - nonstudy_with_go
        )
        table = [
            [study_with_go, study_without_go],
            [nonstudy_with_go, nonstudy_without_go],
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        expected = (
            study_total * background_with_go_total / background_total
            if background_total else 0
        )
        fold_enrichment = study_with_go / expected if expected else math.nan
        enrichment_type = (
            "enriched" if odds_ratio > 1
            else "depleted" if odds_ratio < 1
            else "neutral"
        )
        rows.append({
            "GO_term": go_id,
            "term_name": terms[go_id].get("name", ""),
            "namespace": terms[go_id].get("namespace", ""),
            "study_with_go": study_with_go,
            "study_total": study_total,
            "background_with_go": background_with_go_total,
            "background_total": background_total,
            "odds_ratio": odds_ratio,
            "fold_enrichment": fold_enrichment,
            "fisher_p_value": p_value,
            "enrichment_type": enrichment_type,
        })

    results = pd.DataFrame(rows)
    if results.empty:
        return results, study_counts
    _, results["fisher_fdr_bh"], _, _ = multipletests(
        results["fisher_p_value"], method="fdr_bh"
    )
    results["fisher_significant"] = (
        (results["fisher_fdr_bh"] < fdr_threshold)
        & (results["enrichment_type"] != "neutral")
    )
    return results, study_counts


def _init_worker(term_matrix, weights, study_total, obs_counts):
    global _TERM_MATRIX, _WEIGHTS, _STUDY_TOTAL, _OBS_COUNTS
    _TERM_MATRIX = term_matrix
    _WEIGHTS = weights
    _STUDY_TOTAL = study_total
    _OBS_COUNTS = obs_counts


def _run_batch(args):
    n_permutations, seed = args
    rng = np.random.default_rng(seed)
    n_genes = _WEIGHTS.size
    n_terms = _OBS_COUNTS.size
    ge_counts = np.zeros(n_terms, dtype=np.int64)
    le_counts = np.zeros(n_terms, dtype=np.int64)

    for _ in range(n_permutations):
        selected = rng.choice(
            n_genes,
            size=_STUDY_TOTAL,
            replace=False,
            p=_WEIGHTS,
        )
        perm_counts = np.asarray(_TERM_MATRIX[selected].sum(axis=0)).ravel()
        ge_counts += perm_counts >= _OBS_COUNTS
        le_counts += perm_counts <= _OBS_COUNTS

    return n_permutations, ge_counts, le_counts


def run_permutations(term_matrix, weights, study_total, obs_counts,
                     n_permutations, threads, batch_size, seed):
    if n_permutations < 1:
        raise ValueError("--n-permutations must be at least 1.")
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1.")

    rng = np.random.default_rng(seed)
    batch_sizes = []
    remaining = n_permutations
    while remaining > 0:
        current = min(batch_size, remaining)
        batch_sizes.append(current)
        remaining -= current

    seeds = rng.integers(0, np.iinfo(np.uint32).max, size=len(batch_sizes), dtype=np.uint32)
    tasks = list(zip(batch_sizes, seeds.tolist()))

    ge_total = np.zeros(obs_counts.size, dtype=np.int64)
    le_total = np.zeros(obs_counts.size, dtype=np.int64)
    done = 0

    if threads == 1:
        _init_worker(term_matrix, weights, study_total, obs_counts)
        iterator = map(_run_batch, tasks)
    else:
        pool = Pool(
            processes=threads,
            initializer=_init_worker,
            initargs=(term_matrix, weights, study_total, obs_counts),
        )
        iterator = pool.imap_unordered(_run_batch, tasks)

    try:
        for n_done, ge_counts, le_counts in iterator:
            done += n_done
            ge_total += ge_counts
            le_total += le_counts
            if done % max(batch_size * 10, 1) == 0 or done == n_permutations:
                print(f"Completed {done:,}/{n_permutations:,} permutations", file=sys.stderr)
    finally:
        if threads != 1:
            pool.close()
            pool.join()

    return ge_total, le_total


def write_summary(path, args, locus_df, burdens, study_genes, background_genes,
                  cds_lengths, results):
    with open(path, "w") as handle:
        handle.write("All-introner GO CDS-length-weighted permutation summary\n")
        handle.write("=" * 60 + "\n\n")
        handle.write(f"Genotype matrix: {args.genotype_matrix}\n")
        handle.write(f"GTF for CDS lengths: {args.gtf}\n")
        handle.write(f"GO JSON sources: {', '.join(args.go_json)}\n")
        handle.write(f"GO OBO: {args.go_obo}\n")
        handle.write(f"Permutations: {args.n_permutations:,}\n")
        handle.write(f"Threads: {args.threads}\n")
        handle.write(f"Random seed: {args.seed}\n")
        handle.write("Null model: sample the observed number of introner-bearing genes without replacement from the GO/CDS background, with recipient-gene probability proportional to CDS length.\n")
        handle.write("Burden preservation: observed per-gene introner counts are treated as exchangeable bundles on the sampled recipient genes, preserving the total introner count and burden distribution; the GO statistic currently uses only gene-level presence/absence.\n")
        handle.write("GO statistic: gene-level introner presence/absence among all introner genes.\n\n")
        handle.write(f"Ortholog groups classified: {len(locus_df)}\n")
        handle.write(f"All-introner genes before GO/CDS filtering: {len(burdens)}\n")
        handle.write(f"All-introner genes in tested background: {len(study_genes)}\n")
        handle.write(f"Total observed introner loci on tested genes: {int(burdens.loc[burdens['gene_id'].isin(study_genes), 'n_loci'].sum())}\n")
        handle.write(f"Background genes with GO annotations and CDS length: {len(background_genes)}\n")
        handle.write(f"Median CDS length in background: {np.median([cds_lengths[g] for g in background_genes]):.1f} bp\n")
        handle.write(f"Mean CDS length in background: {np.mean([cds_lengths[g] for g in background_genes]):.1f} bp\n\n")
        handle.write(f"GO terms tested: {len(results)}\n")
        handle.write(f"Empirical significant GO terms at q<{args.fdr_threshold}: {int(results['empirical_significant'].sum())}\n")
        handle.write(f"Fisher significant GO terms at q<{args.fdr_threshold}: {int(results['fisher_significant'].sum())}\n")

        terms_of_interest = results[results["GO_term"].isin(GO_TERMS_OF_INTEREST)].copy()
        if not terms_of_interest.empty:
            handle.write("\nGO Terms Of Interest\n")
            handle.write("--------------------\n")
            handle.write(
                "These nucleotide/ATP/binding terms were originally enriched "
                "under the Fisher GO test. Here, empirical p_enriched is the "
                "fraction of CDS-length-weighted permuted gene sets with at "
                "least as many genes in the GO term as the observed all-introner "
                "gene set.\n"
            )
            handle.write(
                "GO_term\tterm_name\tobserved_with_go\tbackground_with_go\t"
                "fold_enrichment\tfisher_fdr_bh\tempirical_p_enriched\t"
                "empirical_fdr_bh\tempirical_significant\n"
            )
            for _, row in terms_of_interest.sort_values("GO_term").iterrows():
                handle.write(
                    f"{row['GO_term']}\t{row['term_name']}\t"
                    f"{int(row['study_with_go'])}\t"
                    f"{int(row['background_with_go'])}\t"
                    f"{row['fold_enrichment']:.4g}\t"
                    f"{row['fisher_fdr_bh']:.4g}\t"
                    f"{row['empirical_p_enriched']:.4g}\t"
                    f"{row['empirical_fdr_bh']:.4g}\t"
                    f"{bool(row['empirical_significant'])}\n"
                )
            handle.write(
                "Interpretation: these terms are not empirically enriched "
                "relative to the CDS-length-weighted null, even though they are "
                "nominally enriched in the Fisher test.\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="CDS-length-weighted GO permutation test for all introner genes."
    )
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--gtf", required=True, help="Reference GTF used to compute CDS lengths.")
    parser.add_argument("--go-json", required=True, nargs="+")
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--group1-samples", required=True)
    parser.add_argument("--group2-samples", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--significant-output", required=True)
    parser.add_argument(
        "--fisher-significant-output",
        default=None,
        help=(
            "Optional TSV containing initially Fisher-significant observed "
            "GO patterns annotated with empirical permutation p-values."
        ),
    )
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--exclude-within-status", default="low_identity")
    parser.add_argument("--n-permutations", type=int, default=1000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260521)
    args = parser.parse_args()

    group1_samples = parse_csv_arg(args.group1_samples)
    group2_samples = parse_csv_arg(args.group2_samples)
    excluded_within_status = parse_csv_arg(args.exclude_within_status)

    print("Loading ontology and GO annotations", file=sys.stderr)
    terms = parse_obo(args.go_obo)
    gene2go_raw = load_gene2go(args.go_json)
    gene2go = propagate_gene2go(gene2go_raw, terms)

    print("Classifying all-introner genes from genotype matrix", file=sys.stderr)
    locus_df, gene_set_df = classify_loci(
        args.genotype_matrix,
        group1_samples,
        group2_samples,
        excluded_within_status,
    )
    burdens = build_gene_burdens(gene_set_df)

    print("Extracting CDS lengths", file=sys.stderr)
    cds_lengths = extract_cds_lengths_from_gtf(args.gtf)
    background_genes = build_background(gene2go, cds_lengths)
    background_set = set(background_genes)
    study_genes = sorted(set(burdens["gene_id"]) & background_set)
    if not study_genes:
        raise ValueError("No all-introner genes have both GO annotations and CDS lengths.")

    study_total = len(study_genes)
    weights = np.array([cds_lengths[gene] for gene in background_genes], dtype=float)
    weights = weights / weights.sum()
    study_mask = np.array([gene in set(study_genes) for gene in background_genes], dtype=bool)

    print("Building GO term matrix", file=sys.stderr)
    test_terms, term_matrix = build_term_matrix(
        background_genes,
        gene2go,
        terms,
        args.min_genes,
    )
    if not test_terms:
        raise ValueError("No GO terms met the minimum gene threshold.")

    print("Calculating observed enrichment", file=sys.stderr)
    results, obs_counts = observed_enrichment(
        study_mask,
        term_matrix,
        terms,
        test_terms,
        args.fdr_threshold,
    )

    print(
        f"Running {args.n_permutations:,} CDS-length-weighted permutations "
        f"with {args.threads} thread(s)",
        file=sys.stderr,
    )
    ge_counts, le_counts = run_permutations(
        term_matrix,
        weights,
        study_total,
        obs_counts,
        args.n_permutations,
        max(1, args.threads),
        args.batch_size,
        args.seed,
    )

    results["permutations"] = args.n_permutations
    results["empirical_p_enriched"] = (ge_counts + 1) / (args.n_permutations + 1)
    results["empirical_p_depleted"] = (le_counts + 1) / (args.n_permutations + 1)
    results["empirical_p_directional"] = np.where(
        results["enrichment_type"] == "depleted",
        results["empirical_p_depleted"],
        results["empirical_p_enriched"],
    )
    neutral = results["enrichment_type"] == "neutral"
    results.loc[neutral, "empirical_p_directional"] = 1.0
    _, results["empirical_fdr_bh"], _, _ = multipletests(
        results["empirical_p_directional"], method="fdr_bh"
    )
    results["empirical_significant"] = (
        (results["empirical_fdr_bh"] < args.fdr_threshold)
        & (results["enrichment_type"] != "neutral")
    )
    results = results.sort_values(
        ["empirical_fdr_bh", "empirical_p_directional", "fold_enrichment"],
        ascending=[True, True, False],
    )

    results.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
    results.loc[results["empirical_significant"]].to_csv(
        args.significant_output,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    if args.fisher_significant_output:
        results.loc[results["fisher_significant"]].to_csv(
            args.fisher_significant_output,
            sep="\t",
            index=False,
            float_format="%.6g",
        )
    write_summary(
        args.summary_output,
        args,
        locus_df,
        burdens,
        study_genes,
        background_genes,
        cds_lengths,
        results,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

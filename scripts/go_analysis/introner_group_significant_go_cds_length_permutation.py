#!/usr/bin/env python3
"""
CDS-length-weighted permutation test for significant introner-group GO patterns.

This script starts from introner_group_go_enrichment.significant.tsv and asks,
for each observed significant group-term enrichment/depletion pattern, how often
a CDS-length-weighted random gene set of the same size produces an equally
extreme count in the same direction.
"""

import argparse
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from statsmodels.stats.multitest import multipletests

from introner_group_go_enrichment import (
    clean_gene_id,
    load_gene2go,
    parse_obo,
    propagate_gene2go,
)


_TERM_MATRIX = None
_WEIGHTS = None
_STUDY_TOTAL = None
_OBS_COUNTS = None


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


def build_background(gene2go, cds_lengths):
    genes = sorted(
        gene for gene, go_terms in gene2go.items()
        if go_terms and cds_lengths.get(gene, 0) > 0
    )
    if not genes:
        raise ValueError("No genes have both propagated GO annotations and CDS lengths.")
    return genes


def load_gene_sets(path, background_genes):
    background = set(background_genes)
    gene_sets = {}
    gene_set_df = pd.read_csv(path, sep="\t")
    for group, group_df in gene_set_df.groupby("introner_group", sort=False):
        genes = {
            gene for gene in group_df["gene_id"].map(clean_gene_id)
            if gene and gene in background
        }
        gene_sets[group] = genes
    return gene_sets


def build_term_matrix(genes, gene2go, selected_terms):
    term_index = {go_id: i for i, go_id in enumerate(selected_terms)}
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
        shape=(len(genes), len(selected_terms)),
        dtype=np.int16,
    )
    return matrix


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


def run_group_permutations(term_matrix, weights, study_total, obs_counts,
                           n_permutations, threads, batch_size, seed):
    rng = np.random.default_rng(seed)
    batch_sizes = []
    remaining = n_permutations
    while remaining > 0:
        current = min(batch_size, remaining)
        batch_sizes.append(current)
        remaining -= current

    seeds = rng.integers(
        0, np.iinfo(np.uint32).max, size=len(batch_sizes), dtype=np.uint32
    )
    tasks = list(zip(batch_sizes, seeds.tolist()))

    ge_total = np.zeros(obs_counts.size, dtype=np.int64)
    le_total = np.zeros(obs_counts.size, dtype=np.int64)
    done = 0

    if threads == 1:
        _init_worker(term_matrix, weights, study_total, obs_counts)
        iterator = map(_run_batch, tasks)
        pool = None
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
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    return ge_total, le_total


def write_summary(path, args, observed, output, background_genes, gene_sets):
    with open(path, "w") as handle:
        handle.write("Introner Group Significant GO CDS-length-weighted Permutation Summary\n")
        handle.write("=" * 75 + "\n\n")
        handle.write(f"Observed significant GO input: {args.observed_significant}\n")
        handle.write(f"Gene sets: {args.gene_sets}\n")
        handle.write(f"GTF for CDS lengths: {args.gtf}\n")
        handle.write(f"GO JSON sources: {', '.join(args.go_json)}\n")
        handle.write(f"Permutations per group: {args.n_permutations:,}\n")
        handle.write(f"Threads: {args.threads}\n")
        handle.write(f"Random seed: {args.seed}\n")
        handle.write("Null model: for each introner group, sample the observed number of GO/CDS-eligible study genes without replacement from the GO/CDS background, with recipient-gene probability proportional to CDS length.\n")
        handle.write("Empirical p-value: enriched terms use P(permuted_count >= observed_count); depleted terms use P(permuted_count <= observed_count), each with +1 correction.\n\n")
        handle.write(f"Background genes with GO annotations and CDS length: {len(background_genes)}\n")
        handle.write(f"Observed significant group-term rows tested: {len(observed)}\n")
        handle.write(f"Empirical significant rows at global q<{args.fdr_threshold}: {int(output['empirical_significant_global'].sum())}\n\n")

        for group in sorted(observed["introner_group"].unique()):
            group_rows = output[output["introner_group"] == group]
            empirical_sig = int(group_rows["empirical_significant_global"].sum())
            enriched = int(((group_rows["empirical_significant_global"]) & (group_rows["enrichment_type"] == "enriched")).sum())
            depleted = int(((group_rows["empirical_significant_global"]) & (group_rows["enrichment_type"] == "depleted")).sum())
            handle.write(
                f"{group}: {len(gene_sets.get(group, set()))} GO/CDS-eligible genes, "
                f"{len(group_rows)} observed significant GO rows tested, "
                f"{empirical_sig} empirical significant rows "
                f"({enriched} enriched, {depleted} depleted)\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Permute significant introner-group GO enrichment patterns."
    )
    parser.add_argument("--observed-significant", required=True)
    parser.add_argument("--gene-sets", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--go-json", required=True, nargs="+")
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--empirical-significant-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--n-permutations", type=int, default=1000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260521)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    args = parser.parse_args()

    if args.n_permutations < 1:
        raise ValueError("--n-permutations must be at least 1.")

    print("Loading observed significant GO patterns", file=sys.stderr)
    observed = pd.read_csv(args.observed_significant, sep="\t")
    if observed.empty:
        raise ValueError("Observed significant GO table is empty.")
    required = {"introner_group", "GO_term", "study_with_go", "enrichment_type"}
    missing = sorted(required - set(observed.columns))
    if missing:
        raise ValueError("Observed table missing columns: " + ", ".join(missing))

    print("Loading ontology and GO annotations", file=sys.stderr)
    terms = parse_obo(args.go_obo)
    gene2go = propagate_gene2go(load_gene2go(args.go_json), terms)
    cds_lengths = extract_cds_lengths_from_gtf(args.gtf)
    background_genes = build_background(gene2go, cds_lengths)
    weights = np.array([cds_lengths[gene] for gene in background_genes], dtype=float)
    weights = weights / weights.sum()
    gene_sets = load_gene_sets(args.gene_sets, background_genes)

    outputs = []
    master_rng = np.random.default_rng(args.seed)
    for group, group_rows in observed.groupby("introner_group", sort=False):
        study_genes = gene_sets.get(group, set())
        if not study_genes:
            raise ValueError(f"No GO/CDS-eligible study genes found for {group}.")
        selected_terms = list(group_rows["GO_term"])
        obs_counts = group_rows["study_with_go"].astype(int).to_numpy()
        term_matrix = build_term_matrix(background_genes, gene2go, selected_terms)

        group_seed = int(master_rng.integers(0, np.iinfo(np.uint32).max))
        print(
            f"Running {args.n_permutations:,} permutations for {group} "
            f"({len(study_genes)} genes, {len(selected_terms)} GO rows)",
            file=sys.stderr,
        )
        ge_counts, le_counts = run_group_permutations(
            term_matrix,
            weights,
            len(study_genes),
            obs_counts,
            args.n_permutations,
            max(1, args.threads),
            args.batch_size,
            group_seed,
        )

        group_output = group_rows.copy()
        group_output["permutations"] = args.n_permutations
        group_output["study_total_permutation"] = len(study_genes)
        group_output["background_total_permutation"] = len(background_genes)
        group_output["empirical_p_enriched"] = (ge_counts + 1) / (args.n_permutations + 1)
        group_output["empirical_p_depleted"] = (le_counts + 1) / (args.n_permutations + 1)
        group_output["empirical_p_directional"] = np.where(
            group_output["enrichment_type"] == "depleted",
            group_output["empirical_p_depleted"],
            group_output["empirical_p_enriched"],
        )
        group_output.loc[
            ~group_output["enrichment_type"].isin(["enriched", "depleted"]),
            "empirical_p_directional",
        ] = 1.0
        outputs.append(group_output)

    output = pd.concat(outputs, ignore_index=True)
    _, output["empirical_fdr_bh_global"], _, _ = multipletests(
        output["empirical_p_directional"], method="fdr_bh"
    )
    output["empirical_fdr_bh_by_group"] = np.nan
    for group, idx in output.groupby("introner_group").groups.items():
        _, corrected, _, _ = multipletests(
            output.loc[idx, "empirical_p_directional"], method="fdr_bh"
        )
        output.loc[idx, "empirical_fdr_bh_by_group"] = corrected
    output["empirical_significant_global"] = (
        output["empirical_fdr_bh_global"] < args.fdr_threshold
    )
    output["empirical_significant_by_group"] = (
        output["empirical_fdr_bh_by_group"] < args.fdr_threshold
    )

    output = output.sort_values(
        ["empirical_fdr_bh_global", "empirical_p_directional", "introner_group", "GO_term"]
    )
    output.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
    output.loc[output["empirical_significant_global"]].to_csv(
        args.empirical_significant_output,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    write_summary(args.summary_output, args, observed, output, background_genes, gene_sets)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

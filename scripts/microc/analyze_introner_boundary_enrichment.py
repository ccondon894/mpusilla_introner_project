#!/usr/bin/env python3
"""Test present-introner insulation scores against gene-conditioned permutations."""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PRESENT = 1


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--insulation", required=True)
    parser.add_argument("--genes-gtf", required=True)
    parser.add_argument("--exclude-chroms", default="")
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=1545)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--chrom-summary-tsv", required=True)
    parser.add_argument("--loci-tsv", required=True)
    parser.add_argument("--null-tsv", required=True)
    parser.add_argument("--png", required=True)
    parser.add_argument("--pdf", required=True)
    return parser.parse_args()


def parse_excluded_chroms(value):
    return {chrom for chrom in value.split(",") if chrom}


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def exclude_chroms(chrom_sizes, excluded_chroms):
    if not excluded_chroms:
        return chrom_sizes
    return {
        chrom: size
        for chrom, size in chrom_sizes.items()
        if chrom not in excluded_chroms
    }


def parse_gtf_attributes(attributes):
    parsed = {}
    for item in attributes.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            parsed[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value.strip().strip('"')
    return parsed


def load_gene_intervals(path, chrom_sizes):
    gene_rows = []
    transcript_rows = []
    feature_rows_by_gene = {}

    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            feature = fields[2]
            start = max(0, int(fields[3]) - 1)
            end = min(chrom_sizes[chrom], int(fields[4]))
            if end <= start:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("ID")
            if not gene_id:
                gene_id = attrs.get("transcript_id", f"{chrom}:{start}-{end}")
            row = {"chrom": chrom, "start": start, "end": end, "gene_id": gene_id}
            if feature == "gene":
                gene_rows.append(row)
            elif feature == "transcript":
                transcript_rows.append(row)
            elif feature in {"exon", "CDS"}:
                feature_rows_by_gene.setdefault(gene_id, []).append(row)

    if gene_rows:
        genes = pd.DataFrame(gene_rows).drop_duplicates("gene_id")
    elif transcript_rows:
        genes = pd.DataFrame(transcript_rows).drop_duplicates("gene_id")
    else:
        rows = []
        for gene_id, records in feature_rows_by_gene.items():
            chroms = {record["chrom"] for record in records}
            if len(chroms) != 1:
                continue
            rows.append(
                {
                    "chrom": records[0]["chrom"],
                    "start": min(record["start"] for record in records),
                    "end": max(record["end"] for record in records),
                    "gene_id": gene_id,
                }
            )
        genes = pd.DataFrame(rows)
    return genes[["chrom", "start", "end", "gene_id"]]


def merge_intervals(intervals):
    by_chrom = {}
    for row in intervals.itertuples(index=False):
        by_chrom.setdefault(row.chrom, []).append((int(row.start), int(row.end)))

    merged_by_chrom = {}
    total_intervals = 0
    total_bp = 0
    for chrom, chrom_intervals in by_chrom.items():
        chrom_intervals.sort()
        merged = []
        for start, end in chrom_intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        merged_by_chrom[chrom] = [(start, end) for start, end in merged]
        total_intervals += len(merged)
        total_bp += sum(end - start for start, end in merged)
    return merged_by_chrom, total_intervals, total_bp


def load_present_introners(path, sample, chrom_sizes):
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["sample"].astype(str) == sample].copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["presence"] = df["presence"].astype(int)
    df = df.loc[df["contig"].isin(chrom_sizes)].copy()
    df = df.loc[df["presence"] == PRESENT].copy()
    df["length"] = df["end"] - df["start"]
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype(int)
    return df[
        [
            "ortholog_id",
            "contig",
            "start",
            "end",
            "length",
            "midpoint",
            "family",
            "gene",
        ]
    ]


def build_shuffle_mask(genes, loci, chrom_sizes):
    rows = genes[["chrom", "start", "end"]].copy()
    observed_rows = pd.DataFrame(
        {
            "chrom": loci["contig"],
            "start": loci["start"].clip(lower=0),
            "end": [
                min(chrom_sizes[chrom], int(end))
                for chrom, end in zip(loci["contig"], loci["end"])
            ],
        }
    )
    rows = pd.concat([rows, observed_rows], ignore_index=True)
    rows = rows.loc[rows["end"] > rows["start"]].copy()
    return merge_intervals(rows)


def insulation_score_column(insulation, window):
    exact = f"log2_insulation_score_{window}"
    if exact in insulation.columns:
        return exact
    matches = [
        column for column in insulation.columns
        if column.startswith("log2_insulation_score")
    ]
    if not matches:
        raise ValueError("Could not find a log2_insulation_score column")
    return matches[0]


def build_insulation_lookup(path, chrom_sizes, window):
    insulation = pd.read_csv(path, sep="\t")
    score_col = insulation_score_column(insulation, window)
    lookup = {}
    for chrom, data in insulation.groupby("chrom", sort=False):
        if chrom not in chrom_sizes:
            continue
        data = data.sort_values("start")
        lookup[chrom] = {
            "starts": data["start"].to_numpy(dtype=int),
            "ends": data["end"].to_numpy(dtype=int),
            "scores": data[score_col].to_numpy(dtype=float),
        }
    return lookup, score_col


def scores_at_midpoints(chroms, midpoints, insulation_lookup):
    scores = np.full(len(midpoints), np.nan, dtype=float)
    for idx, (chrom, midpoint) in enumerate(zip(chroms, midpoints)):
        chrom_lookup = insulation_lookup.get(chrom)
        if chrom_lookup is None:
            continue
        bin_idx = np.searchsorted(chrom_lookup["starts"], int(midpoint), side="right") - 1
        if bin_idx < 0 or bin_idx >= len(chrom_lookup["starts"]):
            continue
        if int(midpoint) >= chrom_lookup["ends"][bin_idx]:
            continue
        scores[idx] = chrom_lookup["scores"][bin_idx]
    return scores


def prepare_shuffle_choices(loci, chrom_sizes, shuffle_mask):
    choices = []
    for row in loci.itertuples(index=False):
        chrom_size = chrom_sizes[row.contig]
        length = max(1, min(int(row.length), chrom_size))
        eligible = [
            (start, end)
            for start, end in shuffle_mask.get(row.contig, [])
            if end - start >= length
        ]
        if eligible:
            weights = np.array([end - start - length + 1 for start, end in eligible])
            weights = weights / weights.sum()
        else:
            weights = None
        choices.append(
            {
                "chrom": row.contig,
                "length": length,
                "eligible": eligible,
                "weights": weights,
                "chrom_size": chrom_size,
            }
        )
    return choices


def shuffled_midpoints(shuffle_choices, rng):
    chroms = []
    midpoints = []
    for choice_data in shuffle_choices:
        length = choice_data["length"]
        eligible = choice_data["eligible"]
        if eligible:
            choice = int(rng.choice(len(eligible), p=choice_data["weights"]))
            interval_start, interval_end = eligible[choice]
            max_start = interval_end - length
            start = int(rng.integers(interval_start, max_start + 1))
        else:
            max_start = max(0, choice_data["chrom_size"] - length)
            start = int(rng.integers(0, max_start + 1)) if max_start else 0
        chroms.append(choice_data["chrom"])
        midpoints.append(start + length // 2)
    return chroms, np.array(midpoints, dtype=int)


def chromosome_number(chrom):
    match = chrom.rsplit("scaffold_", 1)
    if len(match) == 2 and match[1].isdigit():
        return int(match[1])
    fields = chrom.split("#")
    for field in reversed(fields):
        if field.isdigit():
            return int(field)
    return None


def chromosome_sort_key(chrom):
    number = chromosome_number(chrom)
    return (number is None, number if number is not None else chrom)


def permutation_test(loci, chrom_sizes, shuffle_mask, insulation_lookup, permutations, seed):
    rng = np.random.default_rng(seed)
    choices = prepare_shuffle_choices(loci, chrom_sizes, shuffle_mask)
    null_means = np.full(permutations, np.nan, dtype=float)
    chroms = sorted(loci["contig"].unique(), key=chromosome_sort_key)
    chrom_to_idx = {chrom: idx for idx, chrom in enumerate(chroms)}
    chrom_null_means = np.full((permutations, len(chroms)), np.nan, dtype=float)

    for idx in range(permutations):
        shuffled_chroms, midpoints = shuffled_midpoints(choices, rng)
        scores = scores_at_midpoints(shuffled_chroms, midpoints, insulation_lookup)
        null_means[idx] = np.nanmean(scores)
        shuffled_chroms = np.array(shuffled_chroms, dtype=object)
        for chrom in chroms:
            chrom_scores = scores[shuffled_chroms == chrom]
            chrom_null_means[idx, chrom_to_idx[chrom]] = np.nanmean(chrom_scores)
    return null_means, chroms, chrom_null_means


def summarize_chromosomes(loci, chroms, chrom_null_means, sample, permutations):
    rows = []
    for idx, chrom in enumerate(chroms):
        chrom_loci = loci.loc[loci["contig"] == chrom].copy()
        observed_scores = chrom_loci["insulation_score"].to_numpy(dtype=float)
        observed_scores = observed_scores[np.isfinite(observed_scores)]
        observed = float(observed_scores.mean()) if len(observed_scores) else np.nan
        null_values = chrom_null_means[:, idx]
        null_values = null_values[np.isfinite(null_values)]
        expected = float(null_values.mean()) if len(null_values) else np.nan
        null_sd = float(null_values.std(ddof=1)) if len(null_values) > 1 else np.nan
        rows.append(
            {
                "sample": sample,
                "chrom": chrom,
                "state": "present",
                "statistic": "mean_insulation_score_at_introner_midpoints",
                "n_loci": len(chrom_loci),
                "n_loci_scored": len(observed_scores),
                "observed": observed,
                "expected": expected,
                "delta": observed - expected,
                "linear_observed": 2 ** observed if np.isfinite(observed) else np.nan,
                "linear_expected": 2 ** expected if np.isfinite(expected) else np.nan,
                "linear_percent_change": (
                    ((2 ** observed) - (2 ** expected)) / (2 ** expected) * 100
                    if np.isfinite(observed) and np.isfinite(expected)
                    else np.nan
                ),
                "z_score": (
                    (observed - expected) / null_sd
                    if null_sd and np.isfinite(observed)
                    else np.nan
                ),
                "empirical_p_upper": (
                    (np.sum(null_values >= observed) + 1) / (len(null_values) + 1)
                    if len(null_values) and np.isfinite(observed)
                    else np.nan
                ),
                "empirical_p_lower": (
                    (np.sum(null_values <= observed) + 1) / (len(null_values) + 1)
                    if len(null_values) and np.isfinite(observed)
                    else np.nan
                ),
                "permutations": permutations,
            }
        )
    return pd.DataFrame(rows)


def plot_null(null_means, observed_mean, sample, png, pdf):
    finite_null = null_means[np.isfinite(null_means)]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(finite_null, bins=50, color="#8d99ae", edgecolor="white", linewidth=0.4)
    ax.axvline(observed_mean, color="#2d6cdf", linewidth=2.2, label="Observed")
    ax.axvline(
        finite_null.mean(),
        color="#293241",
        linewidth=1.8,
        linestyle="--",
        label="Permutation mean",
    )
    ax.set_xlabel("Mean insulation score at present-introner midpoints")
    ax.set_ylabel("Permutations")
    ax.set_title(f"{sample} gene-conditioned introner insulation")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


def main():
    args = parse_args()
    for path in (
        args.summary_tsv,
        args.chrom_summary_tsv,
        args.loci_tsv,
        args.null_tsv,
        args.png,
        args.pdf,
    ):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    excluded_chroms = parse_excluded_chroms(args.exclude_chroms)
    chrom_sizes = exclude_chroms(read_chrom_sizes(args.chrom_sizes), excluded_chroms)
    loci = load_present_introners(args.genotype_matrix, args.sample, chrom_sizes)
    genes = load_gene_intervals(args.genes_gtf, chrom_sizes)
    shuffle_mask, shuffle_mask_intervals, shuffle_mask_bp = build_shuffle_mask(
        genes,
        loci,
        chrom_sizes,
    )
    insulation_lookup, score_col = build_insulation_lookup(
        args.insulation,
        chrom_sizes,
        args.window,
    )

    observed_scores = scores_at_midpoints(
        loci["contig"].to_list(),
        loci["midpoint"].to_numpy(dtype=int),
        insulation_lookup,
    )
    loci = loci.copy()
    loci["insulation_score"] = observed_scores
    scored_loci = loci.loc[np.isfinite(loci["insulation_score"])].copy()
    observed_mean = float(scored_loci["insulation_score"].mean())

    null_means, chroms, chrom_null_means = permutation_test(
        loci,
        chrom_sizes,
        shuffle_mask,
        insulation_lookup,
        args.permutations,
        args.seed,
    )
    finite_null = null_means[np.isfinite(null_means)]
    expected = float(finite_null.mean()) if len(finite_null) else np.nan
    null_sd = float(finite_null.std(ddof=1)) if len(finite_null) > 1 else np.nan

    summary = pd.DataFrame(
        [
            {
                "sample": args.sample,
                "state": "present",
                "statistic": "mean_insulation_score_at_introner_midpoints",
                "score_column": score_col,
                "n_loci": len(loci),
                "n_loci_scored": len(scored_loci),
                "observed": observed_mean,
                "expected": expected,
                "z_score": (
                    (observed_mean - expected) / null_sd
                    if null_sd and np.isfinite(observed_mean)
                    else np.nan
                ),
                "empirical_p_upper": (
                    (np.sum(finite_null >= observed_mean) + 1) / (len(finite_null) + 1)
                    if len(finite_null)
                    else np.nan
                ),
                "empirical_p_lower": (
                    (np.sum(finite_null <= observed_mean) + 1) / (len(finite_null) + 1)
                    if len(finite_null)
                    else np.nan
                ),
                "permutations": args.permutations,
                "null_model": "same_chrom_length_preserving_shuffle_within_gene_bodies",
                "excluded_chroms": ",".join(sorted(excluded_chroms)) or ".",
                "shuffle_mask_intervals": shuffle_mask_intervals,
                "shuffle_mask_bp": shuffle_mask_bp,
                "n_genes_in_mask": len(genes),
            }
        ]
    )

    null = pd.DataFrame(
        {
            "sample": args.sample,
            "permutation": np.arange(1, len(null_means) + 1),
            "mean_insulation_score": null_means,
        }
    )
    chrom_summary = summarize_chromosomes(
        loci,
        chroms,
        chrom_null_means,
        args.sample,
        args.permutations,
    )
    chrom_summary.insert(4, "score_column", score_col)
    chrom_summary["null_model"] = "same_chrom_length_preserving_shuffle_within_gene_bodies"
    chrom_summary["excluded_chroms"] = ",".join(sorted(excluded_chroms)) or "."

    summary.to_csv(args.summary_tsv, sep="\t", index=False)
    chrom_summary.to_csv(args.chrom_summary_tsv, sep="\t", index=False)
    loci.to_csv(args.loci_tsv, sep="\t", index=False)
    null.to_csv(args.null_tsv, sep="\t", index=False)
    plot_null(null_means, observed_mean, args.sample, args.png, args.pdf)


if __name__ == "__main__":
    main()

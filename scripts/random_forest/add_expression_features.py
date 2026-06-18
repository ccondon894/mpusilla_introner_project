#!/usr/bin/env python3
"""Add CCMP1545 expression and optional introner PSI features to an RF matrix."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument(
        "--counts",
        default="/scratch1/chris/mpusilla_introner_project/results/expression/counts/merged_counts_matrix.csv",
    )
    parser.add_argument(
        "--boundary-support",
        default="/scratch1/chris/mpusilla_introner_project/results/expression/splice_junctions/introner_boundary_support/CCMP1545.per_locus.tsv",
    )
    parser.add_argument(
        "--sample-prefix",
        default="834",
        help="Column prefix for CCMP1545 RNA-seq replicates in the counts matrix.",
    )
    parser.add_argument(
        "--skip-psi",
        action="store_true",
        help="Only add gene-level TPM features; skip introner PSI/splice support.",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args()


def collapse_boundary_support(path):
    df = pd.read_csv(path, sep="\t")
    key = ["contig", "body_start", "body_end"]
    df = (
        df.sort_values(
            key + ["best_junction_score", "max_abs_boundary_delta", "sum_abs_boundary_delta"],
            ascending=[True, True, True, False, True, True],
        )
        .drop_duplicates(key, keep="first")
    )
    keep = key + [
        "retention_psi",
        "spliced_fraction",
        "covered_for_psi",
        "psi_junction_signal_within_2bp",
        "retained_signal",
        "retained_aligned_bases",
        "retained_mean_depth",
        "best_junction_score",
        "within_2bp_junction_support",
        "within_5bp_junction_support",
        "within_10bp_junction_support",
    ]
    out = df[keep].rename(columns={
        "retention_psi": "expr_splice_retention_psi",
        "spliced_fraction": "expr_splice_spliced_fraction",
        "covered_for_psi": "expr_splice_covered_for_psi",
        "psi_junction_signal_within_2bp": "expr_splice_junction_signal_within_2bp",
        "retained_signal": "expr_splice_retained_signal",
        "retained_aligned_bases": "expr_splice_retained_aligned_bases",
        "retained_mean_depth": "expr_splice_retained_mean_depth",
        "best_junction_score": "expr_splice_best_junction_score",
        "within_2bp_junction_support": "expr_splice_within_2bp_junction_support",
        "within_5bp_junction_support": "expr_splice_within_5bp_junction_support",
        "within_10bp_junction_support": "expr_splice_within_10bp_junction_support",
    })
    out["expr_splice_covered_for_psi"] = out["expr_splice_covered_for_psi"].map({
        "yes": 1,
        "no": 0,
    })
    return out


def calculate_tpm(counts_path, sample_prefix):
    counts = pd.read_csv(counts_path)
    sample_cols = [
        col for col in counts.columns
        if col.startswith(f"{sample_prefix}_")
    ]
    if not sample_cols:
        raise ValueError(f"No columns found for sample prefix {sample_prefix!r}")

    length_kb = counts["Length"].astype(float) / 1000.0
    tpm_df = counts[["Geneid"]].copy()
    tpm_values = []
    for col in sample_cols:
        rpk = counts[col].astype(float) / length_kb
        scale = rpk.sum() / 1_000_000.0
        tpm = rpk / scale if scale else rpk * np.nan
        tpm_df[f"expr_tpm_{col}"] = tpm
        tpm_values.append(tpm)

    tpm_matrix = pd.concat(tpm_values, axis=1)
    tpm_df["expr_tpm_mean"] = tpm_matrix.mean(axis=1)
    tpm_df["expr_tpm_median"] = tpm_matrix.median(axis=1)
    tpm_df["expr_tpm_sd"] = tpm_matrix.std(axis=1)
    tpm_df["expr_tpm_max"] = tpm_matrix.max(axis=1)
    tpm_df["expr_tpm_log1p_mean"] = np.log1p(tpm_df["expr_tpm_mean"])
    tpm_df["expr_tpm_cv"] = tpm_df["expr_tpm_sd"] / tpm_df["expr_tpm_mean"].replace(0, np.nan)
    return tpm_df


def summarize(df):
    expr_cols = [
        col for col in df.columns
        if col.startswith("expr_") and pd.api.types.is_numeric_dtype(df[col])
    ]
    rows = [
        {"metric": "n_rows", "value": len(df)},
        {"metric": "n_expr_numeric_features", "value": len(expr_cols)},
    ]
    for col in expr_cols:
        values = df[col]
        rows.append({
            "metric": f"{col}_missing_fraction",
            "value": float(values.isna().mean()),
        })
        for label, group in df.groupby("label", dropna=False):
            rows.append({
                "metric": f"{col}_mean_label_{label}",
                "value": group[col].mean(),
            })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    output = Path(args.output)
    summary = Path(args.summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)

    print("Loading matrix...")
    matrix_df = pd.read_csv(args.matrix, sep="\t", low_memory=False)

    out_df = matrix_df
    if args.skip_psi:
        print("Skipping introner PSI features.")
    else:
        print("Adding introner PSI features...")
        boundary_df = collapse_boundary_support(args.boundary_support)
        out_df = matrix_df.merge(
            boundary_df,
            left_on=["contig", "start", "end"],
            right_on=["contig", "body_start", "body_end"],
            how="left",
        ).drop(columns=["body_start", "body_end"], errors="ignore")

    print("Calculating TPM features...")
    tpm_df = calculate_tpm(args.counts, args.sample_prefix)
    out_df = out_df.merge(
        tpm_df,
        left_on="gene",
        right_on="Geneid",
        how="left",
    ).drop(columns=["Geneid"], errors="ignore")

    out_df.to_csv(output, sep="\t", index=False)
    summarize(out_df).to_csv(summary, sep="\t", index=False)
    print(f"Wrote {output}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()

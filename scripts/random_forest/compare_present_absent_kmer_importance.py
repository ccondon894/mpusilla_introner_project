#!/usr/bin/env python3
"""Compare CCMP1545 and RCC1749 present/absent k-mer permutation importance."""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ccmp1545-importance", required=True)
    parser.add_argument("--rcc1749-importance", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_importance(path, prefix):
    df = pd.read_csv(path, sep="\t")
    keep = [
        "kmer",
        "n_features",
        "roc_auc_importance_mean",
        "roc_auc_importance_ci95",
        "roc_auc_importance_rank",
        "pr_auc_importance_mean",
        "pr_auc_importance_ci95",
        "pr_auc_importance_rank",
        "baseline_roc_auc_mean",
        "baseline_pr_auc_mean",
    ]
    df = df[keep].copy()
    return df.rename(columns={col: f"{prefix}_{col}" for col in keep if col != "kmer"})


def main():
    args = parse_args()
    ccmp = read_importance(args.ccmp1545_importance, "ccmp1545")
    rcc = read_importance(args.rcc1749_importance, "rcc1749")
    merged = ccmp.merge(rcc, on="kmer", how="outer", validate="one_to_one")
    merged["delta_roc_auc_importance_mean"] = (
        merged["rcc1749_roc_auc_importance_mean"]
        - merged["ccmp1545_roc_auc_importance_mean"]
    )
    merged["delta_pr_auc_importance_mean"] = (
        merged["rcc1749_pr_auc_importance_mean"]
        - merged["ccmp1545_pr_auc_importance_mean"]
    )
    merged["roc_auc_rank_shift_rcc1749_minus_ccmp1545"] = (
        merged["rcc1749_roc_auc_importance_rank"]
        - merged["ccmp1545_roc_auc_importance_rank"]
    )
    merged["pr_auc_rank_shift_rcc1749_minus_ccmp1545"] = (
        merged["rcc1749_pr_auc_importance_rank"]
        - merged["ccmp1545_pr_auc_importance_rank"]
    )
    merged = merged.sort_values(
        ["delta_roc_auc_importance_mean", "rcc1749_roc_auc_importance_mean"],
        ascending=[False, False],
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, sep="\t", index=False)
    print(merged.head(25).to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

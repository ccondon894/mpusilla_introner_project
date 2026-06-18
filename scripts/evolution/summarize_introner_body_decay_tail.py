#!/usr/bin/env python3
"""Summarize introner-body pi tails and frequency-level decay patterns."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METRICS = [
    "pi_introner_body",
    "segregating_sites_per_valid_bp",
    "gap_fraction",
    "max_low_concordance_fraction",
    "min_mean_major_allele_fraction",
    "max_softclip_read_fraction",
    "max_depth_ratio",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-locus", required=True, type=Path)
    parser.add_argument("--tail-loci", required=True, type=Path)
    parser.add_argument("--frequency-summary", required=True, type=Path)
    parser.add_argument("--tail-summary", required=True, type=Path)
    parser.add_argument("--plot", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--tail-quantile", type=float, default=0.95)
    return parser.parse_args()


def numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in METRICS + [
        "group_present_count",
        "n_aligned_sequences",
        "n_present_rows",
        "pairwise_comparisons",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def summarize_group(group: pd.DataFrame, value_col: str) -> pd.Series:
    values = group[value_col].dropna()
    out = {
        "n_loci": len(group),
        f"n_{value_col}": len(values),
        f"mean_{value_col}": values.mean() if len(values) else np.nan,
        f"median_{value_col}": values.median() if len(values) else np.nan,
        f"q75_{value_col}": values.quantile(0.75) if len(values) else np.nan,
        f"q90_{value_col}": values.quantile(0.90) if len(values) else np.nan,
        f"q95_{value_col}": values.quantile(0.95) if len(values) else np.nan,
        f"max_{value_col}": values.max() if len(values) else np.nan,
    }
    return pd.Series(out)


def tail_summary(df: pd.DataFrame, tail_quantile: float) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(["analysis_scope", "analysis_class", "family"], dropna=False):
        scope, cls, family = keys
        values = group["pi_introner_body"].dropna()
        threshold = values.quantile(tail_quantile) if len(values) else np.nan
        tail = group[group["pi_introner_body"] >= threshold] if len(values) else group.iloc[0:0]
        top5_n = max(1, int(np.ceil(0.05 * len(values)))) if len(values) else 0
        top5 = values.sort_values(ascending=False).head(top5_n)
        total = values.sum()
        rows.append(
            {
                "analysis_scope": scope,
                "analysis_class": cls,
                "family": family,
                "n_loci": len(group),
                "n_pi": len(values),
                "mean_pi": values.mean() if len(values) else np.nan,
                "median_pi": values.median() if len(values) else np.nan,
                "trimmed_mean_pi_drop_top5pct": (
                    values[~values.index.isin(top5.index)].mean()
                    if len(values) > len(top5)
                    else np.nan
                ),
                "tail_quantile": tail_quantile,
                "tail_threshold_pi": threshold,
                "n_tail_loci": len(tail),
                "tail_sum_pi_fraction": (
                    tail["pi_introner_body"].sum() / total if total > 0 else np.nan
                ),
                "top5pct_sum_pi_fraction": top5.sum() / total if total > 0 else np.nan,
                "median_tail_gap_fraction": tail["gap_fraction"].median()
                if len(tail)
                else np.nan,
                "median_tail_low_concordance_fraction": tail[
                    "max_low_concordance_fraction"
                ].median()
                if len(tail)
                else np.nan,
                "median_tail_n_aligned_sequences": tail["n_aligned_sequences"].median()
                if len(tail)
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def frequency_summary(df: pd.DataFrame) -> pd.DataFrame:
    group1 = df[df["analysis_scope"] == "group1_primary"].copy()
    group1["frequency_bin"] = group1["group_present_count"].astype("Int64").astype(str)
    rows = []
    for keys, group in group1.groupby(["analysis_class", "family", "frequency_bin"], dropna=False):
        cls, family, frequency_bin = keys
        row = {
            "analysis_class": cls,
            "family": family,
            "group1_present_count": frequency_bin,
            "n_loci": len(group),
            "n_qc_pass": int((group["qc_drop_reason"] == "PASS").sum()),
            "n_with_pi": int(group["pi_introner_body"].notna().sum()),
            "median_n_aligned_sequences": group["n_aligned_sequences"].median(),
        }
        for metric in METRICS:
            if metric in group.columns:
                vals = group[metric].dropna()
                row[f"mean_{metric}"] = vals.mean() if len(vals) else np.nan
                row[f"median_{metric}"] = vals.median() if len(vals) else np.nan
                row[f"q90_{metric}"] = vals.quantile(0.90) if len(vals) else np.nan
                row[f"max_{metric}"] = vals.max() if len(vals) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    out["group1_present_count_sort"] = pd.to_numeric(
        out["group1_present_count"], errors="coerce"
    )
    return out.sort_values(["analysis_class", "family", "group1_present_count_sort"]).drop(
        columns=["group1_present_count_sort"]
    )


def top_tail_loci(df: pd.DataFrame, top_n: int, tail_quantile: float) -> pd.DataFrame:
    group1_poly = df[
        (df["analysis_scope"] == "group1_primary")
        & (df["analysis_class"] == "polymorphic")
        & df["pi_introner_body"].notna()
    ].copy()
    threshold = group1_poly["pi_introner_body"].quantile(tail_quantile)
    group1_poly["is_tail_by_quantile"] = group1_poly["pi_introner_body"] >= threshold
    cols = [
        "ortholog_id",
        "family",
        "group_present_count",
        "n_present_rows",
        "n_aligned_sequences",
        "pi_introner_body",
        "segregating_sites_per_valid_bp",
        "gap_fraction",
        "max_low_concordance_fraction",
        "min_mean_major_allele_fraction",
        "max_softclip_read_fraction",
        "max_depth_ratio",
        "mean_callable_fraction",
        "qc_drop_reason",
        "present_samples",
        "body_len_min",
        "body_len_max",
        "pairwise_comparisons",
        "pairwise_valid_sites",
        "alignment_path",
        "is_tail_by_quantile",
    ]
    cols = [col for col in cols if col in group1_poly.columns]
    return group1_poly.sort_values("pi_introner_body", ascending=False)[cols].head(top_n)


def make_plot(df: pd.DataFrame, out_path: Path) -> None:
    group1 = df[df["analysis_scope"] == "group1_primary"].copy()
    group1["frequency"] = group1["group_present_count"].astype(int)
    grouped = (
        group1.groupby(["analysis_class", "frequency"])
        .agg(
            n=("ortholog_id", "size"),
            n_pi=("pi_introner_body", lambda s: s.notna().sum()),
            median_pi=("pi_introner_body", "median"),
            mean_pi=("pi_introner_body", "mean"),
            q90_pi=("pi_introner_body", lambda s: s.quantile(0.90)),
            median_gap=("gap_fraction", "median"),
            median_low_concordance=("max_low_concordance_fraction", "median"),
        )
        .reset_index()
    )

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = {"fixed_present": "#4C78A8", "polymorphic": "#F58518"}
    for cls, sub in grouped.groupby("analysis_class"):
        sub = sub.sort_values("frequency")
        label = cls.replace("_", " ")
        axes[0, 0].plot(sub["frequency"], sub["median_pi"], marker="o", label=label, color=colors.get(cls))
        axes[0, 0].plot(sub["frequency"], sub["mean_pi"], marker="s", linestyle="--", color=colors.get(cls), alpha=0.8)
        axes[0, 1].plot(sub["frequency"], sub["q90_pi"], marker="o", label=label, color=colors.get(cls))
        axes[1, 0].plot(sub["frequency"], sub["median_gap"], marker="o", label=label, color=colors.get(cls))
        axes[1, 1].bar(
            sub["frequency"] + (-0.18 if cls == "fixed_present" else 0.18),
            sub["n_pi"],
            width=0.34,
            label=label,
            color=colors.get(cls),
            alpha=0.85,
        )

    axes[0, 0].set_title("Group 1 body pi by present count")
    axes[0, 0].set_ylabel("pi")
    axes[0, 0].text(0.02, 0.94, "solid=median, dashed=mean", transform=axes[0, 0].transAxes, fontsize=9)
    axes[0, 1].set_title("Group 1 upper tail pi")
    axes[0, 1].set_ylabel("90th percentile pi")
    axes[1, 0].set_title("Group 1 gap fraction")
    axes[1, 0].set_ylabel("median gap fraction")
    axes[1, 1].set_title("Loci with estimable pi")
    axes[1, 1].set_ylabel("n loci")
    for ax in axes.flat:
        ax.set_xlabel("Group 1 present count")
        ax.set_xticks(range(1, 12))
        ax.grid(axis="y", alpha=0.25)
    axes[0, 1].legend(frameon=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    df = numeric_columns(pd.read_csv(args.per_locus, sep="\t"))
    args.tail_loci.parent.mkdir(parents=True, exist_ok=True)
    args.frequency_summary.parent.mkdir(parents=True, exist_ok=True)
    args.tail_summary.parent.mkdir(parents=True, exist_ok=True)

    top_tail_loci(df, args.top_n, args.tail_quantile).to_csv(
        args.tail_loci, sep="\t", index=False
    )
    frequency_summary(df).to_csv(args.frequency_summary, sep="\t", index=False)
    tail_summary(df, args.tail_quantile).to_csv(args.tail_summary, sep="\t", index=False)
    make_plot(df, args.plot)

    print(f"Wrote tail loci to {args.tail_loci}")
    print(f"Wrote frequency summary to {args.frequency_summary}")
    print(f"Wrote tail summary to {args.tail_summary}")
    print(f"Wrote plot to {args.plot}")


if __name__ == "__main__":
    main()

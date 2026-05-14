#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


GROUP_ORDER = [
    "all_introners",
    "ancestral",
    "independent_insertion",
    "polymorphic_within_group1",
    "consistent_group1",
    "consistent_group2",
]
GROUP_LABELS = {
    "all_introners": "All introners",
    "ancestral": "Ancestral",
    "independent_insertion": "Independent insertion",
    "polymorphic_within_group1": "Polymorphic within Group 1",
    "consistent_group1": "Consistent Group 1",
    "consistent_group2": "Consistent Group 2",
}


def save_figure_with_png(fig, path):
    """Save a figure and emit a sibling PNG for PDF outputs."""
    path = Path(path)
    fig.savefig(path, bbox_inches="tight")
    if path.suffix.lower() == ".pdf":
        fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")


def load_results(path):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return df
    df = df[df["study_with_go"] > 0].copy()
    df["display_group"] = df["introner_group"].map(GROUP_LABELS).fillna(df["introner_group"])
    df["display_term"] = df.apply(
        lambda row: f"{row['term_name']} ({row['GO_term']})",
        axis=1,
    )
    df["neg_log10_fdr"] = -np.log10(np.clip(df["fdr_bh_global"], 1e-300, 1.0))
    return df


def write_empty_plot(path, title, message):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    save_figure_with_png(fig, path)
    plt.close(fig)


def select_terms_for_heatmap(df, max_terms):
    significant = df[df["significant"]].copy()
    source = significant if not significant.empty else df.copy()
    if source.empty:
        return []

    source = source.sort_values(["fdr_bh_global", "p_value", "fold_enrichment"])
    selected = []
    seen = set()
    for _, row in source.iterrows():
        go_term = row["GO_term"]
        if go_term in seen:
            continue
        seen.add(go_term)
        selected.append(go_term)
        if len(selected) >= max_terms:
            break
    return selected


def plot_heatmap(df, output, max_terms):
    terms = select_terms_for_heatmap(df, max_terms)
    if not terms:
        write_empty_plot(
            output,
            "Introner Group GO Enrichment",
            "No GO terms with study genes were available to plot.",
        )
        return

    plot_df = df[df["GO_term"].isin(terms)].copy()
    term_order = (
        plot_df.groupby("display_term")["fdr_bh_global"]
        .min()
        .sort_values()
        .index
        .tolist()
    )
    group_labels = [GROUP_LABELS[group] for group in GROUP_ORDER]

    heatmap_data = plot_df.pivot_table(
        index="display_term",
        columns="display_group",
        values="fold_enrichment",
        aggfunc="max",
    ).reindex(index=term_order, columns=group_labels)

    annot_data = plot_df.pivot_table(
        index="display_term",
        columns="display_group",
        values="fdr_bh_global",
        aggfunc="min",
    ).reindex(index=term_order, columns=group_labels)

    annot = annot_data.copy().astype(object)
    for idx in annot.index:
        for col in annot.columns:
            value = annot_data.loc[idx, col]
            if pd.isna(value):
                annot.loc[idx, col] = ""
            elif value < 0.001:
                annot.loc[idx, col] = "***"
            elif value < 0.01:
                annot.loc[idx, col] = "**"
            elif value < 0.05:
                annot.loc[idx, col] = "*"
            else:
                annot.loc[idx, col] = f"{value:.2g}"

    fig_height = max(6, 0.35 * len(heatmap_data))
    fig, ax = plt.subplots(figsize=(11, fig_height))
    sns.heatmap(
        heatmap_data.fillna(0),
        annot=annot,
        fmt="",
        cmap="viridis",
        linewidths=0.4,
        cbar_kws={"label": "Fold enrichment"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Introner Group GO Enrichment")
    fig.tight_layout()
    save_figure_with_png(fig, output)
    plt.close(fig)


def plot_top_terms(df, output, max_terms):
    if df.empty:
        write_empty_plot(
            output,
            "Top GO Terms by Introner Group",
            "No GO terms with study genes were available to plot.",
        )
        return

    rows = []
    for group in GROUP_ORDER:
        subset = df[df["introner_group"] == group].copy()
        if subset.empty:
            continue
        subset = subset.sort_values(["fdr_bh_global", "p_value", "fold_enrichment"]).head(max_terms)
        rows.append(subset)

    if not rows:
        write_empty_plot(
            output,
            "Top GO Terms by Introner Group",
            "No GO terms with study genes were available to plot.",
        )
        return

    plot_df = pd.concat(rows, ignore_index=True)
    plot_df["label"] = plot_df["display_term"].str.wrap(55)
    plot_df["score"] = -np.log10(np.clip(plot_df["fdr_bh_global"], 1e-300, 1.0))
    plot_df.loc[~np.isfinite(plot_df["score"]), "score"] = 0

    n_groups = plot_df["introner_group"].nunique()
    fig, axes = plt.subplots(n_groups, 1, figsize=(11, max(4, 3.2 * n_groups)), squeeze=False)

    for ax, group in zip(axes.flatten(), [g for g in GROUP_ORDER if g in set(plot_df["introner_group"])]):
        subset = plot_df[plot_df["introner_group"] == group].sort_values("score", ascending=True)
        colors = ["#d95f02" if sig else "#7570b3" for sig in subset["significant"]]
        ax.barh(subset["label"], subset["score"], color=colors)
        ax.set_title(GROUP_LABELS[group])
        ax.set_xlabel("-log10(global BH FDR)")
        ax.grid(axis="x", alpha=0.25)
        if len(subset) == 0:
            ax.text(0.5, 0.5, "No terms", transform=ax.transAxes, ha="center", va="center")

    fig.tight_layout()
    save_figure_with_png(fig, output)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot final-matrix introner group GO enrichment.")
    parser.add_argument("enrichment_tsv")
    parser.add_argument("--heatmap-output", required=True)
    parser.add_argument("--top-terms-output", required=True)
    parser.add_argument("--max-terms", type=int, default=25)
    args = parser.parse_args()

    heatmap_output = Path(args.heatmap_output)
    top_terms_output = Path(args.top_terms_output)
    heatmap_output.parent.mkdir(parents=True, exist_ok=True)
    top_terms_output.parent.mkdir(parents=True, exist_ok=True)
    df = load_results(args.enrichment_tsv)

    plot_heatmap(
        df,
        heatmap_output,
        max_terms=args.max_terms,
    )
    plot_top_terms(
        df,
        top_terms_output,
        max_terms=min(args.max_terms, 10),
    )


if __name__ == "__main__":
    main()

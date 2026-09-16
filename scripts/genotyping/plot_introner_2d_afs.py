#!/usr/bin/env python3
"""
Plot a 2D AFS for introner presence/absence genotypes.

This counts the introner-present allele (presence == 1) within each population.
By default, loci with missing calls, low-identity within-group classifications,
or independent cross-group origins are excluded.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


INDEPENDENT_STATUSES = {"independent", "likely_independent"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute and plot an introner-present 2D AFS."
    )
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--group1-samples", nargs="+", required=True)
    parser.add_argument("--group2-samples", nargs="+", required=True)
    parser.add_argument("--output-pdf", required=True)
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--fig-width", type=float, default=12.0,
                        help="Figure width in inches (default: 12).")
    parser.add_argument("--fig-height", type=float, default=3.0,
                        help="Figure height in inches (default: 3).")
    parser.add_argument("--axis-fontsize", type=float, default=12.0,
                        help="Axis-label font size in points (default: 12).")
    parser.add_argument("--tick-fontsize", type=float, default=10.0,
                        help="Axis-tick font size in points (default: 10).")
    parser.add_argument("--annotation-fontsize", type=float, default=10.0,
                        help="Heatmap cell-label font size in points (default: 10).")
    parser.add_argument("--colorbar-fontsize", type=float, default=10.0,
                        help="Colorbar label and tick font size in points (default: 10).")
    parser.add_argument(
        "--fill-heatmap-height",
        action="store_true",
        help="Allow rectangular cells so the heatmap fills the colorbar height.",
    )
    parser.add_argument(
        "--match-colorbar-height",
        action="store_true",
        help="Size the colorbar to the heatmap height while retaining square cells.",
    )
    parser.add_argument(
        "--drop-absent-monomorphic",
        action="store_true",
        help="Drop loci where present count is 0 in both groups.",
    )
    return parser.parse_args()


def first_non_null(values, default=""):
    values = values.dropna().astype(str).unique()
    if len(values) == 0:
        return default
    return values[0]


def compute_present_afs(df, group1_samples, group2_samples,
                        drop_absent_monomorphic=False):
    group1_set = set(group1_samples)
    group2_set = set(group2_samples)
    spectrum = np.zeros((len(group2_samples) + 1, len(group1_samples) + 1),
                        dtype=int)

    skipped = {
        "missing_or_incomplete": 0,
        "independent": 0,
        "low_identity": 0,
        "absent_monomorphic": 0,
    }
    used = 0
    records = []

    for ortholog_id, group in df.groupby("ortholog_id", sort=True):
        within_status = first_non_null(group.get("within_group_status", pd.Series(dtype=str)))
        if within_status == "low_identity":
            skipped["low_identity"] += 1
            continue

        cross_statuses = set(
            group.get("cross_group_status", pd.Series(dtype=str))
            .dropna()
            .astype(str)
        )
        if cross_statuses & INDEPENDENT_STATUSES:
            skipped["independent"] += 1
            continue

        group1_rows = group[group["sample"].isin(group1_set)]
        group2_rows = group[group["sample"].isin(group2_set)]
        group1_calls = group1_rows["presence"].astype(int).tolist()
        group2_calls = group2_rows["presence"].astype(int).tolist()

        complete = (
            len(group1_calls) == len(group1_samples)
            and len(group2_calls) == len(group2_samples)
            and 3 not in group1_calls
            and 3 not in group2_calls
        )
        if not complete:
            skipped["missing_or_incomplete"] += 1
            continue

        group1_present = group1_calls.count(1)
        group2_present = group2_calls.count(1)

        if drop_absent_monomorphic and group1_present == 0 and group2_present == 0:
            skipped["absent_monomorphic"] += 1
            continue

        spectrum[group2_present, group1_present] += 1
        used += 1
        records.append({
            "ortholog_id": ortholog_id,
            "group1_present_count": group1_present,
            "group2_present_count": group2_present,
            "cross_group_status": ",".join(sorted(cross_statuses)) if cross_statuses else "NA",
            "within_group_status": within_status if within_status else "NA",
        })

    return spectrum, pd.DataFrame(records), used, skipped


def write_spectrum_tsv(spectrum, output_tsv):
    rows = []
    for group2_present in range(spectrum.shape[0]):
        for group1_present in range(spectrum.shape[1]):
            rows.append({
                "group2_present_count": group2_present,
                "group1_present_count": group1_present,
                "count": int(spectrum[group2_present, group1_present]),
            })
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)


def write_summary(summary_path, spectrum, used, skipped):
    with open(summary_path, "w") as handle:
        handle.write("Introner 2D AFS summary\n")
        handle.write("=======================\n\n")
        handle.write("Allele counted: introner present call (presence == 1)\n")
        handle.write(f"Loci used: {used:,}\n")
        for reason, count in skipped.items():
            handle.write(f"Loci skipped ({reason}): {count:,}\n")
        handle.write("\nNonzero cells:\n")
        for group2_present in range(spectrum.shape[0]):
            for group1_present in range(spectrum.shape[1]):
                count = int(spectrum[group2_present, group1_present])
                if count:
                    handle.write(
                        f"  G1 present={group1_present}, "
                        f"G2 present={group2_present}: {count}\n"
                    )


def plot_heatmap(spectrum, group1_n, group2_n, output_file, fig_width=12.0,
                 fig_height=3.0, axis_fontsize=12.0, tick_fontsize=10.0,
                 annotation_fontsize=10.0, colorbar_fontsize=10.0,
                 fill_heatmap_height=False, match_colorbar_height=False):
    log_spectrum = np.log10(spectrum + 1)
    annot = spectrum.astype(str)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax = sns.heatmap(
        log_spectrum,
        ax=ax,
        annot=annot,
        fmt="",
        cmap="viridis",
        cbar=not match_colorbar_height,
        cbar_kws={"label": "log10(count + 1)"},
        linewidths=0.5,
        linecolor="white",
        xticklabels=list(range(group1_n + 1)),
        yticklabels=list(range(group2_n + 1)),
        square=not fill_heatmap_height,
        annot_kws={"fontsize": annotation_fontsize},
    )

    ax.set_xlabel("Population 1\nintroner count", fontsize=axis_fontsize)
    ax.set_ylabel("Population 2\nintroner count", fontsize=axis_fontsize)
    ax.set_title("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0,
                       fontsize=tick_fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,
                       fontsize=tick_fontsize)
    if match_colorbar_height:
        # Reserve space for a colorbar and then match its vertical bounds to
        # the square-cell heatmap after layout has finalized the axes box.
        fig.tight_layout(rect=(0, 0, 0.90, 1))
        fig.canvas.draw()
        heatmap_box = ax.get_position()
        colorbar_ax = fig.add_axes([
            heatmap_box.x1 + 0.025,
            heatmap_box.y0,
            0.020,
            heatmap_box.height,
        ])
        colorbar = fig.colorbar(ax.collections[0], cax=colorbar_ax)
    else:
        colorbar = ax.collections[0].colorbar
    colorbar.set_label("log10(count + 1)", fontsize=colorbar_fontsize)
    colorbar.ax.tick_params(labelsize=colorbar_fontsize)
    ax.invert_yaxis()
    if not match_colorbar_height:
        fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    df = pd.read_csv(args.genotype_matrix, sep="\t")

    spectrum, records, used, skipped = compute_present_afs(
        df,
        args.group1_samples,
        args.group2_samples,
        args.drop_absent_monomorphic,
    )

    for path in [args.output_pdf, args.output_png, args.output_tsv, args.summary]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    write_spectrum_tsv(spectrum, args.output_tsv)
    write_summary(args.summary, spectrum, used, skipped)
    plot_kwargs = {
        "fig_width": args.fig_width,
        "fig_height": args.fig_height,
        "axis_fontsize": args.axis_fontsize,
        "tick_fontsize": args.tick_fontsize,
        "annotation_fontsize": args.annotation_fontsize,
        "colorbar_fontsize": args.colorbar_fontsize,
        "fill_heatmap_height": args.fill_heatmap_height,
        "match_colorbar_height": args.match_colorbar_height,
    }
    plot_heatmap(spectrum, len(args.group1_samples), len(args.group2_samples),
                 args.output_pdf, **plot_kwargs)
    plot_heatmap(spectrum, len(args.group1_samples), len(args.group2_samples),
                 args.output_png, **plot_kwargs)

    print(f"Used {used} loci")
    for reason, count in skipped.items():
        print(f"Skipped {reason}: {count}")
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()

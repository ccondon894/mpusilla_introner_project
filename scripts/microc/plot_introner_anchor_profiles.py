#!/usr/bin/env python3
"""Aggregate Micro-C anchor signal around introner boundary classes."""

import argparse
import bisect
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PRESENT = 1
ABSENT = 2


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot aggregate Micro-C anchor signal around introner boundaries."
    )
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--anchor-bedgraph", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--normal-introns", required=True)
    parser.add_argument("--window-bp", type=int, required=True)
    parser.add_argument("--bin-bp", type=int, required=True)
    parser.add_argument("--absent-max-length", type=int, required=True)
    parser.add_argument("--normal-intron-exclusion-window-bp", type=int, required=True)
    parser.add_argument("--profile-tsv", required=True)
    parser.add_argument("--intron-comparison-profile-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--png", required=True)
    parser.add_argument("--pdf", required=True)
    parser.add_argument("--intron-comparison-png", required=True)
    parser.add_argument("--intron-comparison-pdf", required=True)
    return parser.parse_args()


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def classify_loci(matrix_path, sample, chrom_sizes, absent_max_length):
    df = pd.read_csv(matrix_path, sep="\t")
    df = df.loc[df["sample"].astype(str) == sample].copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["length"] = df["end"] - df["start"]
    df["presence"] = df["presence"].astype(int)

    loci = []
    summary = []
    for state, state_df in (
        ("present", df.loc[df["presence"] == PRESENT]),
        (
            "absent",
            df.loc[(df["presence"] == ABSENT) & (df["length"] <= absent_max_length)],
        ),
        (
            "absent_long_excluded",
            df.loc[(df["presence"] == ABSENT) & (df["length"] > absent_max_length)],
        ),
    ):
        state_df = state_df.loc[state_df["contig"].isin(chrom_sizes)].copy()
        summary.append(
            {
                "state": state,
                "n_loci": len(state_df),
                "min_length": state_df["length"].min() if len(state_df) else np.nan,
                "median_length": state_df["length"].median() if len(state_df) else np.nan,
                "max_length": state_df["length"].max() if len(state_df) else np.nan,
            }
        )
        for _, row in state_df.iterrows():
            loci.append(
                {
                    "state": state,
                    "chrom": row["contig"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "length": int(row["length"]),
                    "ortholog_id": str(row.get("ortholog_id", ".")),
                    "family": str(row.get("family", ".")),
                    "gene": str(row.get("gene", ".")) if pd.notna(row.get("gene", np.nan)) else ".",
                }
            )
    return loci, pd.DataFrame(summary)


def write_bed(loci, path):
    state_order = {"present": 0, "absent": 1, "absent_long_excluded": 2}
    rows = sorted(
        loci,
        key=lambda x: (x["chrom"], x["start"], x["end"], state_order.get(x["state"], 9)),
    )
    with open(path, "w") as handle:
        for locus in rows:
            name = (
                f"{locus['ortholog_id']}|{locus['state']}|"
                f"family_{locus['family']}|len_{locus['length']}"
            )
            handle.write(
                "\t".join(
                    [
                        locus["chrom"],
                        str(locus["start"]),
                        str(locus["end"]),
                        name,
                        "0",
                        ".",
                    ]
                )
                + "\n"
            )


def load_normal_introns(path, chrom_sizes):
    df = pd.read_csv(path, sep="\t")
    introns = []

    if {"contig", "start", "end", "intron_id"}.issubset(df.columns):
        if "presence" in df.columns:
            df = df.loc[df["presence"].astype(str).isin({"1", "1.0", "present"})]
        for _, row in df.iterrows():
            chrom = row["contig"]
            if chrom not in chrom_sizes:
                continue
            start = int(row["start"])
            end = int(row["end"])
            introns.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "intron_id": str(row.get("intron_id", ".")),
                    "gene": str(row.get("gene", ".")),
                }
            )
    elif {"gene_id", "intron_idx", "contig", "start", "end"}.issubset(df.columns):
        for _, row in df.iterrows():
            chrom = row["contig"]
            if chrom not in chrom_sizes:
                continue
            start = int(row["start"])
            end = int(row["end"])
            introns.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "intron_id": f"{row['gene_id']}__intron_{row['intron_idx']}",
                    "gene": str(row["gene_id"]),
                }
            )
    else:
        raise ValueError(
            f"Could not recognize normal intron file columns in {path}: "
            + ", ".join(df.columns)
        )
    return introns


def build_interval_index(loci, state):
    index = {}
    for locus in loci:
        if locus["state"] != state:
            continue
        index.setdefault(locus["chrom"], []).append((locus["start"], locus["end"]))
    for chrom, intervals in index.items():
        intervals.sort()
        index[chrom] = {
            "intervals": intervals,
            "starts": [start for start, _ in intervals],
        }
    return index


def overlaps_index(index, chrom, start, end):
    chrom_index = index.get(chrom)
    if not chrom_index:
        return False
    intervals = chrom_index["intervals"]
    starts = chrom_index["starts"]
    end_idx = bisect.bisect_right(starts, end)
    for intr_start, intr_end in intervals[:end_idx]:
        if intr_end >= start:
            return True
    return False


def filter_introns_away_from_introners(introns, introner_loci, exclusion_window):
    present_index = build_interval_index(introner_loci, "present")
    filtered = []
    excluded = 0
    for intron in introns:
        padded_start = intron["start"] - exclusion_window
        padded_end = intron["end"] + exclusion_window
        if overlaps_index(present_index, intron["chrom"], padded_start, padded_end):
            excluded += 1
            continue
        filtered.append(intron)
    return filtered, excluded


def build_profile_windows(points, chrom_sizes, window_bp, states):
    windows_by_chrom = {}
    state_counts = {state: 0 for state in states}
    for point in points:
        if point["state"] not in state_counts:
            continue
        chrom = point["chrom"]
        chrom_size = chrom_sizes[chrom]
        window_start = point["center"] - window_bp
        window_end = point["center"] + window_bp
        if window_end <= 0 or window_start >= chrom_size:
            continue
        state = point["state"]
        state_idx = state_counts[state]
        state_counts[state] += 1
        windows_by_chrom.setdefault(chrom, []).append(
            {
                "state": state,
                "state_idx": state_idx,
                "start": window_start,
                "end": window_end,
            }
        )

    for chrom, windows in windows_by_chrom.items():
        windows.sort(key=lambda x: x["start"])
        windows_by_chrom[chrom] = {
            "windows": windows,
            "starts": [window["start"] for window in windows],
        }
    return windows_by_chrom, state_counts


def add_interval_to_window_sums(
    state_matrices, interval_start, interval_end, value, window, window_bp, bin_bp
):
    overlap_start = max(interval_start, window["start"])
    overlap_end = min(interval_end, window["end"])
    if overlap_start >= overlap_end:
        return

    rel_start = overlap_start - window["start"]
    rel_end = overlap_end - window["start"]
    first_bin = max(0, rel_start // bin_bp)
    last_bin = min((2 * window_bp) // bin_bp - 1, (rel_end - 1) // bin_bp)

    for bin_idx in range(first_bin, last_bin + 1):
        bin_start = bin_idx * bin_bp
        bin_end = bin_start + bin_bp
        bin_overlap = max(0, min(rel_end, bin_end) - max(rel_start, bin_start))
        if bin_overlap:
            state_matrices[window["state"]][window["state_idx"], bin_idx] += (
                value * bin_overlap
            )


def iter_bedgraph(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            yield chrom, int(start), int(end), float(value)


def aggregate_profiles(points, states, anchor_bedgraph, chrom_sizes, window_bp, bin_bp):
    n_bins = (2 * window_bp) // bin_bp
    rel_starts = np.arange(-window_bp, window_bp, bin_bp)
    windows_by_chrom, state_counts = build_profile_windows(
        points, chrom_sizes, window_bp, states
    )
    state_matrices = {
        state: np.zeros((state_counts[state], n_bins), dtype=float)
        for state in states
    }

    state_n = {state: state_counts[state] for state in states}
    state_n_raw = {
        state: sum(1 for point in points if point["state"] == state)
        for state in states
    }
    for state in states:
        if state_n[state] != state_n_raw[state]:
            print(
                f"Warning: skipped {state_n_raw[state] - state_n[state]} "
                f"{state} points whose windows fall outside chromosome bounds"
            )

    cursor_by_chrom = {chrom: 0 for chrom in windows_by_chrom}

    for chrom, interval_start, interval_end, value in iter_bedgraph(anchor_bedgraph):
        chrom_windows = windows_by_chrom.get(chrom)
        if not chrom_windows:
            continue
        windows = chrom_windows["windows"]
        starts = chrom_windows["starts"]
        cursor = cursor_by_chrom[chrom]
        while cursor < len(windows) and windows[cursor]["end"] <= interval_start:
            cursor += 1
        cursor_by_chrom[chrom] = cursor
        end_idx = bisect.bisect_left(starts, interval_end)
        for idx in range(cursor, end_idx):
            window = windows[idx]
            add_interval_to_window_sums(
                state_matrices,
                interval_start,
                interval_end,
                value,
                window,
                window_bp,
                bin_bp,
            )

    rows = []
    for state in states:
        if state_n[state] == 0:
            continue
        matrix = state_matrices[state] / bin_bp
        means = matrix.mean(axis=0)
        if matrix.shape[0] == 1:
            sems = np.zeros_like(means)
        else:
            sems = matrix.std(axis=0, ddof=1) / np.sqrt(matrix.shape[0])
        for rel_start, mean, sem in zip(rel_starts, means, sems):
            rows.append(
                {
                    "state": state,
                    "rel_start_bp": int(rel_start),
                    "rel_end_bp": int(rel_start + bin_bp),
                    "rel_mid_bp": int(rel_start + bin_bp / 2),
                    "n_loci": state_n[state],
                    "mean_anchor_cpm": mean,
                    "sem_anchor_cpm": sem,
                    "ci95_anchor_cpm": 1.96 * sem,
                }
            )
    return pd.DataFrame(rows)


def profile_points(loci, normal_introns):
    points = []
    for locus in loci:
        if locus["state"] in {"present", "absent"}:
            state_prefix = locus["state"]
            points.extend(
                [
                    {
                        "state": f"{state_prefix}_start_plus_100",
                        "chrom": locus["chrom"],
                        "center": locus["start"] + 100,
                    },
                    {
                        "state": f"{state_prefix}_end_minus_100",
                        "chrom": locus["chrom"],
                        "center": locus["end"] - 100,
                    },
                ]
            )
    points.extend(
        point
        for intron in normal_introns
        for point in (
            {
                "state": "normal_intron_start",
                "chrom": intron["chrom"],
                "center": intron["start"],
            },
            {
                "state": "normal_intron_end",
                "chrom": intron["chrom"],
                "center": intron["end"],
            },
        )
    )
    return points


def stitch_composite_boundary_profile(profile, specs):
    """Build composite traces from start anchors on the left and end anchors on the right."""
    stitched = []
    for output_state, start_state, end_state in specs:
        start_half = profile.loc[
            (profile["state"] == start_state) & (profile["rel_mid_bp"] < 0)
        ].copy()
        end_half = profile.loc[
            (profile["state"] == end_state) & (profile["rel_mid_bp"] >= 0)
        ].copy()
        for source, subset in (
            ("start", start_half),
            ("end", end_half),
        ):
            subset["state"] = output_state
            subset["source_anchor"] = source
            stitched.append(subset)
    return pd.concat(stitched, ignore_index=True).sort_values(["state", "rel_mid_bp"])


def plot_profile(
    profile,
    summary,
    sample,
    window_bp,
    png,
    pdf,
    states,
    colors,
    labels,
    title,
    xlabel,
    note=None,
):
    fig, ax = plt.subplots(figsize=(8.0, 4.6))

    for state in states:
        subset = profile.loc[profile["state"] == state].copy()
        if subset.empty:
            continue
        x = subset["rel_mid_bp"].to_numpy()
        y = subset["mean_anchor_cpm"].to_numpy()
        ci = subset["ci95_anchor_cpm"].to_numpy()
        n_loci = int(subset["n_loci"].iloc[0])
        ax.plot(x, y, color=colors[state], lw=2.0, label=f"{labels[state]} (n={n_loci:,})")
        if not np.isnan(ci).all():
            ax.fill_between(x, y - ci, y + ci, color=colors[state], alpha=0.18, linewidth=0)

    ax.axvline(0, color="#222222", lw=1.0, ls="--", alpha=0.65)
    ax.set_xlim(-window_bp, window_bp)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mean anchor signal (CPM)")
    ax.set_title(f"{sample} {title}")
    if note:
        ax.text(
            0.01,
            0.98,
            note,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            color="#555555",
        )
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


def main():
    args = parse_args()
    if (2 * args.window_bp) % args.bin_bp != 0:
        raise ValueError("--bin-bp must divide 2 * --window-bp")

    for path in (
        args.profile_tsv,
        args.intron_comparison_profile_tsv,
        args.summary_tsv,
        args.bed,
        args.png,
        args.pdf,
        args.intron_comparison_png,
        args.intron_comparison_pdf,
    ):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    loci, summary = classify_loci(
        args.genotype_matrix, args.sample, chrom_sizes, args.absent_max_length
    )
    write_bed(loci, args.bed)
    normal_introns_raw = load_normal_introns(args.normal_introns, chrom_sizes)
    normal_introns, normal_introns_excluded = filter_introns_away_from_introners(
        normal_introns_raw, loci, args.normal_intron_exclusion_window_bp
    )

    aggregate_profile = aggregate_profiles(
        profile_points(loci, normal_introns),
        [
            "present_start_plus_100",
            "present_end_minus_100",
            "absent_start_plus_100",
            "absent_end_minus_100",
            "normal_intron_start",
            "normal_intron_end",
        ],
        args.anchor_bedgraph,
        chrom_sizes,
        args.window_bp,
        args.bin_bp,
    )
    boundary_profile = stitch_composite_boundary_profile(
        aggregate_profile,
        [
            (
                "present_boundary_internal",
                "present_start_plus_100",
                "present_end_minus_100",
            ),
            ("absent_boundary_internal", "absent_start_plus_100", "absent_end_minus_100"),
        ],
    )
    intron_comparison_profile = stitch_composite_boundary_profile(
        aggregate_profile,
        [
            (
                "present_boundary_internal",
                "present_start_plus_100",
                "present_end_minus_100",
            ),
            ("normal_intron_boundary", "normal_intron_start", "normal_intron_end"),
        ],
    )
    boundary_profile.to_csv(args.profile_tsv, sep="\t", index=False)
    intron_comparison_profile.to_csv(
        args.intron_comparison_profile_tsv, sep="\t", index=False
    )
    summary.to_csv(args.summary_tsv, sep="\t", index=False)

    excluded = summary.loc[summary["state"] == "absent_long_excluded", "n_loci"]
    excluded_n = int(excluded.iloc[0]) if len(excluded) else 0
    plot_profile(
        boundary_profile,
        summary,
        args.sample,
        args.window_bp,
        args.png,
        args.pdf,
        ["present_boundary_internal", "absent_boundary_internal"],
        {
            "present_boundary_internal": "#2d6cdf",
            "absent_boundary_internal": "#d65f2d",
        },
        {
            "present_boundary_internal": "Present introners (left=start+100; right=end-100)",
            "absent_boundary_internal": "Absent loci (left=start+100; right=end-100)",
        },
        "Micro-C anchor signal around introner internal boundary positions",
        "Composite coordinate: left=start+100 bp; right=end-100 bp",
        note=f"Excluded absent loci >200 bp: {excluded_n:,}",
    )
    plot_profile(
        intron_comparison_profile,
        summary,
        args.sample,
        args.window_bp,
        args.intron_comparison_png,
        args.intron_comparison_pdf,
        ["present_boundary_internal", "normal_intron_boundary"],
        {"present_boundary_internal": "#2d6cdf", "normal_intron_boundary": "#6f6f6f"},
        {
            "present_boundary_internal": "Present introners (left=start+100; right=end-100)",
            "normal_intron_boundary": "Normal introns (left=start; right=end)",
        },
        "Micro-C anchor signal around introner and normal intron boundaries",
        "Composite coordinate: left=start-like boundary; right=end-like boundary",
        note=(
            "Normal introns excluded within "
            f"{args.normal_intron_exclusion_window_bp:,} bp of present introners: "
            f"{normal_introns_excluded:,}"
        ),
    )


if __name__ == "__main__":
    main()

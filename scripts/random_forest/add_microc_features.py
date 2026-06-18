#!/usr/bin/env python3
"""Add CCMP1545 Micro-C insulation/domain/boundary features to a feature matrix."""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        default="genotype_matrix_with_features.tsv",
        help="Input feature matrix TSV.",
    )
    parser.add_argument(
        "--microc-dir",
        default="/scratch1/chris/mpusilla_introner_project/results/microc",
        help="Micro-C results directory.",
    )
    parser.add_argument("--sample", default="CCMP1545")
    parser.add_argument("--resolution", type=int, default=1000)
    parser.add_argument("--window", type=int, default=25000)
    parser.add_argument("--boundary-flank-bp", type=int, default=5000)
    parser.add_argument(
        "--flank-windows",
        default="25,50,250,500",
        help="Comma-separated flank windows for local Micro-C track summaries.",
    )
    parser.add_argument(
        "--skip-flank-tracks",
        action="store_true",
        help="Only add insulation/domain/boundary features, not local track features.",
    )
    parser.add_argument(
        "--mummer-coords",
        default=None,
        help=(
            "Optional MUMmer show-coords file. Query contigs with mostly reverse "
            "alignments will have left/right Micro-C flank tracks swapped."
        ),
    )
    parser.add_argument(
        "--orientation-threshold",
        type=float,
        default=0.8,
        help="Minimum reverse-aligned fraction for treating a query contig as reversed.",
    )
    parser.add_argument(
        "--output",
        default="genotype_matrix_with_features.microc.tsv",
        help="Output feature matrix TSV.",
    )
    parser.add_argument(
        "--summary",
        default="microc_feature_join_summary.tsv",
        help="Output join summary TSV.",
    )
    return parser.parse_args()


def parse_windows(value):
    return [int(window.strip()) for window in value.split(",") if window.strip()]


def read_reverse_oriented_contigs(coords_path, threshold):
    if not coords_path:
        return set()

    by_contig = {}
    with open(coords_path) as handle:
        for line in handle:
            line = line.strip()
            if (
                not line
                or line.startswith("/")
                or line.startswith("NUCMER")
                or line.startswith("[")
            ):
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            try:
                s2, e2 = int(parts[2]), int(parts[3])
                len1, len2 = int(parts[4]), int(parts[5])
            except ValueError:
                continue
            query_contig = parts[12]
            aligned_bp = min(len1, len2)
            if aligned_bp <= 0:
                continue
            orientation = by_contig.setdefault(query_contig, {"forward": 0, "reverse": 0})
            if s2 > e2:
                orientation["reverse"] += aligned_bp
            else:
                orientation["forward"] += aligned_bp

    reverse_contigs = set()
    for contig, counts in by_contig.items():
        total = counts["forward"] + counts["reverse"]
        if total and counts["reverse"] / total >= threshold:
            reverse_contigs.add(contig)
    return reverse_contigs


def feature_start_end(df):
    start_col = "feature_start" if "feature_start" in df.columns else "start"
    end_col = "feature_end" if "feature_end" in df.columns else "end"
    return df[start_col].astype(int), df[end_col].astype(int)


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def score_column(columns, prefix, window):
    exact = f"{prefix}_{window}"
    if exact in columns:
        return exact
    matches = [column for column in columns if column.startswith(prefix)]
    if not matches:
        raise ValueError(f"Could not find a column starting with {prefix!r}")
    return matches[0]


def load_interval_table(path, names):
    return pd.read_csv(path, sep="\t", header=None, names=names)


def load_dense_bedgraph(path, chrom_sizes):
    arrays = {
        chrom: np.zeros(size, dtype=np.float32)
        for chrom, size in chrom_sizes.items()
    }

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            if chrom not in arrays:
                continue
            start = max(0, int(start))
            end = min(int(end), len(arrays[chrom]))
            if end <= start:
                continue
            arrays[chrom][start:end] = float(value)
    return arrays


def summarize_dense_track(chroms, starts, ends, track_arrays):
    means = np.full(len(starts), np.nan)
    covered_fractions = np.full(len(starts), np.nan)

    for idx, (chrom, start, end) in enumerate(zip(chroms, starts, ends)):
        arr = track_arrays.get(chrom)
        if arr is None:
            continue
        start = max(0, int(start))
        end = min(int(end), len(arr))
        if end <= start:
            continue
        values = arr[start:end]
        means[idx] = float(values.mean())
        covered_fractions[idx] = float(np.count_nonzero(values) / len(values))
    return means, covered_fractions


def add_flank_track_features(df, microc_dir, sample, flank_windows, reverse_contigs=None):
    reverse_contigs = reverse_contigs or set()
    chrom_sizes_path = Path(microc_dir) / "references" / f"{sample}.chrom.sizes"
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)

    track_paths = {
        "anchors": Path(microc_dir) / "tracks" / f"{sample}.anchors.cpm.bedgraph.gz",
        "coverage": Path(microc_dir) / "tracks" / f"{sample}.coverage.cpm.bedgraph.gz",
    }

    chroms = df["contig"].astype(str).to_numpy()
    reverse_flags = df["contig"].isin(reverse_contigs).to_numpy()
    starts, ends = feature_start_end(df)
    locus_starts = starts.to_numpy()
    locus_ends = ends.to_numpy()
    feature_blocks = []

    for track_name, track_path in track_paths.items():
        print(f"Loading {track_name} track from {track_path}")
        track_arrays = load_dense_bedgraph(track_path, chrom_sizes)
        track_features = {}
        prefix = f"microc_{sample}_{track_name}_cpm"

        for window in flank_windows:
            lower_left_starts = locus_starts - window
            lower_left_ends = locus_starts
            upper_right_starts = locus_ends
            upper_right_ends = locus_ends + window

            left_starts = np.where(reverse_flags, upper_right_starts, lower_left_starts)
            left_ends = np.where(reverse_flags, upper_right_ends, lower_left_ends)
            right_starts = np.where(reverse_flags, lower_left_starts, upper_right_starts)
            right_ends = np.where(reverse_flags, lower_left_ends, upper_right_ends)

            left_mean, left_covered = summarize_dense_track(
                chroms,
                left_starts,
                left_ends,
                track_arrays,
            )
            right_mean, right_covered = summarize_dense_track(
                chroms,
                right_starts,
                right_ends,
                track_arrays,
            )
            combined_mean = (left_mean + right_mean) / 2
            combined_covered = (left_covered + right_covered) / 2

            track_features[f"{prefix}_mean_left_{window}"] = left_mean
            track_features[f"{prefix}_mean_right_{window}"] = right_mean
            track_features[f"{prefix}_mean_combined_{window}"] = combined_mean
            track_features[f"{prefix}_mean_delta_{window}"] = left_mean - right_mean
            track_features[f"{prefix}_covered_fraction_left_{window}"] = left_covered
            track_features[f"{prefix}_covered_fraction_right_{window}"] = right_covered
            track_features[f"{prefix}_covered_fraction_combined_{window}"] = combined_covered
            track_features[f"{prefix}_covered_fraction_delta_{window}"] = (
                left_covered - right_covered
            )

        feature_blocks.append(pd.DataFrame(track_features, index=df.index))

    return pd.concat([df] + feature_blocks, axis=1)


def build_insulation_lookup(path, window):
    df = pd.read_csv(path, sep="\t")
    score_col = score_column(df.columns, "log2_insulation_score", window)
    valid_col = score_column(df.columns, "n_valid_pixels", window)
    strength_col = score_column(df.columns, "boundary_strength", window)
    boundary_col = score_column(df.columns, "is_boundary", window)

    lookup = {}
    for chrom, data in df.groupby("chrom", sort=False):
        data = data.sort_values("start").reset_index(drop=True)
        lookup[chrom] = {
            "starts": data["start"].to_numpy(dtype=int),
            "ends": data["end"].to_numpy(dtype=int),
            "scores": data[score_col].to_numpy(dtype=float),
            "n_valid_pixels": data[valid_col].to_numpy(dtype=float),
            "boundary_strengths": data[strength_col].to_numpy(dtype=float),
            "is_boundary": data[boundary_col].fillna(False).astype(bool).to_numpy(),
        }
    return lookup


def build_boundary_lookup(path):
    names = ["chrom", "start", "end", "name", "score", "strand", "strength"]
    df = load_interval_table(path, names)
    lookup = {}
    for chrom, data in df.groupby("chrom", sort=False):
        data = data.sort_values("start").reset_index(drop=True)
        lookup[chrom] = {
            "starts": data["start"].to_numpy(dtype=int),
            "ends": data["end"].to_numpy(dtype=int),
            "strengths": data["strength"].to_numpy(dtype=float),
        }
    return lookup


def build_domain_lookup(path):
    names = ["chrom", "start", "end", "name", "score", "strand"]
    df = load_interval_table(path, names)
    lookup = {}
    for chrom, data in df.groupby("chrom", sort=False):
        data = data.sort_values("start").reset_index(drop=True)
        lookup[chrom] = {
            "starts": data["start"].to_numpy(dtype=int),
            "ends": data["end"].to_numpy(dtype=int),
            "names": data["name"].astype(str).to_numpy(),
        }
    return lookup


def lookup_insulation(chroms, midpoints, lookup):
    n = len(midpoints)
    scores = np.full(n, np.nan)
    valid_pixels = np.full(n, np.nan)
    boundary_strengths = np.full(n, np.nan)
    is_boundary = np.zeros(n, dtype=bool)

    for idx, (chrom, midpoint) in enumerate(zip(chroms, midpoints)):
        chrom_lookup = lookup.get(chrom)
        if chrom_lookup is None:
            continue
        bin_idx = np.searchsorted(chrom_lookup["starts"], int(midpoint), side="right") - 1
        if bin_idx < 0 or bin_idx >= len(chrom_lookup["starts"]):
            continue
        if int(midpoint) >= chrom_lookup["ends"][bin_idx]:
            continue
        scores[idx] = chrom_lookup["scores"][bin_idx]
        valid_pixels[idx] = chrom_lookup["n_valid_pixels"][bin_idx]
        boundary_strengths[idx] = chrom_lookup["boundary_strengths"][bin_idx]
        is_boundary[idx] = bool(chrom_lookup["is_boundary"][bin_idx])

    return scores, valid_pixels, boundary_strengths, is_boundary.astype(int)


def lookup_nearest_boundary(chroms, midpoints, lookup, flank_bp):
    n = len(midpoints)
    distances = np.full(n, np.nan)
    strengths = np.full(n, np.nan)
    within_flank = np.zeros(n, dtype=int)

    for idx, (chrom, midpoint) in enumerate(zip(chroms, midpoints)):
        chrom_lookup = lookup.get(chrom)
        if chrom_lookup is None or len(chrom_lookup["starts"]) == 0:
            continue

        midpoint = int(midpoint)
        starts = chrom_lookup["starts"]
        ends = chrom_lookup["ends"]
        candidate_idxs = set()
        insert_idx = np.searchsorted(starts, midpoint, side="right")
        for candidate_idx in (insert_idx - 2, insert_idx - 1, insert_idx, insert_idx + 1):
            if 0 <= candidate_idx < len(starts):
                candidate_idxs.add(candidate_idx)

        best = None
        for candidate_idx in candidate_idxs:
            start = starts[candidate_idx]
            end = ends[candidate_idx]
            if start <= midpoint < end:
                distance = 0
            elif midpoint < start:
                distance = start - midpoint
            else:
                distance = midpoint - end
            if best is None or distance < best[0]:
                best = (distance, chrom_lookup["strengths"][candidate_idx])

        if best is None:
            continue
        distances[idx] = best[0]
        strengths[idx] = best[1]
        within_flank[idx] = int(best[0] <= flank_bp)

    return distances, strengths, within_flank


def lookup_domain(chroms, midpoints, lookup):
    n = len(midpoints)
    domain_lengths = np.full(n, np.nan)
    relative_positions = np.full(n, np.nan)
    terminality = np.full(n, np.nan)
    dist_left = np.full(n, np.nan)
    dist_right = np.full(n, np.nan)
    dist_nearest_edge = np.full(n, np.nan)
    names = np.full(n, "", dtype=object)

    for idx, (chrom, midpoint) in enumerate(zip(chroms, midpoints)):
        chrom_lookup = lookup.get(chrom)
        if chrom_lookup is None:
            continue
        midpoint = int(midpoint)
        domain_idx = np.searchsorted(chrom_lookup["starts"], midpoint, side="right") - 1
        if domain_idx < 0 or domain_idx >= len(chrom_lookup["starts"]):
            continue
        start = chrom_lookup["starts"][domain_idx]
        end = chrom_lookup["ends"][domain_idx]
        if midpoint >= end:
            continue
        length = end - start
        if length <= 0:
            continue
        left = midpoint - start
        right = end - midpoint
        rel_pos = left / length

        domain_lengths[idx] = length
        relative_positions[idx] = rel_pos
        terminality[idx] = abs(rel_pos - 0.5) * 2
        dist_left[idx] = left
        dist_right[idx] = right
        dist_nearest_edge[idx] = min(left, right)
        names[idx] = chrom_lookup["names"][domain_idx]

    return (
        domain_lengths,
        relative_positions,
        terminality,
        dist_left,
        dist_right,
        dist_nearest_edge,
        names,
    )


def add_microc_features(
    matrix,
    microc_dir,
    sample,
    resolution,
    window,
    boundary_flank_bp,
    flank_windows,
    add_flank_tracks=True,
    reverse_contigs=None,
):
    reverse_contigs = reverse_contigs or set()
    prefix = f"{sample}.{resolution}bp_{window}bp"
    insulation_dir = Path(microc_dir) / "qc" / "insulation"
    insulation_path = insulation_dir / f"{prefix}.insulation.tsv"
    boundaries_path = insulation_dir / f"{prefix}.boundaries.bed"
    domains_path = insulation_dir / f"{prefix}.domains.bed"

    df = pd.read_csv(matrix, sep="\t")
    feature_starts, feature_ends = feature_start_end(df)
    feature_block = pd.DataFrame(index=df.index)
    feature_block["microc_midpoint"] = (
        (feature_starts + feature_ends) / 2
    ).astype(int)

    insulation_lookup = build_insulation_lookup(insulation_path, window)
    boundary_lookup = build_boundary_lookup(boundaries_path)
    domain_lookup = build_domain_lookup(domains_path)

    chroms = df["contig"].astype(str).to_numpy()
    reverse_flags = df["contig"].isin(reverse_contigs).to_numpy()
    midpoints = feature_block["microc_midpoint"].to_numpy(dtype=int)
    col_prefix = f"microc_{sample}_{resolution}bp_{window}bp"

    (
        feature_block[f"{col_prefix}_insulation_score"],
        feature_block[f"{col_prefix}_n_valid_pixels"],
        feature_block[f"{col_prefix}_bin_boundary_strength"],
        feature_block[f"{col_prefix}_bin_is_boundary"],
    ) = lookup_insulation(chroms, midpoints, insulation_lookup)

    (
        feature_block[f"{col_prefix}_nearest_boundary_distance_bp"],
        feature_block[f"{col_prefix}_nearest_boundary_strength"],
        feature_block[f"{col_prefix}_within_{boundary_flank_bp}bp_boundary"],
    ) = lookup_nearest_boundary(chroms, midpoints, boundary_lookup, boundary_flank_bp)

    (
        domain_lengths,
        domain_relative_positions,
        domain_terminality,
        domain_dist_left,
        domain_dist_right,
        domain_dist_nearest_edge,
        domain_names,
    ) = lookup_domain(chroms, midpoints, domain_lookup)
    domain_relative_positions = np.where(
        reverse_flags & ~np.isnan(domain_relative_positions),
        1 - domain_relative_positions,
        domain_relative_positions,
    )
    left_copy = domain_dist_left.copy()
    domain_dist_left = np.where(reverse_flags, domain_dist_right, domain_dist_left)
    domain_dist_right = np.where(reverse_flags, left_copy, domain_dist_right)

    feature_block[f"{col_prefix}_domain_length_bp"] = domain_lengths
    feature_block[f"{col_prefix}_domain_relative_position"] = domain_relative_positions
    feature_block[f"{col_prefix}_domain_terminality"] = domain_terminality
    feature_block[f"{col_prefix}_distance_to_domain_left_bp"] = domain_dist_left
    feature_block[f"{col_prefix}_distance_to_domain_right_bp"] = domain_dist_right
    feature_block[f"{col_prefix}_distance_to_nearest_domain_edge_bp"] = domain_dist_nearest_edge
    feature_block[f"{col_prefix}_domain_name"] = domain_names

    df = pd.concat([df, feature_block], axis=1)
    if add_flank_tracks:
        df = add_flank_track_features(df, microc_dir, sample, flank_windows, reverse_contigs)

    summary = []
    microc_cols = [
        col for col in df.columns
        if col.startswith("microc_")
        and col != "microc_midpoint"
        and pd.api.types.is_numeric_dtype(df[col])
    ]
    for col in microc_cols:
        if col.endswith("_domain_name"):
            continue
        values = df[col]
        summary.append(
            {
                "feature": col,
                "n_rows": len(values),
                "n_missing": int(values.isna().sum()),
                "missing_fraction": float(values.isna().mean()),
                "mean_label_0": values[df["label"] == 0].mean(),
                "mean_label_1": values[df["label"] == 1].mean(),
            }
        )

    return df, pd.DataFrame(summary)


def main():
    args = parse_args()
    reverse_contigs = read_reverse_oriented_contigs(
        args.mummer_coords,
        args.orientation_threshold,
    )
    if reverse_contigs:
        print(f"Orientation-normalizing {len(reverse_contigs)} reverse-oriented contigs.")
    df, summary = add_microc_features(
        args.matrix,
        args.microc_dir,
        args.sample,
        args.resolution,
        args.window,
        args.boundary_flank_bp,
        parse_windows(args.flank_windows),
        add_flank_tracks=not args.skip_flank_tracks,
        reverse_contigs=reverse_contigs,
    )
    df.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)
    print(f"Wrote augmented matrix to {args.output}")
    print(f"Wrote join summary to {args.summary}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

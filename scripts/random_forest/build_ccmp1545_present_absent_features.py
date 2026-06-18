#!/usr/bin/env python3
"""Build per-sample present-vs-absent introner locus features."""

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from build_preinsertion_features import (
    KMERS_3,
    add_window_features,
    fetch_interval,
    kmer_frequencies,
    mask_junction_bases,
    parse_windows,
    scalar_features,
)
from gtf_gene_annotation import GtfGeneAnnotator


MATING_TYPE_CONTIG = "CCMP1545#0#scaffold_2"
REVCOMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BOUNDARY_SUPPORT_COLS = [
    "boundary_support_ortholog_id",
    "boundary_support_best_junction_name",
    "boundary_support_best_junction_score",
    "boundary_support_start_delta",
    "boundary_support_end_delta",
    "boundary_support_max_abs_boundary_delta",
    "boundary_support_sum_abs_boundary_delta",
    "boundary_support_applied",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--present-feature-matrix",
        default=None,
        help="Boundary-corrected CCMP1545 feature matrix containing present introners.",
    )
    parser.add_argument(
        "--genotype-matrix",
        default="../../results/genotyping/genotype_matrix.final.tsv",
    )
    parser.add_argument(
        "--fa",
        default="../../results/assemblies/CCMP1545.vg_paths.fa",
    )
    parser.add_argument(
        "--gtf",
        default="../../results/annotations/CCMP1545.gtf",
    )
    parser.add_argument("--sample", default="CCMP1545")
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument(
        "--introner-boundary-support",
        default=None,
        help="Optional regtools per-locus boundary support table for present introner loci.",
    )
    parser.add_argument(
        "--boundary-max-abs-delta",
        type=int,
        default=25,
        help="Apply supported introner junctions only when max_abs_boundary_delta is <= this value.",
    )
    parser.add_argument(
        "--junction-mask-bp",
        type=int,
        default=3,
        help="Mask this many bases adjacent to present boundaries or absent midpoint.",
    )
    parser.add_argument(
        "--exclude-mating-type-contig",
        action="store_true",
        default=True,
        help="Exclude the CCMP1545 mating-type contig from present and absent loci.",
    )
    parser.add_argument(
        "--exclude-contig",
        action="append",
        default=[],
        help="Additional contig to exclude; can be provided more than once.",
    )
    parser.add_argument(
        "--mummer-coords",
        default=None,
        help=(
            "Optional MUMmer show-coords file. Query contigs with mostly reverse "
            "alignments will have left/right flanks reverse-complemented and swapped."
        ),
    )
    parser.add_argument(
        "--orientation-threshold",
        type=float,
        default=0.8,
        help="Minimum reverse-aligned fraction for treating a query contig as reversed.",
    )
    parser.add_argument(
        "--max-absent-span",
        type=int,
        default=None,
        help="Optional maximum absent interval span to retain for class 0 loci.",
    )
    parser.add_argument(
        "--output",
        default="rf_results/ccmp1545_present_absent_features.mask3.w25_50_100.tsv",
    )
    parser.add_argument(
        "--summary",
        default="rf_results/ccmp1545_present_absent_features.mask3.w25_50_100.summary.tsv",
    )
    return parser.parse_args()


def sample_status(genotype_df, sample):
    keep_cols = [
        "ortholog_id",
        "sequence_id",
        "contig",
        "start",
        "end",
        "family",
        "gene",
        "presence",
        "group1_pattern",
        "group1_present_count",
        "group1_absent_count",
        "group1_missing_count",
        "group1_callable_count",
        "group1_n_samples",
        "within_group_status",
        "within_group_identity",
    ]
    status = genotype_df.loc[
        genotype_df["sample"].eq(sample),
        [col for col in keep_cols if col in genotype_df.columns],
    ].copy()
    status["body_start"] = status["start"].astype(int) + 100
    status["body_end"] = status["end"].astype(int) - 100
    return status


def read_boundary_support(path, max_abs_delta):
    if not path:
        return None

    boundary_df = pd.read_csv(path, sep="\t")
    key_cols = ["contig", "body_start", "body_end"]
    sort_cols = key_cols + ["best_junction_score", "max_abs_boundary_delta", "sum_abs_boundary_delta"]
    boundary_df = (
        boundary_df.sort_values(
            sort_cols,
            ascending=[True, True, True, False, True, True],
        )
        .drop_duplicates(key_cols, keep="first")
        .copy()
    )

    boundary_df["boundary_support_applied"] = (
        boundary_df["max_abs_boundary_delta"].le(max_abs_delta).astype(int)
    )
    boundary_df["boundary_feature_start"] = boundary_df["body_start"].where(
        boundary_df["boundary_support_applied"].eq(0),
        boundary_df["best_junction_start"],
    )
    boundary_df["boundary_feature_end"] = boundary_df["body_end"].where(
        boundary_df["boundary_support_applied"].eq(0),
        boundary_df["best_junction_end"],
    )

    keep = key_cols + [
        "ortholog_id",
        "best_junction_name",
        "best_junction_score",
        "start_delta",
        "end_delta",
        "max_abs_boundary_delta",
        "sum_abs_boundary_delta",
        "boundary_support_applied",
        "boundary_feature_start",
        "boundary_feature_end",
    ]
    boundary_df = boundary_df[keep].rename(columns={
        "ortholog_id": "boundary_support_ortholog_id",
        "best_junction_name": "boundary_support_best_junction_name",
        "best_junction_score": "boundary_support_best_junction_score",
        "start_delta": "boundary_support_start_delta",
        "end_delta": "boundary_support_end_delta",
        "max_abs_boundary_delta": "boundary_support_max_abs_boundary_delta",
        "sum_abs_boundary_delta": "boundary_support_sum_abs_boundary_delta",
    })
    return boundary_df


def apply_boundary_support(present, boundary_support):
    present["feature_start"] = present["start"]
    present["feature_end"] = present["end"]
    if boundary_support is None:
        present["boundary_support_applied"] = 0
        for col in BOUNDARY_SUPPORT_COLS:
            if col not in present:
                present[col] = pd.NA
        return present

    present = present.drop(columns=BOUNDARY_SUPPORT_COLS, errors="ignore")
    out = present.merge(
        boundary_support,
        left_on=["contig", "start", "end"],
        right_on=["contig", "body_start", "body_end"],
        how="left",
    )
    out["feature_start"] = out["boundary_feature_start"].fillna(out["start"]).astype(int)
    out["feature_end"] = out["boundary_feature_end"].fillna(out["end"]).astype(int)
    out["boundary_support_applied"] = out["boundary_support_applied"].fillna(0).astype(int)
    for col in BOUNDARY_SUPPORT_COLS:
        if col not in out:
            out[col] = pd.NA
    return out.drop(
        columns=[
            col for col in ["body_start", "body_end", "boundary_feature_start", "boundary_feature_end"]
            if col in out
        ],
        errors="ignore",
    )


def reverse_complement(seq):
    return seq.translate(REVCOMP_TABLE)[::-1].upper()


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


def fetch_oriented_boundary_flanks(
    fasta_index,
    contig,
    start,
    end,
    window,
    reverse_oriented=False,
):
    if reverse_oriented:
        left = reverse_complement(fetch_interval(fasta_index, contig, end, end + window))
        right = reverse_complement(fetch_interval(fasta_index, contig, start - window, start))
    else:
        left = fetch_interval(fasta_index, contig, start - window, start)
        right = fetch_interval(fasta_index, contig, end, end + window)
    return left, right


def fetch_oriented_site_flanks(
    fasta_index,
    contig,
    site,
    window,
    reverse_oriented=False,
):
    if reverse_oriented:
        left = reverse_complement(fetch_interval(fasta_index, contig, site, site + window))
        right = reverse_complement(fetch_interval(fasta_index, contig, site - window, site))
    else:
        left = fetch_interval(fasta_index, contig, site - window, site)
        right = fetch_interval(fasta_index, contig, site, site + window)
    return left, right


def features_from_flanks(left, right, window):
    features = {}
    combined = left + right

    left_scalars = scalar_features(left)
    right_scalars = scalar_features(right)
    combined_scalars = scalar_features(combined)
    for name, value in left_scalars.items():
        features[f"{name}_left_{window}"] = value
    for name, value in right_scalars.items():
        features[f"{name}_right_{window}"] = value
    for name, value in combined_scalars.items():
        features[f"{name}_combined_{window}"] = value
    for name in left_scalars:
        features[f"{name}_delta_{window}"] = left_scalars[name] - right_scalars[name]

    left_kmers = kmer_frequencies(left)
    right_kmers = kmer_frequencies(right)
    for kmer in KMERS_3:
        features[f"kmer_{kmer}_left_{window}"] = left_kmers[kmer]
        features[f"kmer_{kmer}_right_{window}"] = right_kmers[kmer]
        features[f"kmer_{kmer}_delta_{window}"] = left_kmers[kmer] - right_kmers[kmer]
    return features


def add_present_window_features(row, fasta_index, windows, junction_mask_bp):
    contig = row["contig"]
    start = int(row["feature_start"])
    end = int(row["feature_end"])
    reverse_oriented = bool(row.get("reverse_oriented_contig", 0))
    features = {}
    for window in windows:
        left, right = fetch_oriented_boundary_flanks(
            fasta_index,
            contig,
            start,
            end,
            window,
            reverse_oriented=reverse_oriented,
        )
        left, right = mask_junction_bases(left, right, n=junction_mask_bp)
        features.update(features_from_flanks(left, right, window))
    return features


def add_absent_window_features(row, fasta_index, windows, junction_mask_bp):
    contig = row["contig"]
    site = int(row["site"])
    reverse_oriented = bool(row.get("reverse_oriented_contig", 0))
    features = {}
    for window in windows:
        left, right = fetch_oriented_site_flanks(
            fasta_index,
            contig,
            site,
            window,
            reverse_oriented=reverse_oriented,
        )
        left, right = mask_junction_bases(left, right, n=junction_mask_bp)
        features.update(features_from_flanks(left, right, window))
    return features


def present_rows_from_features(feature_matrix, genotype_df, sample, excluded_contigs):
    features = pd.read_csv(feature_matrix, sep="\t", low_memory=False)
    present = features.loc[features["label"].eq(1)].copy()
    if excluded_contigs:
        present = present.loc[~present["contig"].isin(excluded_contigs)].copy()

    status = sample_status(genotype_df, sample)
    status = status.loc[status["presence"].eq(1)].copy()
    metadata_cols = [
        "ortholog_id",
        "sequence_id",
        "body_start",
        "body_end",
        "group1_pattern",
        "group1_present_count",
        "group1_absent_count",
        "group1_missing_count",
        "group1_callable_count",
        "group1_n_samples",
        "within_group_status",
        "within_group_identity",
    ]
    metadata_cols = [col for col in metadata_cols if col in status.columns]
    present = present.merge(
        status[metadata_cols].drop_duplicates("sequence_id"),
        on="sequence_id",
        how="left",
        suffixes=("", "_from_genotype_matrix"),
    )
    present = present.drop(columns=["body_start", "body_end"], errors="ignore")
    present["sample"] = sample
    present["class_source"] = f"{sample}_present_introner"
    present["sample_presence_status"] = "present"
    present["label"] = 1
    return present


def present_rows_from_genotype(
    genotype_df,
    sample,
    fasta_index,
    windows,
    boundary_support,
    junction_mask_bp,
    excluded_contigs,
    reverse_contigs,
):
    rows = genotype_df.loc[
        genotype_df["sample"].eq(sample) & genotype_df["presence"].eq(1)
    ].copy()
    rows = rows.loc[rows["contig"].notna()].copy()
    if excluded_contigs:
        rows = rows.loc[~rows["contig"].isin(excluded_contigs)].copy()

    rows["start"] = rows["start"].astype(int) + 100
    rows["end"] = rows["end"].astype(int) - 100
    rows = rows.loc[rows["start"].lt(rows["end"])].copy()
    present = rows[[
        col for col in [
            "contig",
            "start",
            "end",
            "sequence_id",
            "gene",
            "family",
            "ortholog_id",
            "sample",
            "group1_pattern",
            "group1_present_count",
            "group1_absent_count",
            "group1_missing_count",
            "group1_callable_count",
            "group1_n_samples",
            "within_group_status",
            "within_group_identity",
        ]
        if col in rows.columns
    ]].copy()
    present = apply_boundary_support(present, boundary_support)
    present["reverse_oriented_contig"] = present["contig"].isin(reverse_contigs).astype(int)

    max_window = max(windows)
    valid_rows = []
    feature_rows = []
    for row in present.to_dict("records"):
        contig = row["contig"]
        start = int(row["feature_start"])
        end = int(row["feature_end"])
        if contig not in fasta_index:
            continue
        if start - max_window < 0 or end + max_window > len(fasta_index[contig]):
            continue
        valid_rows.append(row)
        feature_rows.append(
            add_present_window_features(
                row,
                fasta_index,
                windows,
                junction_mask_bp=junction_mask_bp,
            )
        )

    out = pd.concat(
        [pd.DataFrame(valid_rows), pd.DataFrame(feature_rows)],
        axis=1,
    )
    out["class_source"] = f"{sample}_present_introner"
    out["sample_presence_status"] = "present"
    out["label"] = 1
    return out


def absent_rows(
    genotype_df,
    sample,
    fasta_index,
    windows,
    junction_mask_bp,
    max_absent_span,
    excluded_contigs,
    reverse_contigs,
):
    rows = genotype_df.loc[
        genotype_df["sample"].eq(sample) & genotype_df["presence"].eq(2)
    ].copy()
    rows = rows.loc[rows["contig"].notna()].copy()
    if excluded_contigs:
        rows = rows.loc[~rows["contig"].isin(excluded_contigs)].copy()

    rows["absent_span"] = rows["end"].astype(int) - rows["start"].astype(int)
    if max_absent_span is not None:
        rows = rows.loc[rows["absent_span"].le(max_absent_span)].copy()
    rows["site"] = ((rows["start"].astype(int) + rows["end"].astype(int)) // 2).astype(int)

    max_window = max(windows)
    metadata_rows = []
    feature_rows = []
    for row in rows.to_dict("records"):
        contig = row["contig"]
        site = int(row["site"])
        if pd.isna(contig) or contig not in fasta_index:
            continue
        if site - max_window < 0 or site + max_window > len(fasta_index[contig]):
            continue

        metadata_rows.append({
            "contig": contig,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "sequence_id": row.get("sequence_id", pd.NA),
            "gene": row.get("gene", pd.NA),
            "family": row.get("family", pd.NA),
            "feature_start": site,
            "feature_end": site,
            "site": site,
            "absent_span": int(row["absent_span"]),
            "ortholog_id": row.get("ortholog_id", pd.NA),
            "sample": sample,
            "class_source": f"{sample}_absent_introner_site",
            "sample_presence_status": "absent",
            "group1_pattern": row.get("group1_pattern", pd.NA),
            "group1_present_count": row.get("group1_present_count", pd.NA),
            "group1_absent_count": row.get("group1_absent_count", pd.NA),
            "group1_missing_count": row.get("group1_missing_count", pd.NA),
            "group1_callable_count": row.get("group1_callable_count", pd.NA),
            "group1_n_samples": row.get("group1_n_samples", pd.NA),
            "within_group_status": row.get("within_group_status", pd.NA),
            "within_group_identity": row.get("within_group_identity", pd.NA),
            "reverse_oriented_contig": int(contig in reverse_contigs),
            "label": 0,
        })
        feature_rows.append(
            add_absent_window_features(
                {"contig": contig, "site": site, "reverse_oriented_contig": int(contig in reverse_contigs)},
                fasta_index,
                windows,
                junction_mask_bp,
            )
        )

    return pd.concat(
        [pd.DataFrame(metadata_rows), pd.DataFrame(feature_rows)],
        axis=1,
    )


def summarize(df):
    rows = [
        {"metric": "n_rows", "value": len(df)},
        {"metric": "n_present", "value": int(df["label"].eq(1).sum())},
        {"metric": "n_absent", "value": int(df["label"].eq(0).sum())},
        {"metric": "n_genes", "value": df["gene"].nunique(dropna=True)},
        {"metric": "n_contigs", "value": df["contig"].nunique(dropna=True)},
    ]
    for label, group in df.groupby("label", dropna=False):
        rows.append({
            "metric": f"label_{label}_n_missing_gene",
            "value": int(group["gene"].isna().sum()),
        })
        rows.append({
            "metric": f"label_{label}_n_genes",
            "value": group["gene"].nunique(dropna=True),
        })
    absent = df.loc[df["label"].eq(0)]
    if "absent_span" in absent:
        rows.extend([
            {"metric": "absent_median_span", "value": absent["absent_span"].median()},
            {"metric": "absent_max_span", "value": absent["absent_span"].max()},
            {"metric": "absent_n_span_le_200", "value": int(absent["absent_span"].le(200).sum())},
            {"metric": "absent_n_span_le_225", "value": int(absent["absent_span"].le(225).sum())},
        ])
    return pd.DataFrame(rows)


def excluded_contigs(args):
    contigs = set(args.exclude_contig)
    if args.exclude_mating_type_contig:
        contigs.add(MATING_TYPE_CONTIG)
    return contigs


def main():
    args = parse_args()
    windows = parse_windows(args.windows)
    output = Path(args.output)
    summary = Path(args.summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)

    print("Loading genotype matrix...")
    genotype_df = pd.read_csv(args.genotype_matrix, sep="\t", low_memory=False)
    genotype_df["start"] = genotype_df["start"].astype(int)
    genotype_df["end"] = genotype_df["end"].astype(int)

    print("Loading FASTA...")
    fasta_index = SeqIO.index(args.fa, "fasta")
    boundary_support = read_boundary_support(
        args.introner_boundary_support,
        args.boundary_max_abs_delta,
    )
    excluded = excluded_contigs(args)
    reverse_contigs = read_reverse_oriented_contigs(
        args.mummer_coords,
        args.orientation_threshold,
    )
    if reverse_contigs:
        print(f"Orientation-normalizing {len(reverse_contigs)} reverse-oriented contigs.")

    print("Selecting present class...")
    if args.present_feature_matrix:
        present = present_rows_from_features(
            args.present_feature_matrix,
            genotype_df,
            args.sample,
            excluded,
        )
    else:
        present = present_rows_from_genotype(
            genotype_df,
            args.sample,
            fasta_index,
            windows,
            boundary_support,
            args.junction_mask_bp,
            excluded,
            reverse_contigs,
        )

    print("Building absent class...")
    absent = absent_rows(
        genotype_df,
        args.sample,
        fasta_index,
        windows,
        args.junction_mask_bp,
        args.max_absent_span,
        excluded,
        reverse_contigs,
    )

    out_df = pd.concat([present, absent], ignore_index=True, sort=False)

    print("Replacing gene annotations from GTF intervals...")
    out_df = out_df.rename(columns={"gene": "gene_from_genotype_matrix"})
    annotator = GtfGeneAnnotator(args.gtf)
    out_df = annotator.annotate_frame(
        out_df,
        start_col="feature_start",
        end_col="feature_end",
    )

    out_df.to_csv(output, sep="\t", index=False)
    summarize(out_df).to_csv(summary, sep="\t", index=False)
    print(f"Wrote {output}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()

"""Merge RNA-seq counts with locus-aware introner and gene features."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, help="Merged featureCounts CSV")
    parser.add_argument("--isoform-data", required=True, help="Locus-aware model-data TSV")
    parser.add_argument("--gc-content", required=True, help="Gene-by-strain GC-content CSV")
    parser.add_argument("--output", required=True, help="Output merged CSV")
    return parser.parse_args()


def replicate_to_strain(replicate):
    replicate = str(replicate)
    if replicate.startswith("834"):
        return "CCMP1545"
    if replicate.startswith("1614"):
        return "RCC1614"
    if replicate.startswith("1749"):
        return "RCC1749"
    return np.nan


def main():
    args = parse_args()

    counts = pd.read_csv(args.counts)
    if "Geneid" not in counts.columns:
        raise ValueError("Count matrix must contain a Geneid column")

    annotation_columns = {"Geneid", "Length", "Chr", "Start", "End", "Strand"}
    sample_columns = [column for column in counts.columns if column not in annotation_columns]
    if not sample_columns:
        raise ValueError("No RNA-seq replicate columns found in count matrix")

    library_sizes = counts[sample_columns].sum(axis=0)
    count_long = counts[["Geneid", *sample_columns]].melt(
        id_vars="Geneid", var_name="replicate", value_name="raw_count"
    ).rename(columns={"Geneid": "gene_id"})
    count_long["library_size"] = count_long["replicate"].map(library_sizes)
    count_long["log_library_size"] = np.log(count_long["library_size"])
    count_long["strain"] = count_long["replicate"].map(replicate_to_strain)
    count_long = count_long.dropna(subset=["strain"])

    isoform_data = pd.read_csv(args.isoform_data, sep="\t")
    feature_columns = [
        "common_gene_id",
        "strain",
        "reference_relative_gain_count",
        "reference_relative_loss_count",
        "current_introner_count",
        "log_gene_span",
        "n_isoforms",
        "n_extra_isoforms",
    ]
    missing = sorted(set(feature_columns) - set(isoform_data.columns))
    if missing:
        raise ValueError(f"Isoform model table is missing columns: {', '.join(missing)}")

    isoform_features = (
        isoform_data[feature_columns]
        .drop_duplicates()
        .rename(columns={"common_gene_id": "gene_id"})
    )
    duplicate_keys = isoform_features.duplicated(["gene_id", "strain"], keep=False)
    if duplicate_keys.any():
        raise ValueError("Isoform feature table has conflicting gene-by-strain rows")

    gc_content = pd.read_csv(args.gc_content)
    isoform_features = isoform_features.merge(
        gc_content, on=["gene_id", "strain"], how="left", validate="one_to_one"
    ).dropna(subset=["GC_content"])

    merged = count_long.merge(
        isoform_features,
        on=["gene_id", "strain"],
        how="inner",
        validate="many_to_one",
    )
    if merged.empty:
        raise ValueError("Expression counts and locus-aware features have no matching rows")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    print(
        f"Wrote {len(merged)} gene-replicate rows for "
        f"{merged['gene_id'].nunique()} genes to {output}"
    )


if __name__ == "__main__":
    main()

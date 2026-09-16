#!/usr/bin/env python3
"""Prepare an expression model table directly from per-replicate featureCounts files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


REPLICATE_TO_STRAIN = {
    "834_A1": "CCMP1545",
    "834_A2": "CCMP1545",
    "834_A3": "CCMP1545",
    "834_B1": "CCMP1545",
    "1614_A1": "RCC1614",
    "1614_A2": "RCC1614",
    "1614_A3": "RCC1614",
    "1614_A4": "RCC1614",
    "1749_A1": "RCC1749",
    "1749_A2": "RCC1749",
    "1749_A4": "RCC1749",
    "1749_B1": "RCC1749",
}


def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for key, quoted, bare in re.findall(r"([^;\s=]+)[=\s]+(?:\"([^\"]*)\"|([^;\s]+))", text):
        attrs[key] = quoted or bare
    return attrs


def reference_gene_table(gtf_path: Path) -> pd.DataFrame:
    bounds: dict[tuple[str, str], list[int]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attributes(fields[8])
            gene = attrs.get("gene_id") or attrs.get("gene")
            if not gene:
                continue
            key = (fields[0], gene)
            start, end = int(fields[3]), int(fields[4])
            if key not in bounds:
                bounds[key] = [start, end]
            else:
                bounds[key][0] = min(bounds[key][0], start)
                bounds[key][1] = max(bounds[key][1], end)
    return pd.DataFrame(
        [
            {"contig": contig, "gene_id": gene, "gene_start": start, "gene_end": end}
            for (contig, gene), (start, end) in bounds.items()
        ]
    )


def read_featurecounts(path: Path, replicate: str, strain: str) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t", skiprows=1)
    count_column = data.columns[-1]
    result = data[["Geneid", "Length", count_column]].rename(
        columns={"Geneid": "gene_id", "Length": "effective_exon_length", count_column: "raw_count"}
    )
    result["replicate"] = replicate
    result["strain"] = strain
    return result


def median_ratio_size_factors(counts: pd.DataFrame) -> tuple[pd.Series, int]:
    matrix = counts.pivot(index="gene_id", columns="replicate", values="raw_count")
    complete_positive = matrix.notna().all(axis=1) & (matrix > 0).all(axis=1)
    normalization_genes = matrix.loc[complete_positive]
    if len(normalization_genes) < 100:
        raise ValueError("Too few complete positive genes for median-ratio normalization")
    geometric_means = np.exp(np.log(normalization_genes).mean(axis=1))
    ratios = normalization_genes.div(geometric_means, axis=0)
    factors = ratios.median(axis=0)
    factors = factors / np.exp(np.log(factors).mean())
    return factors, len(normalization_genes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    count_dir = root / "results" / "expression" / "counts"
    frames = []
    for replicate, strain in REPLICATE_TO_STRAIN.items():
        frames.append(read_featurecounts(count_dir / f"{replicate}.counts.txt", replicate, strain))
    counts = pd.concat(frames, ignore_index=True)

    reference_genes = reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")
    mt_genes = set(
        reference_genes.loc[
            (reference_genes["contig"] == "CCMP1545#0#scaffold_2")
            & (reference_genes["gene_start"] <= 1730591)
            & (reference_genes["gene_end"] >= 49808),
            "gene_id",
        ]
    )
    reference_gene_ids = set(reference_genes["gene_id"]) - mt_genes
    counts = counts[counts["gene_id"].isin(reference_gene_ids)].copy()

    length_conflicts = (
        counts.groupby(["gene_id", "strain"])["effective_exon_length"].nunique().gt(1).sum()
    )
    if length_conflicts:
        raise ValueError(f"{length_conflicts} gene-strain pairs have replicate-specific exon lengths")

    size_factors, n_normalization_genes = median_ratio_size_factors(counts)
    library_sizes = counts.groupby("replicate")["raw_count"].sum()
    normalization = pd.DataFrame(
        {
            "replicate": list(REPLICATE_TO_STRAIN),
            "strain": [REPLICATE_TO_STRAIN[x] for x in REPLICATE_TO_STRAIN],
        }
    )
    normalization["median_ratio_size_factor"] = normalization["replicate"].map(size_factors)
    normalization["library_size"] = normalization["replicate"].map(library_sizes)
    normalization["library_size_factor"] = normalization["library_size"] / np.exp(
        np.log(normalization["library_size"]).mean()
    )
    normalization.to_csv(output_dir / "normalization_factors.tsv", sep="\t", index=False)

    counts["median_ratio_size_factor"] = counts["replicate"].map(size_factors)
    counts["library_size"] = counts["replicate"].map(library_sizes)
    counts["library_size_factor"] = counts["replicate"].map(
        normalization.set_index("replicate")["library_size_factor"]
    )

    introner_path = (
        root
        / "results"
        / "expression"
        / "isoform_analysis"
        / "turnover"
        / "introner_features_by_gene_strain.tsv"
    )
    features = pd.read_csv(introner_path, sep="\t").rename(columns={"common_gene_id": "gene_id"})
    feature_columns = [
        "gene_id",
        "strain",
        "reference_relative_gain_count",
        "reference_relative_loss_count",
        "current_introner_count",
        "uncallable_locus_count",
        "current_missing_locus_count",
        "gene_locus_count",
        "event_callable_fraction",
        "mating_type_gene",
    ]
    counts = counts.merge(
        features[feature_columns], on=["gene_id", "strain"], how="left", validate="many_to_one"
    )
    zero_columns = [
        "reference_relative_gain_count",
        "reference_relative_loss_count",
        "current_introner_count",
        "uncallable_locus_count",
        "current_missing_locus_count",
        "gene_locus_count",
    ]
    for column in zero_columns:
        counts[column] = counts[column].fillna(0).astype(int)
    counts["event_callable_fraction"] = counts["event_callable_fraction"].fillna(1.0)
    counts["mating_type_gene"] = counts["mating_type_gene"].eq(True)
    counts["complete_event_callability"] = counts["event_callable_fraction"].eq(1.0)

    gc = pd.read_csv(root / "results" / "expression" / "glm_modeling" / "gc_content_by_gene_strain.csv")
    counts = counts.merge(gc, on=["gene_id", "strain"], how="left", validate="many_to_one")
    counts = counts.dropna(subset=["GC_content"]).copy()
    counts["effective_exon_kb"] = counts["effective_exon_length"] / 1000.0
    counts["offset_median_ratio_length"] = np.log(
        counts["median_ratio_size_factor"] * counts["effective_exon_kb"]
    )
    counts["offset_library_length"] = np.log(
        counts["library_size_factor"] * counts["effective_exon_kb"]
    )
    counts["normalized_count_per_kb"] = counts["raw_count"] / (
        counts["median_ratio_size_factor"] * counts["effective_exon_kb"]
    )
    counts["log1p_normalized_count_per_kb"] = np.log1p(counts["normalized_count_per_kb"])

    counts.to_csv(output_dir / "expression_model_data.tsv", sep="\t", index=False)
    gene_strain = counts.drop_duplicates(["gene_id", "strain"])
    summary = [
        "Negative-binomial expression pilot data preparation",
        "===================================================",
        f"Gene-replicate rows: {len(counts)}",
        f"Genes: {counts['gene_id'].nunique()}",
        f"Gene-strain rows: {len(gene_strain)}",
        f"Genes used for median-ratio normalization: {n_normalization_genes}",
        f"Mating-type genes excluded: {len(mt_genes)}",
        f"Complete-callability gene-strain rows: {int(gene_strain['complete_event_callability'].sum())}",
        f"Gene-strain rows with gains: {int((gene_strain['reference_relative_gain_count'] > 0).sum())}",
        f"Gene-strain rows with losses: {int((gene_strain['reference_relative_loss_count'] > 0).sum())}",
        f"Gene-strain rows with current introners: {int((gene_strain['current_introner_count'] > 0).sum())}",
        "",
        "The expression universe is built from featureCounts and exact reference gene IDs;",
        "SQANTI isoform detection is not used as an inclusion criterion.",
    ]
    (output_dir / "data_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()

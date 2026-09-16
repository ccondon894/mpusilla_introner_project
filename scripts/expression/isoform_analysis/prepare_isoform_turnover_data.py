#!/usr/bin/env python3
"""Prepare locus-aware data for introner turnover and isoform-richness models."""

from __future__ import annotations

import argparse
import bisect
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
REFERENCE = "CCMP1545"
SQANTI_CATEGORIES = {"full-splice_match", "novel_in_catalog", "novel_not_in_catalog"}
CALL_LABELS = {1: "present", 2: "absent", 3: "missing"}


def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for key, quoted, bare in re.findall(r"([^;\s=]+)[=\s]+(?:\"([^\"]*)\"|([^;\s]+))", text):
        attrs[key] = quoted or bare
    return attrs


def read_gtf(path: Path) -> tuple[pd.DataFrame, set[str]]:
    gene_bounds: dict[tuple[str, str], list[int]] = {}
    transcripts: dict[tuple[str, str], set[str]] = defaultdict(set)
    transcript_exons: dict[tuple[str, str, str], int] = defaultdict(int)

    with path.open() as handle:
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
            contig = fields[0]
            start, end = int(fields[3]), int(fields[4])
            key = (contig, gene)
            if key not in gene_bounds:
                gene_bounds[key] = [start, end]
            else:
                gene_bounds[key][0] = min(gene_bounds[key][0], start)
                gene_bounds[key][1] = max(gene_bounds[key][1], end)
            transcript = attrs.get("transcript_id") or attrs.get("transcript")
            if transcript:
                transcripts[key].add(transcript)
                if fields[2] == "exon":
                    transcript_exons[(contig, gene, transcript)] += 1

    rows = []
    for (contig, gene), (start, end) in gene_bounds.items():
        txs = transcripts[(contig, gene)]
        exon_counts = [transcript_exons[(contig, gene, tx)] for tx in txs]
        rows.append(
            {
                "contig": contig,
                "native_gene_id": gene,
                "gene_start": start,
                "gene_end": end,
                "gene_span": end - start + 1,
                "annotated_transcript_count": len(txs),
                "max_annotated_exons": max(exon_counts, default=0),
            }
        )
    frame = pd.DataFrame(rows)
    return frame, set(frame["native_gene_id"])


class GeneIndex:
    def __init__(self, genes: pd.DataFrame):
        self.by_contig: dict[str, tuple[list[int], list[tuple[int, int, str]]]] = {}
        for contig, group in genes.groupby("contig", sort=False):
            records = sorted(
                zip(group["gene_start"], group["gene_end"], group["native_gene_id"]),
                key=lambda row: row[0],
            )
            self.by_contig[contig] = ([row[0] for row in records], records)

    def locate(self, contig: str, start: float, end: float) -> tuple[str | None, int, str]:
        if contig not in self.by_contig or pd.isna(start) or pd.isna(end):
            return None, 0, "unmapped"
        midpoint = (float(start) + float(end)) / 2.0 + 0.5
        starts, records = self.by_contig[contig]
        stop = bisect.bisect_right(starts, midpoint)
        candidates = [row for row in records[:stop] if row[1] >= midpoint]
        if not candidates:
            return None, 0, "unmapped"
        candidates.sort(key=lambda row: ((row[1] - row[0] + 1), row[2]))
        method = "midpoint_unique" if len(candidates) == 1 else "midpoint_smallest_span"
        return candidates[0][2], len(candidates), method


def shannon_effective(counts: pd.Series) -> float:
    values = counts.to_numpy(dtype=float)
    total = values.sum()
    if total <= 0:
        return float("nan")
    probabilities = values[values > 0] / total
    return float(np.exp(-(probabilities * np.log(probabilities)).sum()))


def event_label(reference_call: int, sample_call: int, sample: str) -> str:
    if sample == REFERENCE:
        return "reference_self"
    if reference_call == 3 or sample_call == 3:
        return "uncallable"
    if reference_call == 2 and sample_call == 1:
        return "reference_relative_gain"
    if reference_call == 1 and sample_call == 2:
        return "reference_relative_loss"
    if reference_call == 1 and sample_call == 1:
        return "stable_present"
    if reference_call == 2 and sample_call == 2:
        return "stable_absent"
    return "uncallable"


def audit_and_summarize_sqanti(
    sample: str,
    path: Path,
    sample_genes: set[str],
    common_genes: set[str],
) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(path, sep="\t")
    data = data[data["structural_category"].isin(SQANTI_CATEGORIES)].copy()
    data = data[data["associated_gene"].notna()].copy()
    data["native_gene_id"] = data["associated_gene"].astype(str)
    data["sample"] = sample
    data["in_sample_gtf"] = data["native_gene_id"].isin(sample_genes)
    data["common_gene_id"] = data["native_gene_id"].where(
        data["native_gene_id"].isin(common_genes)
    )
    data["crosswalk_status"] = np.select(
        [
            data["in_sample_gtf"] & data["common_gene_id"].notna(),
            data["in_sample_gtf"],
        ],
        ["exact_common_gene_id", "sample_gtf_only"],
        default="not_in_current_sample_gtf",
    )

    count_columns = [column for column in data.columns if re.fullmatch(r"\d+_[A-Z]\d+", column)]
    for column in count_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce").fillna(0)
    data["isoform_read_count"] = data[count_columns].sum(axis=1)
    data["exons"] = pd.to_numeric(data["exons"], errors="coerce").fillna(0).astype(int)

    audit = {
        "sample": sample,
        "sqanti_path": str(path),
        "qualifying_isoform_rows": len(data),
        "qualifying_unique_genes": data["native_gene_id"].nunique(),
        "genes_matching_sample_gtf": data.loc[data["in_sample_gtf"], "native_gene_id"].nunique(),
        "genes_with_exact_common_id": data["common_gene_id"].nunique(),
        "sample_gtf_gene_match_fraction": data.loc[data["in_sample_gtf"], "native_gene_id"].nunique()
        / max(data["native_gene_id"].nunique(), 1),
        "exact_common_gene_fraction": data["common_gene_id"].nunique()
        / max(data["native_gene_id"].nunique(), 1),
        "replicate_count_columns": ",".join(count_columns),
        "has_transcript_structure_coordinates": False,
    }

    crosswalk = (
        data[["sample", "native_gene_id", "common_gene_id", "in_sample_gtf", "crosswalk_status"]]
        .drop_duplicates()
        .sort_values(["sample", "native_gene_id"])
    )

    mapped = data[data["common_gene_id"].notna()].copy()
    richness_rows = []
    for common_gene, group in mapped.groupby("common_gene_id", sort=False):
        isoforms = group.drop_duplicates("isoform")
        read_counts = isoforms["isoform_read_count"]
        richness_rows.append(
            {
                "common_gene_id": common_gene,
                "strain": sample,
                "n_isoforms": isoforms["isoform"].nunique(),
                "n_extra_isoforms": max(isoforms["isoform"].nunique() - 1, 0),
                "has_alternative_isoform": int(isoforms["isoform"].nunique() > 1),
                "n_multi_exon_isoforms": isoforms.loc[isoforms["exons"] > 1, "isoform"].nunique(),
                "n_novel_splice_isoforms": isoforms.loc[
                    isoforms["structural_category"].isin({"novel_in_catalog", "novel_not_in_catalog"}),
                    "isoform",
                ].nunique(),
                "n_full_splice_match_isoforms": isoforms.loc[
                    isoforms["structural_category"] == "full-splice_match", "isoform"
                ].nunique(),
                "total_long_read_count": float(read_counts.sum()),
                "median_isoform_read_count": float(read_counts.median()),
                "n_isoforms_ge2_reads": int((read_counts >= 2).sum()),
                "n_isoforms_ge5_reads": int((read_counts >= 5).sum()),
                "read_weighted_effective_isoforms": shannon_effective(read_counts),
                "splice_chain_richness": np.nan,
                "splice_chain_status": "unavailable_without_r2c2_transcript_gtf",
            }
        )
    return data, audit, crosswalk, pd.DataFrame(richness_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--genotype-matrix", type=Path)
    parser.add_argument("--ccmp1545-gtf", type=Path)
    parser.add_argument("--rcc1614-gtf", type=Path)
    parser.add_argument("--rcc1749-gtf", type=Path)
    parser.add_argument("--ccmp1545-sqanti", type=Path)
    parser.add_argument("--rcc1614-sqanti", type=Path)
    parser.add_argument("--rcc1749-sqanti", type=Path)
    parser.add_argument("--mt-scaffold", default="CCMP1545#0#scaffold_2")
    parser.add_argument("--mt-start", type=int, default=49808)
    parser.add_argument("--mt-end", type=int, default=1730591)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "results" / "expression" / "isoform_analysis" / "turnover",
    )
    args = parser.parse_args()

    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    gtf_paths = {
        "CCMP1545": args.ccmp1545_gtf or root / "results" / "annotations" / "CCMP1545.gtf",
        "RCC1614": args.rcc1614_gtf or root / "results" / "annotations" / "RCC1614.gtf",
        "RCC1749": args.rcc1749_gtf or root / "results" / "annotations" / "RCC1749.gtf",
    }
    sqanti_paths = {
        "CCMP1545": args.ccmp1545_sqanti
        or root / "data" / "sqanti3_output_834" / "834_isoforms_classification.filtered.txt",
        "RCC1614": args.rcc1614_sqanti
        or root / "data" / "sqanti3_output_1614" / "1614_isoforms_classification.filtered.txt",
        "RCC1749": args.rcc1749_sqanti
        or root / "data" / "sqanti3_output_1749" / "1749_isoforms_classification.filtered.txt",
    }
    genotype_path = (
        args.genotype_matrix or root / "results" / "genotyping" / "genotype_matrix.final.tsv"
    )

    structures: dict[str, pd.DataFrame] = {}
    gene_sets: dict[str, set[str]] = {}
    for sample, path in gtf_paths.items():
        structures[sample], gene_sets[sample] = read_gtf(path)
        structures[sample]["strain"] = sample
    common_genes = gene_sets[REFERENCE]

    audits, crosswalks, richness_frames = [], [], []
    for sample in SAMPLES:
        _, audit, crosswalk, richness = audit_and_summarize_sqanti(
            sample, sqanti_paths[sample], gene_sets[sample], common_genes
        )
        audits.append(audit)
        crosswalks.append(crosswalk)
        richness_frames.append(richness)

    pd.DataFrame(audits).to_csv(output_dir / "sqanti_gtf_audit.tsv", sep="\t", index=False)
    pd.concat(crosswalks, ignore_index=True).to_csv(
        output_dir / "sample_gene_crosswalk.tsv", sep="\t", index=False
    )
    richness = pd.concat(richness_frames, ignore_index=True)
    richness.to_csv(output_dir / "isoform_richness_by_gene_strain.tsv", sep="\t", index=False)

    genotype = pd.read_csv(genotype_path, sep="\t", low_memory=False)
    genotype = genotype[genotype["sample"].isin(SAMPLES)].copy()
    genotype["presence"] = pd.to_numeric(genotype["presence"], errors="coerce").fillna(3).astype(int)
    indexes = {sample: GeneIndex(structures[sample]) for sample in SAMPLES}

    mapped_rows = []
    for row in genotype.itertuples(index=False):
        native_gene, n_candidates, method = indexes[row.sample].locate(
            row.contig, row.start, row.end
        )
        mapped_rows.append(
            {
                "ortholog_id": row.ortholog_id,
                "sample": row.sample,
                "presence": row.presence,
                "call_label": CALL_LABELS.get(row.presence, "missing"),
                "contig": row.contig,
                "start": row.start,
                "end": row.end,
                "native_gene_id": native_gene,
                "gene_candidate_count": n_candidates,
                "gene_mapping_method": method,
            }
        )
    locus_map = pd.DataFrame(mapped_rows)

    ortholog_common = {}
    ortholog_status = {}
    for ortholog, group in locus_map.groupby("ortholog_id", sort=False):
        candidates = sorted(set(group["native_gene_id"].dropna()) & common_genes)
        if len(candidates) == 1:
            ortholog_common[ortholog] = candidates[0]
            ortholog_status[ortholog] = "one_exact_common_gene"
        elif not candidates:
            ortholog_common[ortholog] = None
            ortholog_status[ortholog] = "no_exact_common_gene"
        else:
            ortholog_common[ortholog] = None
            ortholog_status[ortholog] = "conflicting_common_genes"
    locus_map["common_gene_id"] = locus_map["ortholog_id"].map(ortholog_common)
    locus_map["ortholog_gene_mapping_status"] = locus_map["ortholog_id"].map(ortholog_status)

    mt_genes = set(
        structures[REFERENCE].loc[
            (structures[REFERENCE]["contig"] == args.mt_scaffold)
            & (structures[REFERENCE]["gene_start"] <= args.mt_end)
            & (structures[REFERENCE]["gene_end"] >= args.mt_start),
            "native_gene_id",
        ]
    )
    locus_map["mating_type_gene"] = locus_map["common_gene_id"].isin(mt_genes)
    locus_map.to_csv(output_dir / "locus_gene_map.tsv", sep="\t", index=False)

    pivot = locus_map.pivot_table(
        index=["ortholog_id", "common_gene_id", "mating_type_gene"],
        columns="sample",
        values="presence",
        aggfunc="first",
    ).reset_index()
    native_gene_lookup = locus_map.set_index(["ortholog_id", "sample"])["native_gene_id"].to_dict()
    event_rows = []
    for row in pivot.itertuples(index=False):
        reference_call = int(getattr(row, REFERENCE))
        for sample in SAMPLES:
            sample_call = int(getattr(row, sample))
            label = event_label(reference_call, sample_call, sample)
            if label == "reference_relative_gain":
                event_gene_link_confident = (
                    native_gene_lookup.get((row.ortholog_id, sample)) == row.common_gene_id
                )
            elif label == "reference_relative_loss":
                event_gene_link_confident = (
                    native_gene_lookup.get((row.ortholog_id, REFERENCE)) == row.common_gene_id
                )
            else:
                event_gene_link_confident = True
            event_rows.append(
                {
                    "ortholog_id": row.ortholog_id,
                    "common_gene_id": row.common_gene_id,
                    "strain": sample,
                    "reference_call": reference_call,
                    "sample_call": sample_call,
                    "event_status": label,
                    "gain_count": int(label == "reference_relative_gain"),
                    "loss_count": int(label == "reference_relative_loss"),
                    "event_gene_link_confident": bool(event_gene_link_confident),
                    "current_present_count": int(sample_call == 1),
                    "current_missing_count": int(sample_call == 3),
                    "callable_comparison": int(label not in {"uncallable", "reference_self"}),
                    "mating_type_gene": row.mating_type_gene,
                }
            )
    events = pd.DataFrame(event_rows)
    events.to_csv(output_dir / "locus_level_reference_comparisons.tsv", sep="\t", index=False)

    def aggregate_events(group: pd.DataFrame) -> pd.Series:
        nonself = group["event_status"] != "reference_self"
        denominator = int(nonself.sum())
        return pd.Series(
            {
                "reference_relative_gain_count": int(
                    (group["gain_count"] * group["event_gene_link_confident"].astype(int)).sum()
                ),
                "reference_relative_loss_count": int(
                    (group["loss_count"] * group["event_gene_link_confident"].astype(int)).sum()
                ),
                "raw_reference_relative_gain_count": int(group["gain_count"].sum()),
                "raw_reference_relative_loss_count": int(group["loss_count"].sum()),
                "uncertain_event_gene_link_count": int(
                    (
                        ((group["gain_count"] == 1) | (group["loss_count"] == 1))
                        & ~group["event_gene_link_confident"]
                    ).sum()
                ),
                "current_introner_count": int(group["current_present_count"].sum()),
                "stable_present_count": int((group["event_status"] == "stable_present").sum()),
                "stable_absent_count": int((group["event_status"] == "stable_absent").sum()),
                "uncallable_locus_count": int((group["event_status"] == "uncallable").sum()),
                "current_missing_locus_count": int(group["current_missing_count"].sum()),
                "gene_locus_count": len(group),
                "event_callable_fraction": (
                    float(group.loc[nonself, "callable_comparison"].sum() / denominator)
                    if denominator else 1.0
                ),
                "mating_type_gene": bool(group["mating_type_gene"].any()),
            }
        )

    gene_events = (
        events.groupby(["common_gene_id", "strain"], sort=False)
        .apply(aggregate_events, include_groups=False)
        .reset_index()
    )
    gene_events.to_csv(output_dir / "introner_features_by_gene_strain.tsv", sep="\t", index=False)

    structure = pd.concat(structures.values(), ignore_index=True).rename(
        columns={"native_gene_id": "common_gene_id"}
    )
    structure = structure[structure["common_gene_id"].isin(common_genes)].copy()
    structure = structure.sort_values("gene_span").drop_duplicates(["common_gene_id", "strain"])

    model_data = richness.merge(gene_events, on=["common_gene_id", "strain"], how="left")
    zero_feature_columns = [
        "reference_relative_gain_count",
        "reference_relative_loss_count",
        "raw_reference_relative_gain_count",
        "raw_reference_relative_loss_count",
        "uncertain_event_gene_link_count",
        "current_introner_count",
        "stable_present_count",
        "stable_absent_count",
        "uncallable_locus_count",
        "current_missing_locus_count",
        "gene_locus_count",
    ]
    for column in zero_feature_columns:
        model_data[column] = model_data[column].fillna(0).astype(int)
    model_data["event_callable_fraction"] = model_data["event_callable_fraction"].fillna(1.0)
    model_data["mating_type_gene"] = (
        model_data["mating_type_gene"].where(model_data["mating_type_gene"].notna(), False).astype(bool)
    )
    model_data = model_data.merge(
        structure[
            [
                "common_gene_id",
                "strain",
                "gene_span",
                "annotated_transcript_count",
                "max_annotated_exons",
            ]
        ],
        on=["common_gene_id", "strain"],
        how="left",
    )
    model_data["log1p_long_read_count"] = np.log1p(model_data["total_long_read_count"])
    model_data["log_gene_span"] = np.log(model_data["gene_span"])
    model_data["high_detection_support"] = model_data["total_long_read_count"] >= 10
    model_data["complete_event_callability"] = model_data["event_callable_fraction"] == 1.0
    model_data = model_data[~model_data["mating_type_gene"]].copy()
    model_data.to_csv(output_dir / "isoform_introner_model_data.tsv", sep="\t", index=False)

    summary = [
        "Isoform turnover data preparation",
        "=======================================",
        f"Model rows: {len(model_data)}",
        f"Common genes: {model_data['common_gene_id'].nunique()}",
        f"Rows with >=10 long-read counts: {int(model_data['high_detection_support'].sum())}",
        f"Rows with complete event callability: {int(model_data['complete_event_callability'].sum())}",
        "",
        "Exact splice-chain richness was not computed because the imported SQANTI",
        "classification tables do not contain junction coordinates and the matching",
        "R2C2 transcript GTF files are not present in the project.",
    ]
    (output_dir / "data_preparation_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()

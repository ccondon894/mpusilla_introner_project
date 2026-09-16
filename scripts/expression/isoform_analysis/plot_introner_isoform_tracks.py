#!/usr/bin/env python3
"""Render a custom isoform-model figure for an introner-polymorphic gene."""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/scratch1/chris/tmp/.cache")

import pandas as pd


DEFAULT_GENE = "MicpuC2.estExt_fgenesh1_pg.C_170101.3.0.228"
DEFAULT_INTRONER_INTRON_TABLE = Path(
    "analysis/ccmp1545_gtf_liftover_test/current_genotype_introner_vs_annotated_intron_size.tsv"
)
GROUP2_SAMPLES = {"RCC1749", "RCC3052"}
HIGH_CONF_CATEGORIES = {
    "full-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog",
}
STRAIN_TO_SQANTI = {
    "CCMP1545": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/sqanti3_output_834/834_isoforms_classification.filtered.txt",
    "RCC1614": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/sqanti3_output_1614/1614_isoforms_classification.filtered.txt",
    "RCC1749": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/sqanti3_output_1749/1749_isoforms_classification.filtered.txt",
}
STRAIN_TO_R2C2_GTF = {
    "CCMP1545": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf",
    "RCC1614": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/09092025_1614_Isoforms.filtered.clean.sorted.gtf",
    "RCC1749": "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms/09092025_1749_Isoforms.filtered.clean.sorted.gtf",
}


@dataclass(frozen=True)
class Feature:
    contig: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attrs: dict[str, str]


@dataclass
class IsoformModel:
    sample: str
    transcript_id: str
    label: str
    strand: str
    exons: list[tuple[int, int]]
    cds_start: int | None
    cds_end: int | None
    category: str
    subcategory: str


@dataclass(frozen=True)
class LocusAnnotation:
    ortholog_id: str
    contig: str
    body_start: int
    body_end: int


@dataclass(frozen=True)
class LocusInterval:
    ortholog_id: str
    contig: str
    start: int
    end: int
    presence: int | None


@dataclass
class SamplePanel:
    sample: str
    contig: str
    origin: int
    region_start: int
    region_end: int
    transcript: Feature
    exons: list[Feature]
    isoforms: list[IsoformModel]


class MissingSampleGene(ValueError):
    """Raised when a requested gene is not annotated in a sample GTF."""


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def normalize_gene_id(gene_id: object) -> str:
    value = "" if pd.isna(gene_id) else str(gene_id)
    if value.endswith(".3.0.228"):
        value = value[: -len(".3.0.228")]
    if value.startswith("novelGene_") and value.endswith("_AS"):
        value = value[len("novelGene_") : -len("_AS")]
        if value.endswith(".3.0.228"):
            value = value[: -len(".3.0.228")]
    return value


def chromosome_label(contig: str) -> str:
    named_contig_match = re.search(r"(?:scaffold|intronerless_contig)_(\d+)", contig)
    if named_contig_match:
        return named_contig_match.group(1)

    parts = contig.split("#")
    if len(parts) >= 3 and parts[-1] == "0" and parts[-2].isdigit():
        return parts[-2]
    for part in reversed(parts):
        if part.isdigit():
            return part
    return contig


def read_gtf_features(gtf_path: Path, gene_id: str) -> list[Feature]:
    features: list[Feature] = []
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("gene_id") != gene_id:
                continue
            features.append(
                Feature(
                    contig=fields[0],
                    source=fields[1],
                    feature=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    score=fields[5],
                    strand=fields[6],
                    frame=fields[7],
                    attrs=attrs,
                )
            )
    return features


def read_gtf_features_normalized(gtf_path: Path, gene_id: str) -> list[Feature]:
    target = normalize_gene_id(gene_id)
    features: list[Feature] = []
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attrs(fields[8])
            gene_value = attrs.get("gene_id") or attrs.get("gene_name")
            if normalize_gene_id(gene_value) != target:
                continue
            features.append(
                Feature(
                    contig=fields[0],
                    source=fields[1],
                    feature=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    score=fields[5],
                    strand=fields[6],
                    frame=fields[7],
                    attrs=attrs,
                )
            )
    return features


def read_gene_model(gtf_path: Path, gene_id: str) -> tuple[list[Feature], list[Feature]]:
    features = read_gtf_features(gtf_path, gene_id)
    if not features:
        raise MissingSampleGene(f"No GTF features for {gene_id} in {gtf_path}")

    transcript_features = [f for f in features if f.feature == "transcript"]
    exon_features = [f for f in features if f.feature == "exon"]
    if not transcript_features or not exon_features:
        raise MissingSampleGene(
            f"Expected transcript and exon features for {gene_id} in {gtf_path}"
        )

    # Use the first transcript as the coordinate skeleton. Current annotations
    # for this locus have one miniprot transcript per sample.
    tx_id = transcript_features[0].attrs.get("transcript_id")
    tx_exons = [
        f for f in exon_features if f.attrs.get("transcript_id") == tx_id
    ] or exon_features
    tx_exons = sorted(tx_exons, key=lambda f: (f.start, f.end))
    return [transcript_features[0]], tx_exons


def group1_samples(genotype: pd.DataFrame) -> list[str]:
    samples = sorted(genotype["sample"].dropna().unique())
    return [sample for sample in samples if sample not in GROUP2_SAMPLES]


def read_gene_locus_annotations(
    path: Path,
    gene_id: str,
    min_body_in_intron_frac: float,
) -> dict[str, LocusAnnotation]:
    if not path.exists():
        raise SystemExit(f"Introner/intron overlap table not found: {path}")

    table = pd.read_csv(path, sep="\t")
    required = {"ortholog_id", "gene", "contig", "body_start", "body_end"}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(
            "Introner/intron overlap table is missing columns: "
            f"{', '.join(sorted(missing))}"
        )

    target = normalize_gene_id(gene_id)
    gene_match = table["gene"].map(normalize_gene_id).eq(target)
    if "intron_gene" in table.columns:
        gene_match = gene_match | table["intron_gene"].map(normalize_gene_id).eq(
            target
        )
    linked = table.loc[gene_match].copy()
    if "body_in_intron_frac" in linked.columns:
        linked = linked.loc[
            linked["body_in_intron_frac"].fillna(0) >= min_body_in_intron_frac
        ].copy()
    if linked.empty:
        return {}

    linked = linked.sort_values(["contig", "body_start", "body_end", "ortholog_id"])
    annotations: dict[str, LocusAnnotation] = {}
    for _, row in linked.iterrows():
        oid = str(row["ortholog_id"])
        if oid in annotations:
            continue
        annotations[oid] = LocusAnnotation(
            ortholog_id=oid,
            contig=str(row["contig"]),
            body_start=int(row["body_start"]),
            body_end=int(row["body_end"]),
        )
    return annotations


def group1_polymorphic_loci(
    genotype: pd.DataFrame,
    linked_loci: list[str],
    group1: list[str],
    require_complete: bool,
) -> list[str]:
    result: list[str] = []
    for oid in linked_loci:
        calls = (
            genotype.loc[
                genotype["ortholog_id"].eq(oid) & genotype["sample"].isin(group1),
                ["sample", "presence"],
            ]
            .drop_duplicates("sample")
            .set_index("sample")["presence"]
            .reindex(group1)
        )
        if require_complete and not calls.isin([1, 2]).all():
            continue
        if calls.eq(1).any() and calls.eq(2).any():
            result.append(oid)
    return result


def group1_callable_loci(
    genotype: pd.DataFrame,
    linked_loci: list[str],
    group1: list[str],
) -> list[str]:
    result: list[str] = []
    for oid in linked_loci:
        calls = genotype.loc[
            genotype["ortholog_id"].eq(oid) & genotype["sample"].isin(group1),
            "presence",
        ]
        if calls.isin([1, 2]).any():
            result.append(oid)
    return result


def select_state_loci(
    genotype: pd.DataFrame,
    linked_loci: list[str],
    group1: list[str],
) -> tuple[list[str], str]:
    complete = group1_polymorphic_loci(
        genotype, linked_loci, group1, require_complete=True
    )
    if complete:
        return complete, "complete_group1_polymorphic"

    incomplete = group1_polymorphic_loci(
        genotype, linked_loci, group1, require_complete=False
    )
    if incomplete:
        return incomplete, "group1_polymorphic_allowing_missing"

    callable_loci = group1_callable_loci(genotype, linked_loci, group1)
    if callable_loci:
        return callable_loci, "gene_linked_callable"

    return [], "none"


def observed_group1_states(
    genotype: pd.DataFrame, loci: list[str], group1: list[str]
) -> pd.DataFrame:
    matrix = (
        genotype.loc[genotype["ortholog_id"].isin(loci)]
        .pivot(index="sample", columns="ortholog_id", values="presence")
        .reindex(group1)
        .reindex(columns=loci)
    )
    labels = matrix.apply(
        lambda col: col.map({1: "P", 2: "A", 3: "M"}).fillna("M")
    )
    states = labels.apply(lambda row: "".join(row.values), axis=1)
    rows = []
    for state, sample_index in states.groupby(states).groups.items():
        samples = sorted(sample_index)
        rows.append(
            {
                "state": state,
                "n_samples": len(samples),
                "samples": ",".join(samples),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_samples", "state"], ascending=[False, True]
    )


def read_sqanti_table(path: Path, gene_id: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    sqanti = pd.read_csv(path, sep="\t")
    if sqanti.empty:
        return sqanti
    sqanti["norm_gene"] = sqanti["associated_gene"].map(normalize_gene_id)
    sub = sqanti[
        sqanti["norm_gene"].eq(normalize_gene_id(gene_id))
        & sqanti["structural_category"].isin(HIGH_CONF_CATEGORIES)
    ].copy()
    return sub


def has_coding_exon(
    exons: list[tuple[int, int]],
    cds_start: int | None,
    cds_end: int | None,
) -> bool:
    if cds_start is None or cds_end is None:
        return False
    return any(max(start, cds_start) <= min(end, cds_end) for start, end in exons)


def read_r2c2_isoforms(
    gtf_path: Path,
    sqanti_path: Path,
    gene_id: str,
    sample: str,
) -> list[IsoformModel]:
    features = read_gtf_features_normalized(gtf_path, gene_id)
    if not features:
        return []
    exons_by_tx: dict[str, list[Feature]] = defaultdict(list)
    tx_strand: dict[str, str] = {}
    for feature in features:
        tx = feature.attrs.get("transcript_id")
        if not tx:
            continue
        if feature.feature == "transcript":
            tx_strand[tx] = feature.strand
        elif feature.feature == "exon":
            exons_by_tx[tx].append(feature)

    sqanti = read_sqanti_table(sqanti_path, gene_id)
    sqanti_by_tx = sqanti.set_index("isoform").to_dict("index") if not sqanti.empty else {}

    isoforms: list[IsoformModel] = []
    for tx, exons in sorted(exons_by_tx.items()):
        exons_sorted = sorted((exon.start, exon.end) for exon in exons)
        row = sqanti_by_tx.get(tx)
        if row is None:
            continue
        start_raw = row.get("CDS_genomic_start")
        end_raw = row.get("CDS_genomic_end")
        cds_start = cds_end = None
        if start_raw is not None and end_raw is not None and not pd.isna(start_raw) and not pd.isna(end_raw):
            cds_start = int(min(start_raw, end_raw))
            cds_end = int(max(start_raw, end_raw))
        if not has_coding_exon(exons_sorted, cds_start, cds_end):
            continue
        subcategory = str(row.get("subcategory", ""))
        category = str(row.get("structural_category", ""))
        retained = "intron_retention" in subcategory.lower()
        label_suffix = "IR" if retained else category.replace("_", " ")
        label = f"{tx} ({label_suffix})"
        isoforms.append(
            IsoformModel(
                sample=sample,
                transcript_id=tx,
                label=label,
                strand=tx_strand.get(tx, exons[0].strand),
                exons=exons_sorted,
                cds_start=cds_start,
                cds_end=cds_end,
                category=category,
                subcategory=subcategory,
            )
        )
    return isoforms


def merged_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
            continue
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def introns_from_intervals(
    intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    introns: list[tuple[int, int]] = []
    for left, right in zip(sorted(intervals), sorted(intervals)[1:]):
        start = left[1]
        end = right[0]
        if start < end:
            introns.append((start, end))
    return introns


def transform_pos(pos: int, deletions: list[tuple[int, int]]) -> int:
    offset = 0
    for start, end in deletions:
        if pos >= end:
            offset += end - start
            continue
        if pos > start:
            return start - offset
        break
    return pos - offset


def transform_intervals(
    intervals: list[tuple[int, int]],
    deletions: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    if not deletions:
        return merged_intervals(intervals)

    transformed: list[tuple[int, int]] = []
    for interval_start, interval_end in intervals:
        pieces = [(interval_start, interval_end)]
        for del_start, del_end in deletions:
            next_pieces: list[tuple[int, int]] = []
            for piece_start, piece_end in pieces:
                if del_end <= piece_start or del_start >= piece_end:
                    next_pieces.append((piece_start, piece_end))
                    continue
                if piece_start < del_start:
                    next_pieces.append((piece_start, del_start))
                if del_end < piece_end:
                    next_pieces.append((del_end, piece_end))
            pieces = next_pieces
        for piece_start, piece_end in pieces:
            if piece_start >= piece_end:
                continue
            transformed.append(
                (transform_pos(piece_start, deletions), transform_pos(piece_end, deletions))
            )
    return merged_intervals(transformed)


def local_isoform(isoform: IsoformModel, origin: int) -> IsoformModel:
    cds_start = isoform.cds_start - origin if isoform.cds_start is not None else None
    cds_end = isoform.cds_end - origin if isoform.cds_end is not None else None
    return IsoformModel(
        sample=isoform.sample,
        transcript_id=isoform.transcript_id,
        label=isoform.label,
        strand=isoform.strand,
        exons=merged_intervals(
            [(start - origin, end - origin) for start, end in isoform.exons]
        ),
        cds_start=cds_start,
        cds_end=cds_end,
        category=isoform.category,
        subcategory=isoform.subcategory,
    )


def corrected_isoform(isoform: IsoformModel, deletions: list[tuple[int, int]]) -> IsoformModel:
    if not deletions:
        return isoform
    cds_start = (
        transform_pos(isoform.cds_start, deletions)
        if isoform.cds_start is not None
        else None
    )
    cds_end = (
        transform_pos(isoform.cds_end, deletions)
        if isoform.cds_end is not None
        else None
    )
    return IsoformModel(
        sample=isoform.sample,
        transcript_id=isoform.transcript_id,
        label=isoform.label,
        strand=isoform.strand,
        exons=transform_intervals(isoform.exons, deletions),
        cds_start=cds_start,
        cds_end=cds_end,
        category=isoform.category,
        subcategory=isoform.subcategory,
    )


def gene_model_isoform(
    sample: str,
    transcript: Feature,
    exons: list[Feature],
    origin: int,
    deletions: list[tuple[int, int]],
) -> IsoformModel:
    local_exons = [(exon.start - origin, exon.end - origin) for exon in exons]
    return IsoformModel(
        sample=sample,
        transcript_id=f"{sample} gene",
        label="gene",
        strand=transcript.strand,
        exons=transform_intervals(local_exons, deletions),
        cds_start=transform_pos(min(exon.start for exon in exons) - origin, deletions),
        cds_end=transform_pos(max(exon.end for exon in exons) - origin, deletions),
        category="gene",
        subcategory="",
    )


def sample_locus_intervals(
    genotype: pd.DataFrame,
    sample: str,
    loci: list[str],
    annotations: dict[str, LocusAnnotation],
    genotype_flank_length: int,
) -> dict[str, LocusInterval]:
    intervals: dict[str, LocusInterval] = {}
    sample_rows = (
        genotype.loc[
            genotype["sample"].eq(sample) & genotype["ortholog_id"].isin(loci)
        ]
        .drop_duplicates("ortholog_id")
        .set_index("ortholog_id")
    )
    for oid in loci:
        row = sample_rows.loc[oid] if oid in sample_rows.index else None
        presence = (
            int(row["presence"])
            if row is not None and not pd.isna(row["presence"])
            else None
        )
        if sample == "CCMP1545" and oid in annotations:
            annotation = annotations[oid]
            intervals[oid] = LocusInterval(
                ortholog_id=oid,
                contig=annotation.contig,
                start=annotation.body_start,
                end=annotation.body_end,
                presence=presence,
            )
            continue
        if row is None:
            continue
        body_start = int(row["start"]) + genotype_flank_length
        body_end = int(row["end"]) - genotype_flank_length
        if body_start > body_end:
            midpoint = (int(row["start"]) + int(row["end"])) // 2
            body_start = midpoint
            body_end = midpoint + 1
        intervals[oid] = LocusInterval(
            ortholog_id=oid,
            contig=str(row["contig"]),
            start=body_start,
            end=body_end,
            presence=presence,
        )
    return intervals


def build_sample_panel(
    sample: str,
    gene_id: str,
    annotation_dir: Path,
    flank: int,
) -> SamplePanel:
    gtf_path = annotation_dir / f"{sample}.gtf"
    if not gtf_path.exists():
        raise SystemExit(f"Annotation GTF not found: {gtf_path}")
    transcript_features, exons = read_gene_model(gtf_path, gene_id)
    transcript = transcript_features[0]
    origin = min(exon.start for exon in exons)
    region_start = -flank
    region_end = max(exon.end for exon in exons) - origin + flank
    isoforms = read_r2c2_isoforms(
        Path(STRAIN_TO_R2C2_GTF[sample]),
        Path(STRAIN_TO_SQANTI[sample]),
        gene_id,
        sample,
    )
    return SamplePanel(
        sample=sample,
        contig=transcript.contig,
        origin=origin,
        region_start=region_start,
        region_end=region_end,
        transcript=transcript,
        exons=exons,
        isoforms=[local_isoform(isoform, origin) for isoform in isoforms],
    )


def absent_artifact_deletions(
    panel: SamplePanel,
    intervals: dict[str, LocusInterval],
) -> list[tuple[int, int]]:
    deletions: list[tuple[int, int]] = []
    for interval in intervals.values():
        if interval.presence != 2 or interval.contig != panel.contig:
            continue
        body_start = interval.start - panel.origin
        body_end = interval.end - panel.origin
        host_intron = host_intron_interval(panel, body_start, body_end)
        deletions.append(host_intron or (body_start, body_end))
    return merged_intervals(deletions)


def host_intron_interval(
    panel: SamplePanel,
    local_start: int,
    local_end: int,
) -> tuple[int, int] | None:
    local_gene_exons = [
        (exon.start - panel.origin, exon.end - panel.origin) for exon in panel.exons
    ]
    for intron_start, intron_end in introns_from_intervals(local_gene_exons):
        if interval_overlap(local_start, local_end, intron_start, intron_end):
            return intron_start, intron_end
    return None


def interval_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int
) -> tuple[int, int] | None:
    start = max(a_start, b_start)
    end = min(a_end, b_end)
    if start <= end:
        return start, end
    return None


def draw_segment(
    ax,
    start: int,
    end: int,
    y: float,
    height: float,
    color: str,
    edge: str = "#555555",
) -> None:
    from matplotlib.patches import Rectangle

    ax.add_patch(
        Rectangle(
            (start, y - height / 2),
            max(1, end - start),
            height,
            facecolor=color,
            edgecolor=edge,
            linewidth=0.45,
            zorder=3,
        )
    )


def draw_isoform(
    ax,
    isoform: IsoformModel,
    y: float,
    color: str,
    label_size: int,
    display_label: str | None = None,
) -> None:
    tx_start = min(start for start, _ in isoform.exons)
    tx_end = max(end for _, end in isoform.exons)
    ax.plot([tx_start, tx_end], [y, y], color="#1f1f1f", linewidth=0.55, zorder=1)

    for start, end in isoform.exons:
        cds_piece = None
        if isoform.cds_start is not None and isoform.cds_end is not None:
            cds_piece = interval_overlap(start, end, isoform.cds_start, isoform.cds_end)
        utr_pieces: list[tuple[int, int]] = []
        if cds_piece is None:
            utr_pieces.append((start, end))
        else:
            if start < cds_piece[0]:
                utr_pieces.append((start, cds_piece[0]))
            if cds_piece[1] < end:
                utr_pieces.append((cds_piece[1], end))
            draw_segment(ax, cds_piece[0], cds_piece[1], y, 0.48, color)
        for utr_start, utr_end in utr_pieces:
            draw_segment(ax, utr_start, utr_end, y, 0.22, "#F4A51C", "#B77A12")
    ax.text(
        tx_start - 35,
        y,
        display_label or isoform.transcript_id,
        ha="right",
        va="center",
        fontsize=label_size,
    )


def write_custom_isoform_plot(
    paths: list[Path],
    gene_id: str,
    panels: list[SamplePanel],
    loci: list[str],
    sample_intervals: dict[str, dict[str, LocusInterval]],
    label_size: int,
    fig_width: float,
    fig_height: float | None,
    tick_size: float,
    axis_label_size: float,
    title_size: float,
    dpi: int,
    introner_style: str,
    introner_color: str,
) -> dict[str, int]:
    import matplotlib.pyplot as plt

    row_counts = [1 + len(panel.isoforms) for panel in panels]
    if fig_height is None:
        fig_height = max(3.2, 0.36 * sum(row_counts) + 1.5)
    fig, axes = plt.subplots(
        nrows=len(panels),
        ncols=1,
        figsize=(fig_width, fig_height),
        height_ratios=row_counts,
        squeeze=False,
    )

    sample_colors = {
        "CCMP1545": "#D95F59",
        "RCC1614": "#6A9D73",
        "RCC1749": "#7C6BB0",
    }
    presence_colors = {
        1: "#F4C542",
        2: "#9A9A9A",
        3: "#D6D6D6",
        None: "#D6D6D6",
    }
    presence_labels = {1: "Present", 2: "Absent", 3: "Missing", None: "Missing"}

    isoform_counts: dict[str, int] = {}
    for ax, panel, row_count in zip(axes[:, 0], panels, row_counts):
        deletions = absent_artifact_deletions(
            panel, sample_intervals.get(panel.sample, {})
        )
        corrected_region_start = transform_pos(panel.region_start, deletions)
        corrected_region_end = transform_pos(panel.region_end, deletions)
        ax.set_xlim(corrected_region_start, corrected_region_end)
        ax.set_ylim(-0.8, row_count + 1.15)
        ax.set_yticks([])
        ax.tick_params(axis="x", labelsize=tick_size, length=3, width=0.8)

        gene_y = row_count - 0.35
        y = gene_y
        gene_model = gene_model_isoform(
            panel.sample, panel.transcript, panel.exons, panel.origin, deletions
        )
        draw_isoform(
            ax,
            gene_model,
            y,
            "#5B8CC0",
            label_size,
            display_label="gene",
        )
        ax.axhline(y - 0.55, color="#C8C8C8", linewidth=0.45, zorder=0)

        for isoform_number, isoform in enumerate(panel.isoforms, start=1):
            y -= 1
            plot_isoform = corrected_isoform(isoform, deletions)
            if not plot_isoform.exons or not has_coding_exon(
                plot_isoform.exons, plot_isoform.cds_start, plot_isoform.cds_end
            ):
                continue
            draw_isoform(
                ax,
                plot_isoform,
                y,
                sample_colors.get(panel.sample, "#6C6C6C"),
                label_size,
                display_label=f"Isoform {isoform_number}",
            )
        isoform_counts[panel.sample] = len(panel.isoforms)

        for locus_index, oid in enumerate(loci):
            interval = sample_intervals.get(panel.sample, {}).get(oid)
            if interval is None or interval.contig != panel.contig:
                continue
            local_start = interval.start - panel.origin
            local_end = interval.end - panel.origin
            color = presence_colors.get(interval.presence, "#D6D6D6")
            if interval.presence == 2:
                collapsed_x = transform_pos(local_start, deletions)
                ax.vlines(
                    x=collapsed_x,
                    ymin=-0.8,
                    ymax=row_count + 0.08,
                    color="#6F6F6F",
                    linestyle="--",
                    linewidth=0.7,
                    zorder=2,
                )
                label_x = collapsed_x
            else:
                span_start, span_end = local_start, local_end
                if interval.presence == 1:
                    span_start, span_end = host_intron_interval(
                        panel, local_start, local_end
                    ) or (local_start, local_end)
                plot_start = transform_pos(span_start, deletions)
                plot_end = transform_pos(span_end, deletions)
                if introner_style == "highlight":
                    ax.axvspan(
                        plot_start,
                        plot_end,
                        color=color,
                        alpha=0.2,
                        linewidth=0,
                    )
                else:
                    ax.plot(
                        [plot_start, plot_end],
                        [gene_y, gene_y],
                        color=introner_color,
                        linewidth=4.0,
                        solid_capstyle="butt",
                        zorder=5,
                    )
                label_x = (plot_start + plot_end) / 2
            locus_letter = chr(ord("A") + locus_index)
            label = (
                f"introner {locus_letter}\n"
                f"{presence_labels.get(interval.presence, 'Missing')}"
            )
            ax.text(
                label_x,
                row_count + 0.25,
                label,
                ha="center",
                va="bottom",
                fontsize=max(5, label_size - 1),
            )

        ax.set_xlabel(
            f"{panel.sample} chromosome {chromosome_label(panel.contig)}; "
            "bp from gene start",
            fontsize=axis_label_size,
        )
        for spine in ax.spines.values():
            spine.set_linewidth(0.7)

    fig.suptitle(gene_id, x=0.01, ha="left", fontsize=title_size, weight="bold")
    fig.tight_layout()
    for path in paths:
        fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return isoform_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot introner candidate gene states as an isoform-model panel."
    )
    parser.add_argument("--gene", default=DEFAULT_GENE)
    parser.add_argument("--sample", default="CCMP1545")
    parser.add_argument(
        "--genotype-matrix",
        type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
    )
    parser.add_argument(
        "--annotation-dir", type=Path, default=Path("results/annotations")
    )
    parser.add_argument(
        "--introner-intron-table",
        type=Path,
        default=DEFAULT_INTRONER_INTRON_TABLE,
        help=(
            "Gene-to-introner overlap table used for locus lookup and reference "
            "coordinates. The genotype matrix gene column is not used for this."
        ),
    )
    parser.add_argument(
        "--min-body-in-intron-frac",
        type=float,
        default=0.95,
        help="Minimum introner body fraction overlapping an annotated intron.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/figures/introner_isoform_tracks"),
    )
    parser.add_argument("--flank", type=int, default=400)
    parser.add_argument(
        "--genotype-flank-length",
        type=int,
        default=100,
        help="Flank length to strip from sample-specific genotype matrix loci.",
    )
    parser.add_argument("--top-sample", default="CCMP1545")
    parser.add_argument("--bottom-sample", default="RCC1614")
    parser.add_argument("--third-sample", default="RCC1749")
    parser.add_argument("--label-size", type=int, default=5)
    parser.add_argument(
        "--fig-width",
        type=float,
        default=12.0,
        help="Figure width in inches.",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        help="Optional fixed figure height in inches; otherwise determined by row count.",
    )
    parser.add_argument("--tick-size", type=float, default=7.0)
    parser.add_argument("--axis-label-size", type=float, default=8.0)
    parser.add_argument("--title-size", type=float, default=9.0)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument(
        "--introner-style",
        choices=["highlight", "gene-bar"],
        default="highlight",
        help=(
            "Draw present introners as full-height shaded spans or as compact "
            "bars on the gene-model track."
        ),
    )
    parser.add_argument(
        "--introner-color",
        default="#7C6BB0",
        help="Color used for present introners when --introner-style=gene-bar.",
    )
    parser.add_argument("--format", choices=["pdf", "png", "both"], default="both")
    parser.add_argument(
        "--output-prefix",
        help=(
            "Optional output filename stem, without .pdf/.png. Defaults to "
            "<sample>.<gene>.isoform_model_tracks."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    genotype = pd.read_csv(args.genotype_matrix, sep="\t")

    g1 = group1_samples(genotype)
    locus_annotations = read_gene_locus_annotations(
        args.introner_intron_table,
        args.gene,
        args.min_body_in_intron_frac,
    )
    linked_loci = list(locus_annotations)
    if not linked_loci:
        raise SystemExit(
            "No clean introner loci found for "
            f"{args.gene} in {args.introner_intron_table}"
        )
    loci, locus_mode = select_state_loci(genotype, linked_loci, g1)
    if not loci:
        raise SystemExit(f"No Group-1 callable loci found for {args.gene}")
    locus_annotations = {
        oid: annotation
        for oid, annotation in locus_annotations.items()
        if oid in set(loci)
    }
    states = observed_group1_states(genotype, loci, g1)
    requested_panel_samples = [args.top_sample, args.bottom_sample, args.third_sample]
    panels: list[SamplePanel] = []
    skipped_panels: list[str] = []
    for sample in requested_panel_samples:
        try:
            panels.append(
                build_sample_panel(sample, args.gene, args.annotation_dir, args.flank)
            )
        except MissingSampleGene as exc:
            skipped_panels.append(f"{sample}: {exc}")
    if not panels:
        skipped_text = "; ".join(skipped_panels) if skipped_panels else "none"
        raise SystemExit(
            f"No sample panels could be built for {args.gene}; skipped: {skipped_text}"
        )

    panel_samples = [panel.sample for panel in panels]
    sample_intervals = {
        sample: sample_locus_intervals(
            genotype,
            sample,
            loci,
            locus_annotations,
            args.genotype_flank_length,
        )
        for sample in panel_samples
    }

    output_formats = ["pdf", "png"] if args.format == "both" else [args.format]
    output_prefix = args.output_prefix or f"{args.sample}.{args.gene}.isoform_model_tracks"
    custom_plots = [
        args.outdir / f"{output_prefix}.{fmt}"
        for fmt in output_formats
    ]
    isoform_counts = write_custom_isoform_plot(
        custom_plots,
        args.gene,
        panels,
        loci,
        sample_intervals,
        args.label_size,
        args.fig_width,
        args.fig_height,
        args.tick_size,
        args.axis_label_size,
        args.title_size,
        args.dpi,
        args.introner_style,
        args.introner_color,
    )
    print(f"Gene: {args.gene}")
    print(f"Output prefix sample: {args.sample}")
    print(f"Output filename prefix: {output_prefix}")
    for skipped in skipped_panels:
        print(f"Skipped panel: {skipped}")
    for panel in panels:
        print(
            "Panel: "
            f"{panel.sample} {panel.contig}:"
            f"{panel.origin + panel.region_start}-{panel.origin + panel.region_end}"
        )
    print(f"Gene-linked clean loci: {','.join(linked_loci)}")
    print(f"State locus selection mode: {locus_mode}")
    print(f"State loci order: {','.join(loci)}")
    print(f"Observed Group-1 states: {len(states)}")
    print(states.to_string(index=False))
    counts_text = ", ".join(
        f"{sample}={isoform_counts.get(sample, 0)}" for sample in panel_samples
    )
    print(f"Custom plot R2C2 isoforms: {counts_text}")
    for custom_plot in custom_plots:
        print(f"Wrote plot: {custom_plot}")


if __name__ == "__main__":
    main()

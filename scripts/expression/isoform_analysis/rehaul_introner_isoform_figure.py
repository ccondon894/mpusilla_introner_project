#!/usr/bin/env python3
"""Plot CCMP1545/RCC1749 introner-polymorphic isoform-rich genes.

Each strain is drawn in its own transcript-oriented, gene-local coordinate
system. Candidate genes are supplied explicitly with repeatable ``--gene``
arguments, and callable introner differences are obtained from the complete
spatial locus-to-gene map.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/scratch1/chris/tmp/.cache")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


PROJECT_ROOT = Path(__file__).resolve().parents[3]
ISOFORM_SCRIPT_DIR = PROJECT_ROOT / "scripts/expression/isoform_analysis"
sys.path.insert(0, str(ISOFORM_SCRIPT_DIR))

import plot_introner_isoform_tracks as base


SAMPLES = ("CCMP1545", "RCC1749")
HIGH_CONF_CATEGORIES = {
    "full-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog",
}
SAMPLE_COLORS = {
    "CCMP1545": "#0072B2",
    "RCC1749": "#D55E00",
}
GENE_COLORS = {
    "CCMP1545": "#7DBFE0",
    "RCC1749": "#F0A868",
}
UTR_COLOR = "#999999"
UTR_EDGE = "#666666"
INTRONER_COLOR = "#7B4F7E"
NMD_COLOR = "#9E2A2B"


@dataclass(frozen=True)
class CandidateLocus:
    ortholog_id: str
    family: int
    ccmp_call: int
    rcc1749_call: int


@dataclass
class Candidate:
    gene_id: str
    loci: list[CandidateLocus]
    panels: dict[str, base.SamplePanel]
    nmd_isoforms: dict[str, set[str]]


@dataclass(frozen=True)
class PanelGeometry:
    start: int
    end: int
    strand: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot explicit CCMP1545/RCC1749 introner/isoform candidates."
    )
    parser.add_argument(
        "--genotype-matrix",
        type=Path,
        default=PROJECT_ROOT / "results/genotyping/genotype_matrix.final.tsv",
    )
    parser.add_argument(
        "--locus-gene-map",
        type=Path,
        default=PROJECT_ROOT
        / "results/expression/isoform_analysis/turnover/locus_gene_map.tsv",
        help="Spatial introner-locus to orthologous-gene mapping.",
    )
    parser.add_argument(
        "--annotation-dir",
        type=Path,
        default=PROJECT_ROOT / "results/annotations",
    )
    parser.add_argument(
        "--nmd-table",
        type=Path,
        default=PROJECT_ROOT / "results/expression/sqanti3/parsed_sqanti3_data.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT
        / "results/figures/introner_isoform_tracks_ccmp1545_rcc1749_candidates",
    )
    parser.add_argument(
        "--gene",
        action="append",
        required=True,
        help="Candidate gene to plot; may be supplied multiple times.",
    )
    parser.add_argument(
        "--exclude-locus",
        action="append",
        default=[],
        help=(
            "Orthologous introner locus to omit from every strain track; "
            "may be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--hide-locus-in-sample",
        action="append",
        default=[],
        metavar="SAMPLE:ORTHOLOG_ID",
        help=(
            "Suppress one locus marker and label in one sample while retaining "
            "its annotation in the other sample; may be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--annotation-x-offset",
        action="append",
        default=[],
        metavar="SAMPLE:ORTHOLOG_ID:BP",
        help="Shift one locus marker and label horizontally by the specified bp.",
    )
    parser.add_argument(
        "--annotation-level",
        action="append",
        default=[],
        metavar="SAMPLE:ORTHOLOG_ID:LEVEL",
        help="Set the vertical label level for one locus annotation.",
    )
    parser.add_argument(
        "--collapse-matching-structures",
        action="store_true",
        help=(
            "Collapse isoforms with the same intron chain and CDS, retaining "
            "the representative with the greatest long-read support."
        ),
    )
    parser.add_argument(
        "--single-output",
        type=Path,
        help="Exact output PNG path when plotting one gene.",
    )
    parser.add_argument("--genotype-flank-length", type=int, default=100)
    parser.add_argument("--width", type=float, default=6.5)
    parser.add_argument(
        "--x-padding-bp",
        type=float,
        default=0.0,
        help="Blank genomic-coordinate padding inside each horizontal border.",
    )
    parser.add_argument(
        "--gene-introner-gap",
        type=float,
        default=1.8,
        help="Vertical distance between the gene and introner tracks.",
    )
    parser.add_argument(
        "--hide-track-separators",
        action="store_true",
        help="Omit horizontal rules between introner and isoform tracks.",
    )
    parser.add_argument(
        "--plain-introner-labels",
        action="store_true",
        help="Render strain introner-track labels at regular weight.",
    )
    parser.add_argument("--isoform-dir", type=Path, required=True)
    parser.add_argument("--sqanti-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    for label, sample in [("834", "CCMP1545"), ("1614", "RCC1614"), ("1749", "RCC1749")]:
        suffix = "scaffold_corrected.sorted.gtf" if label == "834" else "sorted.gtf"
        base.STRAIN_TO_R2C2_GTF[sample] = str(args.isoform_dir / f"09092025_{label}_Isoforms.filtered.clean.{suffix}")
        base.STRAIN_TO_SQANTI[sample] = str(args.sqanti_dir / f"sqanti3_output_{label}" / f"{label}_isoforms_classification.filtered.txt")
    if args.x_padding_bp < 0:
        parser.error("--x-padding-bp must be nonnegative")
    if args.gene_introner_gap <= 0:
        parser.error("--gene-introner-gap must be positive")
    if args.single_output and len(args.gene) != 1:
        parser.error("--single-output requires exactly one --gene")
    args.annotation_x_offsets = {}
    for specification in args.annotation_x_offset:
        try:
            sample, ortholog_id, offset = specification.split(":", 2)
            args.annotation_x_offsets[f"{sample}:{ortholog_id}"] = float(offset)
        except ValueError:
            parser.error(
                "--annotation-x-offset must use SAMPLE:ORTHOLOG_ID:BP"
            )
    args.annotation_levels = {}
    for specification in args.annotation_level:
        try:
            sample, ortholog_id, level = specification.split(":", 2)
            args.annotation_levels[f"{sample}:{ortholog_id}"] = int(level)
        except ValueError:
            parser.error(
                "--annotation-level must use SAMPLE:ORTHOLOG_ID:LEVEL"
            )
    return args


def normalize_gene(value: object) -> str:
    return base.normalize_gene_id(value)


def prepare_nmd_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table["norm_gene"] = table["associated_gene"].map(normalize_gene)
    table["is_nmd"] = table["predicted_NMD"].astype(str).str.lower().eq("true")
    table["is_high_conf_coding"] = (
        table["structural_category"].isin(HIGH_CONF_CATEGORIES)
        & table["coding"].eq("coding")
    )
    return table


def isoform_support_lookup(sample: str, gene_id: str) -> dict[str, float]:
    table = base.read_sqanti_table(
        Path(base.STRAIN_TO_SQANTI[sample]),
        gene_id,
    )
    if table.empty:
        return {}
    count_columns = [
        column
        for column in table.columns
        if re.fullmatch(r"\d+_[A-Z]\d+", column)
    ]
    if not count_columns:
        return {}
    support = (
        table[count_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .sum(axis=1)
    )
    return dict(zip(table["isoform"].astype(str), support, strict=True))


def collapse_matching_isoforms(
    panel: base.SamplePanel,
    support: dict[str, float],
) -> base.SamplePanel:
    groups: dict[tuple, list[tuple[int, base.IsoformModel]]] = {}
    order: list[tuple] = []
    for index, isoform in enumerate(panel.isoforms):
        exons = sorted(isoform.exons)
        intron_chain = tuple(
            (left[1] + 1, right[0] - 1)
            for left, right in zip(exons, exons[1:])
        )
        signature = (
            isoform.strand,
            intron_chain,
            isoform.cds_start,
            isoform.cds_end,
        )
        if signature not in groups:
            groups[signature] = []
            order.append(signature)
        groups[signature].append((index, isoform))

    representatives: list[base.IsoformModel] = []
    for signature in order:
        _, representative = max(
            groups[signature],
            key=lambda item: (
                support.get(item[1].transcript_id, 0.0),
                -item[0],
            ),
        )
        representatives.append(representative)
    panel.isoforms = representatives
    return panel


def family_lookup(genotype: pd.DataFrame) -> dict[str, int]:
    present = genotype.loc[
        genotype["presence"].eq(1), ["ortholog_id", "family"]
    ].dropna(subset=["family"])
    return (
        present.drop_duplicates("ortholog_id")
        .set_index("ortholog_id")["family"]
        .astype(int)
        .to_dict()
    )


def nmd_isoform_sets(nmd_table: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
    return (
        nmd_table[
            nmd_table["is_high_conf_coding"] & nmd_table["is_nmd"]
        ]
        .groupby(["norm_gene", "strain"])["isoform"]
        .apply(lambda values: set(values.astype(str)))
        .to_dict()
    )


def build_candidate(
    gene: str,
    args: argparse.Namespace,
    locus_gene_map: pd.DataFrame,
    families: dict[str, int],
    nmd_sets: dict[tuple[str, str], set[str]],
) -> Candidate | None:
    target = normalize_gene(gene)
    mapped = locus_gene_map[
        locus_gene_map["common_gene_id"].map(normalize_gene).eq(target)
        & locus_gene_map["sample"].isin(SAMPLES)
    ].copy()
    if mapped.empty:
        print(f"Skipping {gene}: no loci in the locus-gene map")
        return None

    gene_id = str(mapped["common_gene_id"].dropna().iloc[0])
    calls = (
        mapped.drop_duplicates(["ortholog_id", "sample"])
        .pivot(index="ortholog_id", columns="sample", values="presence")
        .reindex(columns=SAMPLES)
    )
    polymorphic_ids = [
        str(ortholog_id)
        for ortholog_id, row in calls.iterrows()
        if set(row.dropna().astype(int)) == {1, 2}
        and str(ortholog_id) not in set(args.exclude_locus)
    ]
    if not polymorphic_ids:
        print(f"Skipping {gene_id}: no callable CCMP1545/RCC1749 differences")
        return None

    mapped["sample_rank"] = pd.Categorical(
        mapped["sample"], categories=SAMPLES, ordered=True
    )
    positions = (
        mapped[mapped["ortholog_id"].isin(polymorphic_ids)]
        .sort_values(["ortholog_id", "sample_rank"])
        .drop_duplicates("ortholog_id")
        .set_index("ortholog_id")["start"]
        .to_dict()
    )
    loci = [
        CandidateLocus(
            ortholog_id=ortholog_id,
            family=families.get(ortholog_id, -1),
            ccmp_call=int(calls.loc[ortholog_id, "CCMP1545"]),
            rcc1749_call=int(calls.loc[ortholog_id, "RCC1749"]),
        )
        for ortholog_id in polymorphic_ids
    ]
    loci.sort(key=lambda locus: positions.get(locus.ortholog_id, 0))

    panels: dict[str, base.SamplePanel] = {}
    try:
        for sample in SAMPLES:
            panel = base.build_sample_panel(
                sample,
                gene_id,
                args.annotation_dir,
                flank=0,
            )
            if args.collapse_matching_structures:
                panel = collapse_matching_isoforms(
                    panel,
                    isoform_support_lookup(sample, gene_id),
                )
            panels[sample] = panel
    except (base.MissingSampleGene, SystemExit) as error:
        print(f"Skipping {gene_id}: {error}")
        return None
    if any(not panels[sample].isoforms for sample in SAMPLES):
        print(f"Skipping {gene_id}: missing plot-compatible isoforms")
        return None

    return Candidate(
        gene_id=gene_id,
        loci=loci,
        panels=panels,
        nmd_isoforms={
            sample: nmd_sets.get((target, sample), set()) for sample in SAMPLES
        },
    )


def build_candidates(
    args: argparse.Namespace,
    genotype: pd.DataFrame,
    locus_gene_map: pd.DataFrame,
    nmd_table: pd.DataFrame,
) -> list[Candidate]:
    required = {
        "ortholog_id",
        "sample",
        "presence",
        "common_gene_id",
        "start",
    }
    missing = required.difference(locus_gene_map.columns)
    if missing:
        raise SystemExit(
            "Locus-gene map is missing columns: " + ", ".join(sorted(missing))
        )
    families = family_lookup(genotype)
    nmd_sets = nmd_isoform_sets(nmd_table)
    candidates = [
        build_candidate(
            gene,
            args,
            locus_gene_map,
            families,
            nmd_sets,
        )
        for gene in args.gene
    ]
    return [candidate for candidate in candidates if candidate is not None]


def panel_geometry(panel: base.SamplePanel) -> PanelGeometry:
    return PanelGeometry(
        start=min(exon.start for exon in panel.exons),
        end=max(exon.end for exon in panel.exons),
        strand=panel.transcript.strand,
    )


def genomic_to_local_interval(
    start: int,
    end: int,
    geometry: PanelGeometry,
) -> tuple[int, int]:
    """Convert a 1-based inclusive genomic interval to local half-open bases."""
    start, end = min(start, end), max(start, end)
    if geometry.strand == "+":
        local_start = start - geometry.start
        local_end = end - geometry.start + 1
    else:
        local_start = geometry.end - end
        local_end = geometry.end - start + 1
    return (
        max(0, min(local_start, geometry.length)),
        max(0, min(local_end, geometry.length)),
    )


def genomic_to_local_boundary(
    position: int,
    geometry: PanelGeometry,
) -> int:
    local = (
        position - geometry.start
        if geometry.strand == "+"
        else geometry.end - position
    )
    return max(0, min(local, geometry.length))


def draw_rectangles(
    ax,
    intervals: list[tuple[int, int]],
    y: float,
    height: float,
    color: str,
    edge: str,
    linewidth: float = 0.5,
    zorder: int = 3,
) -> None:
    for start, end in intervals:
        ax.add_patch(
            Rectangle(
                (start, y - height / 2),
                max(1, end - start),
                height,
                facecolor=color,
                edgecolor=edge,
                linewidth=linewidth,
                zorder=zorder,
            )
        )


def draw_transcript_model(
    ax,
    exon_intervals: list[tuple[int, int]],
    cds_interval: tuple[int, int] | None,
    y: float,
    color: str,
    geometry: PanelGeometry,
    gene_model: bool = False,
) -> None:
    local_exons = [
        genomic_to_local_interval(start, end, geometry)
        for start, end in exon_intervals
    ]
    if not local_exons:
        return
    ax.plot(
        [min(start for start, _ in local_exons), max(end for _, end in local_exons)],
        [y, y],
        color="#2A2A2A",
        linewidth=0.55,
        zorder=1,
    )

    for exon_start, exon_end in exon_intervals:
        if gene_model or cds_interval is None:
            draw_rectangles(
                ax,
                [genomic_to_local_interval(exon_start, exon_end, geometry)],
                y,
                0.46 if gene_model else 0.22,
                color if gene_model else UTR_COLOR,
                "#555555" if gene_model else UTR_EDGE,
            )
            continue

        cds_start, cds_end = sorted(cds_interval)
        coding_start = max(exon_start, cds_start)
        coding_end = min(exon_end, cds_end)
        if coding_start > coding_end:
            draw_rectangles(
                ax,
                [genomic_to_local_interval(exon_start, exon_end, geometry)],
                y,
                0.22,
                UTR_COLOR,
                UTR_EDGE,
            )
            continue
        if exon_start < coding_start:
            draw_rectangles(
                ax,
                [
                    genomic_to_local_interval(
                        exon_start, coding_start - 1, geometry
                    )
                ],
                y,
                0.22,
                UTR_COLOR,
                UTR_EDGE,
            )
        draw_rectangles(
            ax,
            [genomic_to_local_interval(coding_start, coding_end, geometry)],
            y,
            0.46,
            color,
            "#555555",
        )
        if coding_end < exon_end:
            draw_rectangles(
                ax,
                [
                    genomic_to_local_interval(
                        coding_end + 1, exon_end, geometry
                    )
                ],
                y,
                0.22,
                UTR_COLOR,
                UTR_EDGE,
            )


def isoform_genomic_intervals(
    isoform: base.IsoformModel,
    panel: base.SamplePanel,
) -> tuple[list[tuple[int, int]], tuple[int, int] | None]:
    exons = [
        (start + panel.origin, end + panel.origin) for start, end in isoform.exons
    ]
    cds = None
    if isoform.cds_start is not None and isoform.cds_end is not None:
        cds = (
            isoform.cds_start + panel.origin,
            isoform.cds_end + panel.origin,
        )
    return exons, cds


def introner_x_and_interval(
    interval: base.LocusInterval,
    geometry: PanelGeometry,
) -> tuple[float, tuple[int, int] | None]:
    if interval.presence == 1 and interval.end > interval.start:
        local_interval = genomic_to_local_interval(
            interval.start, interval.end, geometry
        )
        return float(np.mean(local_interval)), local_interval
    midpoint = int(round((interval.start + interval.end) / 2))
    return float(genomic_to_local_boundary(midpoint, geometry)), None


def draw_candidate(
    candidate: Candidate,
    genotype: pd.DataFrame,
    args: argparse.Namespace,
    output_path: Path,
) -> dict[str, int | str]:
    hidden_locus_annotations = set(args.hide_locus_in_sample)
    geometries = {
        sample: panel_geometry(candidate.panels[sample]) for sample in SAMPLES
    }
    locus_ids = [locus.ortholog_id for locus in candidate.loci]
    sample_intervals = {
        sample: base.sample_locus_intervals(
            genotype,
            sample,
            locus_ids,
            annotations={},
            genotype_flank_length=args.genotype_flank_length,
        )
        for sample in SAMPLES
    }

    row_counts = {
        sample: 2 + len(candidate.panels[sample].isoforms) for sample in SAMPLES
    }
    figure_height = max(4.4, 1.35 + 0.35 * sum(row_counts.values()))
    fig, ax = plt.subplots(figsize=(args.width, figure_height))
    fig.subplots_adjust(left=0.23, right=0.985, top=0.94, bottom=0.09)

    group_span = {
        sample: (
            len(candidate.panels[sample].isoforms)
            + args.gene_introner_gap
            + 1.75
        )
        for sample in SAMPLES
    }
    group_base = {
        "RCC1749": 0.0,
        "CCMP1545": group_span["RCC1749"],
    }
    maximum_length = max(geometry.length for geometry in geometries.values())
    gene_y_values = {
        sample: (
            group_base[sample]
            + len(candidate.panels[sample].isoforms)
            + 0.15
            + args.gene_introner_gap
        )
        for sample in SAMPLES
    }
    ax.set_xlim(-args.x_padding_bp, maximum_length + args.x_padding_bp)
    ax.set_ylim(-0.65, max(gene_y_values.values()) + 0.95)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7.5, length=3, width=0.7)
    ax.set_xlabel("bp", fontsize=8.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)

    nmd_counts: dict[str, int] = {}
    for sample in SAMPLES:
        panel = candidate.panels[sample]
        geometry = geometries[sample]
        isoform_count = len(panel.isoforms)
        base_y = group_base[sample]
        introner_y = base_y + isoform_count + 0.15
        gene_y = gene_y_values[sample]

        ax.text(
            -0.015,
            gene_y,
            f"{sample} gene",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=7.2,
        )
        ax.text(
            -0.015,
            introner_y,
            f"{sample} Introners",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=7.2,
            fontweight="normal" if args.plain_introner_labels else "bold",
        )

        draw_transcript_model(
            ax,
            [(feature.start, feature.end) for feature in panel.exons],
            None,
            gene_y,
            GENE_COLORS[sample],
            geometry,
            gene_model=True,
        )
        if not args.hide_track_separators:
            ax.axhline(
                introner_y - 0.55,
                color="#C8C8C8",
                linewidth=0.45,
                zorder=0,
            )

        sample_nmd_count = 0
        for isoform_number, isoform in enumerate(panel.isoforms, start=1):
            y = base_y + isoform_count - isoform_number
            exons, cds = isoform_genomic_intervals(isoform, panel)
            is_nmd = isoform.transcript_id in candidate.nmd_isoforms[sample]
            sample_nmd_count += int(is_nmd)
            label = f"Isoform {isoform_number}" + (" (NMD)" if is_nmd else "")
            ax.text(
                -0.015,
                y,
                label,
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=7.0,
                color=NMD_COLOR if is_nmd else "#202020",
                fontweight="bold" if is_nmd else "normal",
            )
            draw_transcript_model(
                ax,
                exons,
                cds,
                y,
                SAMPLE_COLORS[sample],
                geometry,
            )
        nmd_counts[sample] = sample_nmd_count

        positions: list[tuple[float, CandidateLocus, base.LocusInterval]] = []
        for locus in candidate.loci:
            annotation_key = f"{sample}:{locus.ortholog_id}"
            if annotation_key in hidden_locus_annotations:
                continue
            interval = sample_intervals[sample].get(locus.ortholog_id)
            if interval is None or interval.contig != panel.contig:
                continue
            x, local_interval = introner_x_and_interval(interval, geometry)
            x += args.annotation_x_offsets.get(annotation_key, 0.0)
            positions.append((x, locus, interval))
            if interval.presence == 1 and local_interval is not None:
                ax.plot(
                    list(local_interval),
                    [introner_y, introner_y],
                    color=INTRONER_COLOR,
                    linewidth=5.0,
                    solid_capstyle="butt",
                    zorder=5,
                )
            elif interval.presence == 2:
                ax.scatter(
                    [x],
                    [introner_y],
                    marker="^",
                    s=24,
                    facecolors="white",
                    edgecolors=INTRONER_COLOR,
                    linewidths=1.0,
                    zorder=5,
                )
            else:
                ax.scatter(
                    [x],
                    [introner_y],
                    marker="x",
                    s=22,
                    color="#8B8B8B",
                    linewidths=1.0,
                    zorder=5,
                )

        positions.sort(key=lambda item: item[0])
        last_position_by_level: list[float] = []
        levels: list[int] = []
        minimum_separation = geometry.length * 0.12
        for x, _, _ in positions:
            level = 0
            while (
                level < len(last_position_by_level)
                and x - last_position_by_level[level] < minimum_separation
            ):
                level += 1
            if level == len(last_position_by_level):
                last_position_by_level.append(x)
            else:
                last_position_by_level[level] = x
            levels.append(level)

        for (x, locus, interval), level in zip(positions, levels, strict=True):
            level = args.annotation_levels.get(
                f"{sample}:{locus.ortholog_id}", level
            )
            state = {1: "Present", 2: "Absent", 3: "Missing"}.get(
                interval.presence, "Missing"
            )
            family_label = (
                f"Family {locus.family}" if locus.family >= 0 else "Family unknown"
            )
            ax.text(
                x,
                introner_y + 0.23 + 0.70 * level,
                f"{family_label}\n{state}",
                ha="center",
                va="bottom",
                fontsize=5.8,
                color=(
                    INTRONER_COLOR if interval.presence in {1, 2} else "#666666"
                ),
                linespacing=0.95,
            )

    fig.suptitle(
        candidate.gene_id,
        x=0.23,
        ha="left",
        fontsize=8.7,
        fontweight="bold",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    return {
        "coordinate_system": "strain_gene_local",
        "CCMP1545_region_length": geometries["CCMP1545"].length,
        "RCC1749_region_length": geometries["RCC1749"].length,
        "CCMP1545_nmd_plotted": nmd_counts["CCMP1545"],
        "RCC1749_nmd_plotted": nmd_counts["RCC1749"],
    }


def safe_stem(gene_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", normalize_gene(gene_id))


def candidate_summary_row(
    rank: int,
    candidate: Candidate,
    output: Path,
    plot_stats: dict[str, int | str],
) -> dict[str, object]:
    return {
        "candidate_rank": rank,
        "gene_id": candidate.gene_id,
        "n_discordant_loci": len(candidate.loci),
        "families": ",".join(
            map(str, sorted({locus.family for locus in candidate.loci}))
        ),
        "ortholog_ids": ",".join(
            locus.ortholog_id for locus in candidate.loci
        ),
        "CCMP1545_present_introners": sum(
            locus.ccmp_call == 1 for locus in candidate.loci
        ),
        "RCC1749_present_introners": sum(
            locus.rcc1749_call == 1 for locus in candidate.loci
        ),
        "CCMP1545_isoforms": len(candidate.panels["CCMP1545"].isoforms),
        "RCC1749_isoforms": len(candidate.panels["RCC1749"].isoforms),
        "CCMP1545_nmd_isoforms": sum(
            isoform.transcript_id in candidate.nmd_isoforms["CCMP1545"]
            for isoform in candidate.panels["CCMP1545"].isoforms
        ),
        "RCC1749_nmd_isoforms": sum(
            isoform.transcript_id in candidate.nmd_isoforms["RCC1749"]
            for isoform in candidate.panels["RCC1749"].isoforms
        ),
        "figure": str(output),
        **plot_stats,
    }


def main() -> None:
    args = parse_args()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
        }
    )

    print("Loading genotype, locus-gene, and NMD tables...")
    genotype = pd.read_csv(args.genotype_matrix, sep="\t", low_memory=False)
    genotype["presence"] = pd.to_numeric(
        genotype["presence"], errors="coerce"
    ).fillna(3).astype(int)
    locus_gene_map = pd.read_csv(args.locus_gene_map, sep="\t")
    nmd_table = prepare_nmd_table(args.nmd_table)

    print("Validating requested candidates...")
    candidates = build_candidates(args, genotype, locus_gene_map, nmd_table)
    if not candidates:
        raise SystemExit("No requested candidates could be plotted")
    if args.single_output and len(candidates) != 1:
        raise SystemExit("--single-output requires exactly one valid candidate")

    summary_rows: list[dict[str, object]] = []
    for rank, candidate in enumerate(candidates, start=1):
        output = args.single_output or (
            args.output_dir / f"candidate_{rank}.{safe_stem(candidate.gene_id)}.png"
        )
        print(
            f"Plotting candidate {rank}: {candidate.gene_id} "
            f"({len(candidate.loci)} polymorphic loci; "
            f"{len(candidate.panels['CCMP1545'].isoforms)}/"
            f"{len(candidate.panels['RCC1749'].isoforms)} isoforms)"
        )
        plot_stats = draw_candidate(candidate, genotype, args, output)
        summary_rows.append(
            candidate_summary_row(rank, candidate, output, plot_stats)
        )

    summary = pd.DataFrame(summary_rows)
    summary_path = (
        args.single_output.with_suffix(".summary.tsv")
        if args.single_output
        else args.output_dir / "candidate_summary.tsv"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote candidate summary: {summary_path}")
    for path in summary["figure"]:
        print(f"Wrote plot: {path}")


if __name__ == "__main__":
    main()

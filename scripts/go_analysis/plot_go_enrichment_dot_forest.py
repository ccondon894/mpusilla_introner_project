#!/usr/bin/env python3
"""Curated GO enrichment dot/forest plot for manuscript figures."""

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 13,
        "xtick.labelsize": 11.5,
        "ytick.labelsize": 11.5,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
    }
)


DEFAULT_INPUT = (
    "results/go_enrichment/results/"
    "introner_group_go_enrichment.significant.cds_length_permutation.tsv"
)
DEFAULT_OUTPUT = "results/figures/go_enrichment_dot_forest.png"
DEFAULT_FAMILY1_OUTPUT = "results/figures/go_enrichment_dot_forest.family1.png"

DEFAULT_MAIN_GROUPS = [
    "consistent_group1",
    "consistent_group2",
    "family_2",
]
DEFAULT_FAMILY1_GROUPS = ["family_1"]

GROUP_LABELS = {
    "all_introners": "All introners",
    "consistent_group1": "Pop1 fixed",
    "consistent_group2": "Pop2 fixed",
    "family_1": "Family 1",
    "family_2": "Family 2",
}

GROUP_COLORS = {
    "all_introners": "#0b81a2",
    "consistent_group1": "#326b77",
    "consistent_group2": "#d55e00",
    "family_1": "#80ae9a",
    "family_2": "#7DBFE0",
}

MAIN_THEMES = [
    (
        "Ribosomal and transcriptional machinery",
        [
            "GO:0003735",
            "GO:0005840",
            "GO:1990904",
            "GO:0006351",
            "GO:0140110",
        ],
    ),
    (
        "Vacuolar and vesicle trafficking",
        [
            "GO:0006623",
            "GO:0007034",
            "GO:0035542",
        ],
    ),
    (
        "Metabolism and catalytic activity",
        [
            "GO:0003824",
            "GO:0044281",
            "GO:0006767",
            "GO:0016705",
        ],
    ),
]

FAMILY1_THEMES = [
    (
        "Family 1 cell division signal",
        [
            "GO:0032954",
            "GO:1901891",
            "GO:0000917",
            "GO:0110020",
            "GO:0007265",
            "GO:0007096",
        ],
    ),
]

TERM_LABELS = {
    "GO:0003735": "Ribosomal structural\nconstituent",
    "GO:0044391": "Ribosomal subunit",
    "GO:0005840": "Ribosome",
    "GO:0022626": "Cytosolic ribosome",
    "GO:1990904": "Ribonucleoprotein\ncomplex",
    "GO:0006351": "DNA-templated\ntranscription",
    "GO:0140110": "Transcription regulator\nactivity",
    "GO:0006357": "RNA polymerase II\ntranscription regulation",
    "GO:0006623": "Protein targeting\nto vacuole",
    "GO:0072665": "Protein localization\nto vacuole",
    "GO:0007034": "Vacuolar transport",
    "GO:0035542": "SNARE complex\nassembly regulation",
    "GO:0042144": "Non-autophagic\nvacuole fusion",
    "GO:0031338": "Vesicle fusion\nregulation",
    "GO:0060627": "Vesicle-mediated\ntransport regulation",
    "GO:0003824": "Catalytic activity",
    "GO:0044281": "Small molecule\nmetabolism",
    "GO:0016491": "Oxidoreductase\nactivity",
    "GO:0006767": "Water-soluble vitamin\nmetabolism",
    "GO:0042364": "Water-soluble vitamin\nbiosynthesis",
    "GO:0016705": "O2-dependent\noxidoreductase activity",
    "GO:0032954": "Cytokinetic process\nregulation",
    "GO:1901891": "Cell septum assembly\nregulation",
    "GO:0000917": "Division septum\nassembly",
    "GO:0090529": "Cell septum assembly",
    "GO:0110020": "Actomyosin organization\nregulation",
    "GO:0007265": "Ras protein\nsignal transduction",
    "GO:0007096": "Exit from mitosis\nregulation",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot a curated, CDS-length-permutation-supported GO enrichment "
            "dot/forest panel."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Permutation TSV to plot.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="PNG output path.")
    parser.add_argument(
        "--plotted-data",
        default=None,
        help="Optional TSV with the exact rows plotted in the main panel. Defaults to output stem + .tsv.",
    )
    parser.add_argument(
        "--family1-output",
        default=DEFAULT_FAMILY1_OUTPUT,
        help="PNG output path for the Family 1 companion panel. Use 'none' to skip.",
    )
    parser.add_argument(
        "--family1-plotted-data",
        default=None,
        help="Optional TSV with the exact rows plotted in the Family 1 panel.",
    )
    parser.add_argument(
        "--groups",
        default=",".join(DEFAULT_MAIN_GROUPS),
        help="Comma-separated introner groups to include in the main panel.",
    )
    parser.add_argument(
        "--family1-groups",
        default=",".join(DEFAULT_FAMILY1_GROUPS),
        help="Comma-separated introner groups to include in the Family 1 panel.",
    )
    parser.add_argument(
        "--include-nonsignificant",
        action="store_true",
        help="Include curated rows even if empirical_significant_global is false.",
    )
    return parser.parse_args()


def as_bool(series):
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def wrap_label(value, width=42):
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False))


def display_term_label(go_term, term_name):
    if go_term in TERM_LABELS:
        return TERM_LABELS[go_term]
    return wrap_label(term_name, width=30)


def load_and_filter(path, groups, themes, include_nonsignificant):
    df = pd.read_csv(path, sep="\t")
    needed = {
        "introner_group",
        "GO_term",
        "term_name",
        "enrichment_type",
        "study_with_go",
        "study_total",
        "background_with_go",
        "background_total",
        "fold_enrichment",
        "empirical_fdr_bh_global",
        "empirical_significant_global",
    }
    missing = sorted(needed.difference(df.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {', '.join(missing)}")

    term_to_theme = {term: theme for theme, terms in themes for term in terms}
    term_order = [term for _, terms in themes for term in terms]
    df = df[df["introner_group"].isin(groups) & df["GO_term"].isin(term_to_theme)].copy()
    if not include_nonsignificant:
        df = df[as_bool(df["empirical_significant_global"])].copy()

    numeric_cols = [
        "study_with_go",
        "study_total",
        "background_with_go",
        "background_total",
        "fold_enrichment",
        "empirical_fdr_bh_global",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["theme"] = df["GO_term"].map(term_to_theme)
    df["group_label"] = df["introner_group"].map(GROUP_LABELS).fillna(df["introner_group"])
    df["log2_fold_enrichment"] = np.log2(df["fold_enrichment"].replace(0, np.nan))
    df["zero_observed_term"] = df["fold_enrichment"].eq(0)
    df["term_rank"] = df["GO_term"].map({term: i for i, term in enumerate(term_order)})
    df = df.sort_values(["term_rank", "introner_group"]).copy()
    return df, term_order


def build_y_positions(df, term_order, themes, group_by_theme=True):
    present_terms = [term for term in term_order if term in set(df["GO_term"])]
    y_positions = {}
    y_labels = {}
    theme_headers = {}
    theme_bounds = []
    y = 0

    if not group_by_theme:
        for term in present_terms:
            term_rows = df[df["GO_term"] == term]
            name = term_rows["term_name"].iloc[0]
            y_positions[term] = y
            y_labels[y] = display_term_label(term, name)
            y += 1
        return y_positions, y_labels, theme_headers, theme_bounds

    for theme, terms in themes:
        terms_here = [term for term in terms if term in present_terms]
        if not terms_here:
            continue
        theme_headers[theme] = y
        y += 0.72
        start = y
        for term in terms_here:
            term_rows = df[df["GO_term"] == term]
            name = term_rows["term_name"].iloc[0]
            y_positions[term] = y
            y_labels[y] = display_term_label(term, name)
            y += 1
        end = y - 1
        theme_bounds.append((start - 1.05, end + 0.5))
        y += 0.88

    return y_positions, y_labels, theme_headers, theme_bounds


def size_from_counts(values):
    values = np.asarray(values, dtype=float)
    return 30 + 13 * np.sqrt(np.clip(values, 0, None))


def plot(
    df,
    term_order,
    themes,
    groups,
    output,
    title,
    show_legends=True,
    group_by_theme=True,
):
    if df.empty:
        raise ValueError("No curated rows were available to plot after filtering.")

    y_positions, y_labels, theme_headers, theme_bounds = build_y_positions(
        df,
        term_order,
        themes,
        group_by_theme=group_by_theme,
    )
    offsets = np.linspace(-0.28, 0.28, num=max(len(groups), 1))
    offset_map = {group: offsets[i] for i, group in enumerate(groups)}
    df["y"] = df["GO_term"].map(y_positions) + df["introner_group"].map(offset_map)

    finite_effects = df["log2_fold_enrichment"].dropna()
    max_abs = float(np.nanmax(np.abs(finite_effects))) if not finite_effects.empty else 1.0
    xlim = max(1.75, np.ceil((max_abs + 0.25) * 2) / 2)
    zero_x = -xlim + 0.12
    df["plot_log2_fold_enrichment"] = df["log2_fold_enrichment"].fillna(zero_x)
    fig_height = max(5.4, 0.48 * len(y_labels) + 2.8)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))

    ax.axvspan(-xlim, 0, color="#eef3f5", zorder=0)
    ax.axvspan(0, xlim, color="#fff4ec", zorder=0)
    ax.axvline(0, color="#333333", linewidth=1.0, zorder=1)

    for low, high in theme_bounds:
        ax.axhline(low, color="#d7d7d7", linewidth=0.8, zorder=1)
        ax.axhline(high, color="#d7d7d7", linewidth=0.8, zorder=1)

    for group in groups:
        subset = df[df["introner_group"] == group]
        if subset.empty:
            continue
        for marker, marker_subset in [
            ("o", subset[~subset["zero_observed_term"]]),
            ("<", subset[subset["zero_observed_term"]]),
        ]:
            if marker_subset.empty:
                continue
            ax.scatter(
                marker_subset["plot_log2_fold_enrichment"],
                marker_subset["y"],
                s=size_from_counts(marker_subset["study_with_go"]),
                color=GROUP_COLORS.get(group, "#666666"),
                edgecolor="white",
                linewidth=0.8,
                alpha=0.92,
                marker=marker,
                label=GROUP_LABELS.get(group, group),
                zorder=3,
            )

    ax.set_xlim(-xlim, xlim)
    ax.set_yticks(list(y_labels))
    ax.set_yticklabels([y_labels[y] for y in y_labels], fontsize=11.5)
    ax.invert_yaxis()
    ax.set_xlabel("log2 fold enrichment vs GO-annotated background", fontsize=13)
    ax.set_ylabel("")
    ax.grid(axis="x", color="#cccccc", linewidth=0.6, alpha=0.6)
    ax.tick_params(axis="y", length=4, width=1.0, direction="out", color="#333333")
    ax.tick_params(axis="x", length=4, width=1.0, direction="out", color="#333333")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(1.0)

    for theme, header_y in theme_headers.items():
        ax.text(
            -xlim + 0.06 * (2 * xlim),
            header_y,
            theme,
            ha="left",
            va="center",
            fontsize=11.3,
            fontweight="bold",
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.8},
        )

    ax.set_title(title, fontsize=15, pad=18)

    if show_legends:
        group_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor=GROUP_COLORS.get(group, "#666666"),
                markeredgecolor="white",
                markersize=8,
                label=GROUP_LABELS.get(group, group),
            )
            for group in groups
            if group in set(df["introner_group"])
        ]
        size_values = [10, 50, 200]
        size_handles = [
            plt.scatter([], [], s=size_from_counts([value])[0], color="#777777", alpha=0.75)
            for value in size_values
        ]
        handles = group_handles + size_handles
        labels = [
            handle.get_label()
            for handle in group_handles
        ] + [
            f"{value} observed genes"
            for value in size_values
        ]
        if df["zero_observed_term"].any():
            zero_handle = Line2D(
                [0],
                [0],
                marker="<",
                linestyle="",
                markerfacecolor="#777777",
                markeredgecolor="white",
                markersize=8,
                label="0 observed genes",
            )
            handles.append(zero_handle)
            labels.append("triangle: 0 observed genes")
        ax.legend(
            handles=handles,
            labels=labels,
            title=None,
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
            frameon=True,
            facecolor="white",
            edgecolor="black",
            framealpha=0.88,
            scatterpoints=1,
            labelspacing=0.8,
            handletextpad=0.8,
            borderpad=0.6,
        )

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def write_plotted_data(df, path):
    cols = [
        "theme",
        "introner_group",
        "group_label",
        "GO_term",
        "term_name",
        "enrichment_type",
        "study_with_go",
        "study_total",
        "background_with_go",
        "background_total",
        "fold_enrichment",
        "log2_fold_enrichment",
        "empirical_p_directional",
        "empirical_fdr_bh_global",
        "empirical_significant_global",
        "zero_observed_term",
    ]
    cols = [col for col in cols if col in df.columns]
    path.parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(path, sep="\t", index=False)


def main():
    args = parse_args()
    groups = [group.strip() for group in args.groups.split(",") if group.strip()]
    output = Path(args.output)
    plotted_data = (
        Path(args.plotted_data)
        if args.plotted_data
        else output.with_suffix(".tsv")
    )
    df, term_order = load_and_filter(
        Path(args.input),
        groups,
        MAIN_THEMES,
        include_nonsignificant=args.include_nonsignificant,
    )
    plot(
        df,
        term_order,
        MAIN_THEMES,
        groups,
        output,
        "Representative GO terms supported by CDS-length-weighted permutation tests",
        group_by_theme=False,
    )
    write_plotted_data(df, plotted_data)
    print(f"Wrote {output}")
    print(f"Wrote {plotted_data}")

    if args.family1_output.lower() != "none":
        family1_groups = [
            group.strip()
            for group in args.family1_groups.split(",")
            if group.strip()
        ]
        family1_output = Path(args.family1_output)
        family1_plotted_data = (
            Path(args.family1_plotted_data)
            if args.family1_plotted_data
            else family1_output.with_suffix(".tsv")
        )
        family1_df, family1_term_order = load_and_filter(
            Path(args.input),
            family1_groups,
            FAMILY1_THEMES,
            include_nonsignificant=args.include_nonsignificant,
        )
        plot(
            family1_df,
            family1_term_order,
            FAMILY1_THEMES,
            family1_groups,
            family1_output,
            "Family 1 cell-division GO enrichments",
            show_legends=False,
        )
        write_plotted_data(family1_df, family1_plotted_data)
        print(f"Wrote {family1_output}")
        print(f"Wrote {family1_plotted_data}")


if __name__ == "__main__":
    main()

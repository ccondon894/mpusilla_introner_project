#!/usr/bin/env python3
"""
Render the genome-alignment circos plot with pyCirclize.

This is a first-pass replacement candidate for the Perl Circos render step. It
uses the same generated inputs as the existing plot:
  - final karyotype file
  - merged ribbon/link TSV
  - four histogram tracks

The plotting recipe follows pyCirclize's comparative-genomics pattern: create
sectors from sequence sizes, reverse the query genome sectors, and draw MUMmer
alignment links with ``circos.link``.
"""

import argparse
import colorsys
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from figure_color_guide import (
    DEFAULT_GUIDE_PATH,
    build_synteny_mirrored_chr_colors,
    darken_hex,
    load_color_guide,
)

BASE_COLORS = {
    "black": "#000000",
    "white": "#ffffff",
    "red": "#ff0000",
    "blue": "#0000ff",
    "grey": "#808080",
    "gray": "#808080",
    "dgrey": "#404040",
    "vdgrey": "#202020",
}
HEX_COLOR_RE = re.compile(r"^#?([0-9a-fA-F]{6})$")


def parse_args():
    parser = argparse.ArgumentParser(description="Render genome alignment with pyCirclize")
    parser.add_argument("--karyotype", required=True, help="Final Circos karyotype file")
    parser.add_argument("--ribbons", required=True, help="Merged ribbon/link TSV")
    parser.add_argument("--strain1-name", required=True, help="First strain name, e.g. CCMP1545")
    parser.add_argument("--strain2-name", required=True, help="Second strain name, e.g. RCC1749")
    parser.add_argument("--strain1-introns", required=True, help="Strain1 non-introner intron histogram")
    parser.add_argument("--strain2-introns", required=True, help="Strain2 non-introner intron histogram")
    parser.add_argument("--strain1-introners", required=True, help="Strain1 introner histogram")
    parser.add_argument("--strain2-introners", required=True, help="Strain2 introner histogram")
    parser.add_argument("--output-png", required=True, help="Output PNG path")
    parser.add_argument("--output-svg", default="", help="Optional output SVG path")
    parser.add_argument("--links", required=True, help="Synteny links TSV used for color mirroring")
    parser.add_argument(
        "--color-guide",
        default=str(DEFAULT_GUIDE_PATH),
        help="Master figure color guide TSV",
    )
    parser.add_argument("--start", type=float, default=-358, help="Start angle")
    parser.add_argument("--end", type=float, default=2, help="End angle")
    parser.add_argument("--space", type=float, default=1.5, help="Space between sectors")
    parser.add_argument("--figsize", type=float, default=10.0, help="Figure width/height in inches")
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    parser.add_argument("--link-alpha", type=float, default=0.35, help="Ribbon/link alpha")
    parser.add_argument("--link-gradient-alpha", type=float, default=None, help="Alpha for gradient links")
    parser.add_argument("--link-saturation", type=float, default=1.0, help="Multiplier for link color saturation")
    parser.add_argument("--link-radius", type=float, default=78.0, help="Link radius")
    parser.add_argument(
        "--link-gradient-steps",
        type=int,
        default=0,
        help="If >1, split links into this many color-interpolated segments",
    )
    parser.add_argument("--intron-max", type=float, default=60.0, help="Non-introner intron histogram max")
    parser.add_argument("--introner-max", type=float, default=50.0, help="Introner histogram max")
    parser.add_argument("--left-label", default="", help="Left strain label")
    parser.add_argument("--right-label", default="", help="Right strain label")
    parser.add_argument("--chromosome-label-size", type=float, default=7.5)
    parser.add_argument("--strain-label-size", type=float, default=12.0)
    parser.add_argument("--strain-label-y", type=float, default=0.95, help="Figure y-position for strain labels")
    parser.add_argument("--tick-interval", type=int, default=500000)
    parser.add_argument(
        "--exclude-sectors",
        default="",
        help="Comma/semicolon-separated sector IDs or labels to omit from the plot",
    )
    parser.add_argument(
        "--reverse-rcc-link-pairs",
        default="",
        help=(
            "Comma/semicolon-separated RCC1749:CCMP1545 chromosome pairs whose "
            "RCC1749-side link coordinates should be reversed, e.g. Chr41:Chr4"
        ),
    )
    parser.add_argument("--no-reverse-strain2", action="store_true", help="Do not reverse strain2 sectors")
    return parser.parse_args()


def import_pycirclize():
    try:
        from pycirclize import Circos
    except ImportError as exc:
        raise SystemExit(
            "pycirclize is required for this renderer. Install it in the active "
            "environment, e.g. with conda-forge::pycirclize."
        ) from exc
    return Circos


def chromosome_key(value):
    match = re.search(r"(\d+)$", str(value).strip())
    if match:
        return str(int(match.group(1)))
    return str(value).strip().lower()


def parse_color(value, alpha=None):
    value = str(value).strip()
    if not value:
        return "#808080"

    if value.endswith("_a3"):
        return parse_color(value[:-3], alpha=0.35)
    if value.endswith("_a2"):
        return parse_color(value[:-3], alpha=0.20)

    if value in BASE_COLORS:
        value = BASE_COLORS[value]

    match = HEX_COLOR_RE.match(value)
    if match:
        hex_value = match.group(1)
        rgb = tuple(int(hex_value[i:i + 2], 16) / 255 for i in (0, 2, 4))
        return (*rgb, alpha) if alpha is not None else f"#{hex_value}"

    hue_match = re.match(r"hue(\d{3})$", value)
    if hue_match:
        hue = int(hue_match.group(1)) / 360
        rgb = colorsys.hsv_to_rgb(hue, 1, 1)
        return (*rgb, alpha) if alpha is not None else rgb

    return value


def strain_for_chr_id(chr_id, strain1_name, strain2_name):
    if chr_id.startswith(strain1_name):
        return strain1_name
    if chr_id.startswith(strain2_name):
        return strain2_name
    return None


def read_karyotype_entries(path, strain1_name, strain2_name):
    entries = []
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 7 or parts[0] != "chr" or parts[1] != "-":
                continue

            chr_id = parts[2]
            entries.append(
                {
                    "chr_id": chr_id,
                    "label": parts[3],
                    "size": int(float(parts[5])) - int(float(parts[4])),
                    "strain": strain_for_chr_id(chr_id, strain1_name, strain2_name),
                    "fallback_color": parts[6],
                }
            )
    return entries


def build_karyotype_from_entries(entries, chr_colors):
    sectors = OrderedDict()
    labels = {}
    colors = {}
    strains = {}

    for entry in entries:
        chr_id = entry["chr_id"]
        sectors[chr_id] = entry["size"]
        labels[chr_id] = entry["label"]
        strains[chr_id] = entry["strain"]
        colors[chr_id] = parse_color(chr_colors.get(chr_id, entry["fallback_color"]))

    if not sectors:
        sys.exit("No chromosome entries found in karyotype")

    return sectors, labels, colors, strains


def load_links(path):
    links = []
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 6:
                continue

            chr1, start1, end1, chr2, start2, end2 = parts[:6]
            options = parts[6:]
            color = "grey"
            for option in options:
                if option.startswith("color="):
                    color = option.split("=", 1)[1]
                    break

            links.append(
                (
                    (chr1, int(start1), int(end1)),
                    (chr2, int(start2), int(end2)),
                    parse_color(color),
                )
            )
    return links


def load_histogram(path):
    by_chr = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            chrom, start, end, value = stripped.split()[:4]
            start_i = int(start)
            end_i = int(end)
            by_chr[chrom].append((start_i, end_i, float(value)))
    return by_chr


def normalize_alias(value):
    return re.sub(r"[^a-z0-9_]+", "", str(value).strip().lower())


def sector_aliases(chr_id, label):
    aliases = {normalize_alias(chr_id), normalize_alias(label)}
    key = chromosome_key(label)
    aliases.add(normalize_alias(key))
    aliases.add(normalize_alias(f"chr{key}"))
    return aliases


def parse_exclude_sectors(value):
    if not value:
        return set()
    return {
        normalize_alias(item)
        for item in re.split(r"[,;]", value)
        if item.strip()
    }


def filter_excluded_sectors(sectors, labels, colors, strains, exclude_tokens):
    if not exclude_tokens:
        return sectors, labels, colors, strains

    keep_ids = [
        chr_id
        for chr_id in sectors
        if sector_aliases(chr_id, labels.get(chr_id, chr_id)).isdisjoint(exclude_tokens)
    ]
    return (
        OrderedDict((chr_id, sectors[chr_id]) for chr_id in keep_ids),
        {chr_id: labels[chr_id] for chr_id in keep_ids},
        {chr_id: colors[chr_id] for chr_id in keep_ids},
        {chr_id: strains[chr_id] for chr_id in keep_ids},
    )


def parse_reverse_rcc_link_pairs(value):
    pairs = []
    if not value:
        return pairs

    for item in re.split(r"[,;]", value):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            rcc_token, ccmp_token = item.split(":", 1)
        else:
            fields = re.split(r"\s+to\s+", item, maxsplit=1, flags=re.IGNORECASE)
            if len(fields) != 2:
                sys.exit(f"Invalid reverse RCC link pair '{item}'. Expected RCC1749:CCMP1545")
            rcc_token, ccmp_token = fields

        pairs.append((normalize_alias(rcc_token), normalize_alias(ccmp_token)))

    return pairs


def should_flip_rcc_region(region1, region2, labels, strains, reverse_rcc_pairs, strain2_name):
    if not reverse_rcc_pairs:
        return False

    chr1, chr2 = region1[0], region2[0]
    if strains.get(chr1) == strain2_name:
        rcc_chr = chr1
        ccmp_chr = chr2
    elif strains.get(chr2) == strain2_name:
        rcc_chr = chr2
        ccmp_chr = chr1
    else:
        return False

    rcc_aliases = sector_aliases(rcc_chr, labels.get(rcc_chr, rcc_chr))
    ccmp_aliases = sector_aliases(ccmp_chr, labels.get(ccmp_chr, ccmp_chr))

    return any(
        rcc_token in rcc_aliases and ccmp_token in ccmp_aliases
        for rcc_token, ccmp_token in reverse_rcc_pairs
    )


def flip_region(region, sector_sizes):
    chrom, start, end = region
    size = sector_sizes[chrom]
    return chrom, max(0, size - end), max(0, size - start)


def maybe_flip_rcc_region(region1, region2, sector_sizes, labels, strains, reverse_rcc_pairs, strain2_name):
    if not should_flip_rcc_region(region1, region2, labels, strains, reverse_rcc_pairs, strain2_name):
        return region1, region2

    if strains.get(region1[0]) == strain2_name:
        return flip_region(region1, sector_sizes), region2
    if strains.get(region2[0]) == strain2_name:
        return region1, flip_region(region2, sector_sizes)
    return region1, region2


def plot_histogram_track(circos, sectors, histogram, r_lim, vmax, color, edgecolor):
    for sector in circos.sectors:
        track = sector.add_track(r_lim, r_pad_ratio=0.04)
        track.axis(fc="none", ec="#555555", lw=0.45)
        rows = histogram.get(sector.name, [])
        if not rows:
            continue

        x = [(start + end) / 2 for start, end, _ in rows]
        heights = [min(value, vmax) for _, _, value in rows]
        widths = [max(1, end - start + 1) for start, end, _ in rows]
        width = Counter(widths).most_common(1)[0][0]
        track.bar(
            x=x,
            height=heights,
            width=width,
            vmin=0,
            vmax=vmax,
            color=color,
            ec=edgecolor,
            lw=0.28,
            alpha=0.55,
        )


def add_ideograms(circos, labels, colors, tick_interval, label_size):
    for sector in circos.sectors:
        track = sector.add_track((94, 100), r_pad_ratio=0.02)
        track.axis(fc=colors.get(sector.name, "#cccccc"), ec="#333333", lw=0.35)
        sector.text(labels.get(sector.name, sector.name), r=103, size=label_size)
        if sector.size >= tick_interval:
            track.xticks_by_interval(
                tick_interval,
                label_formatter=lambda v: f"{v / 1_000_000:g} Mb",
                label_orientation="vertical",
                show_label=False,
                tick_length=1,
            )


def add_links(circos, links, valid_sectors, sector_sizes, labels, strains, reverse_rcc_pairs, strain2_name, radius, alpha, link_color):
    for region1, region2, color in links:
        if region1[0] not in valid_sectors or region2[0] not in valid_sectors:
            continue
        region1, region2 = maybe_flip_rcc_region(
            region1,
            region2,
            sector_sizes,
            labels,
            strains,
            reverse_rcc_pairs,
            strain2_name,
        )
        circos.link(region1, region2, color=link_color, alpha=alpha, r1=radius, r2=radius, lw=0)


def saturate_color(color, multiplier):
    if multiplier == 1:
        return color

    from matplotlib.colors import to_rgba

    rgba = to_rgba(color)
    hue, saturation, value = colorsys.rgb_to_hsv(*rgba[:3])
    saturation = max(0, min(1, saturation * multiplier))
    rgb = colorsys.hsv_to_rgb(hue, saturation, value)
    return (*rgb, rgba[3])


def interpolate_color(color1, color2, fraction, saturation_multiplier):
    from matplotlib.colors import to_rgba

    rgba1 = to_rgba(color1)
    rgba2 = to_rgba(color2)
    color = tuple(
        rgba1[i] + (rgba2[i] - rgba1[i]) * fraction
        for i in range(4)
    )
    return saturate_color(color, saturation_multiplier)


def split_region(region, fraction_start, fraction_end):
    chrom, start, end = region
    new_start = start + (end - start) * fraction_start
    new_end = start + (end - start) * fraction_end
    return chrom, new_start, new_end


def add_gradient_links(
    circos,
    links,
    valid_sectors,
    sector_sizes,
    sector_colors,
    labels,
    strains,
    reverse_rcc_pairs,
    strain2_name,
    radius,
    alpha,
    saturation_multiplier,
    steps,
):
    steps = max(1, steps)
    for region1, region2, fallback_color in links:
        if region1[0] not in valid_sectors or region2[0] not in valid_sectors:
            continue
        region1, region2 = maybe_flip_rcc_region(
            region1,
            region2,
            sector_sizes,
            labels,
            strains,
            reverse_rcc_pairs,
            strain2_name,
        )

        color1 = sector_colors.get(region1[0], fallback_color)
        color2 = sector_colors.get(region2[0], fallback_color)

        for i in range(steps):
            fraction_start = i / steps
            fraction_end = (i + 1) / steps
            fraction_mid = (fraction_start + fraction_end) / 2
            color = interpolate_color(color1, color2, fraction_mid, saturation_multiplier)
            circos.link(
                split_region(region1, fraction_start, fraction_end),
                split_region(region2, fraction_start, fraction_end),
                color=color,
                alpha=alpha,
                r1=radius,
                r2=radius,
                lw=0,
            )


def save_figure(fig, output_png, output_svg, dpi):
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    if output_svg:
        Path(output_svg).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_svg, bbox_inches="tight")


def main():
    args = parse_args()
    mpl_config_dir = Path(args.output_png).parent / ".matplotlib"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    Circos = import_pycirclize()
    guide_colors = load_color_guide(args.color_guide)
    entries = read_karyotype_entries(args.karyotype, args.strain1_name, args.strain2_name)
    labels = {entry["chr_id"]: entry["label"] for entry in entries}
    strains = {entry["chr_id"]: entry["strain"] for entry in entries}
    mirrored_chr_colors = build_synteny_mirrored_chr_colors(
        args.links,
        [entry["chr_id"] for entry in entries],
        labels,
        strains,
        args.strain1_name,
        args.strain2_name,
        guide=guide_colors,
    )
    sectors, labels, colors, strains = build_karyotype_from_entries(entries, mirrored_chr_colors)
    sectors, labels, colors, strains = filter_excluded_sectors(
        sectors,
        labels,
        colors,
        strains,
        parse_exclude_sectors(args.exclude_sectors),
    )
    sector2clockwise = {}
    if not args.no_reverse_strain2:
        sector2clockwise = {
            sector: False
            for sector, strain in strains.items()
            if strain == args.strain2_name
        }

    circos = Circos(
        sectors=sectors,
        start=args.start,
        end=args.end,
        space=args.space,
        sector2clockwise=sector2clockwise,
    )

    add_ideograms(circos, labels, colors, args.tick_interval, args.chromosome_label_size)
    intron_color = guide_colors["Intron"]
    introner_color = guide_colors["Introner"]
    link_color = guide_colors["Ideogram link"]
    plot_histogram_track(
        circos,
        sectors,
        load_histogram(args.strain1_introns) | load_histogram(args.strain2_introns),
        (70, 82),
        args.intron_max,
        intron_color,
        darken_hex(intron_color),
    )
    plot_histogram_track(
        circos,
        sectors,
        load_histogram(args.strain1_introners) | load_histogram(args.strain2_introners),
        (82, 94),
        args.introner_max,
        introner_color,
        darken_hex(introner_color),
    )
    links = load_links(args.ribbons)
    reverse_rcc_pairs = parse_reverse_rcc_link_pairs(args.reverse_rcc_link_pairs)
    link_alpha = args.link_gradient_alpha if args.link_gradient_alpha is not None else args.link_alpha
    if args.link_gradient_steps > 1:
        add_gradient_links(
            circos,
            links,
            set(sectors),
            sectors,
            colors,
            labels,
            strains,
            reverse_rcc_pairs,
            args.strain2_name,
            args.link_radius,
            link_alpha,
            args.link_saturation,
            args.link_gradient_steps,
        )
    else:
        add_links(
            circos,
            links,
            set(sectors),
            sectors,
            labels,
            strains,
            reverse_rcc_pairs,
            args.strain2_name,
            args.link_radius,
            args.link_alpha,
            link_color,
        )

    fig = circos.plotfig(figsize=(args.figsize, args.figsize))
    if args.left_label:
        fig.text(0.075, args.strain_label_y, args.left_label, ha="left", va="center", fontsize=args.strain_label_size, weight="bold")
    if args.right_label:
        fig.text(0.925, args.strain_label_y, args.right_label, ha="right", va="center", fontsize=args.strain_label_size, weight="bold")

    save_figure(fig, args.output_png, args.output_svg, args.dpi)


if __name__ == "__main__":
    main()

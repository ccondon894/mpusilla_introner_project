#!/usr/bin/env python3
"""Load figure colors from master_figure_color_guide.tsv."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GUIDE_PATH = REPO_ROOT / "master_figure_color_guide.tsv"

POPULATION_1_DARK = "Population 1 dark"
POPULATION_1_LIGHT = "Population 1 light"
POPULATION_2_DARK = "Population 2 dark"
POPULATION_2_LIGHT = "Population 2 light"


def normalize_hex(value: str) -> str:
    value = str(value).strip()
    if not value.startswith("#"):
        value = f"#{value}"
    return value.upper()


def load_color_guide(path: Path | str | None = None) -> dict[str, str]:
    guide_path = Path(path) if path else DEFAULT_GUIDE_PATH
    colors: dict[str, str] = {}
    with open(guide_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            name = (row.get("name") or "").strip()
            hex_code = (row.get("hex_code") or "").strip()
            if name and hex_code:
                colors[name] = normalize_hex(hex_code)
    return colors


def get_color(name: str, guide: dict[str, str] | None = None, path: Path | str | None = None) -> str:
    colors = guide if guide is not None else load_color_guide(path)
    if name not in colors:
        raise KeyError(f"Color '{name}' not found in figure color guide")
    return colors[name]


def population_chromosome_color(
    strain: str | None,
    anchor_index: int,
    strain1_name: str,
    strain2_name: str,
    guide: dict[str, str],
) -> str:
    """Alternating dark/light blue (pop1) or orange (pop2) by synteny anchor index."""
    if not anchor_index or strain is None:
        return "#808080"

    shade_key = POPULATION_1_DARK if anchor_index % 2 == 1 else POPULATION_1_LIGHT
    if strain == strain2_name:
        shade_key = POPULATION_2_DARK if anchor_index % 2 == 1 else POPULATION_2_LIGHT
    elif strain != strain1_name:
        return "#808080"

    if shade_key in guide:
        return guide[shade_key]

    base = guide.get("Population 1" if strain == strain1_name else "Population 2", "#808080")
    return darken_hex(base) if anchor_index % 2 == 1 else lighten_hex(base)


def lighten_hex(hex_color: str, factor: float = 1.35) -> str:
    import colorsys

    hex_color = normalize_hex(hex_color).lstrip("#")
    rgb = tuple(int(hex_color[i : i + 2], 16) / 255 for i in (0, 2, 4))
    hue, saturation, value = colorsys.rgb_to_hsv(*rgb)
    value = max(0.0, min(1.0, value * factor))
    saturation = max(0.0, min(1.0, saturation * 0.85))
    light_rgb = colorsys.hsv_to_rgb(hue, saturation, value)
    return "#{:02X}{:02X}{:02X}".format(
        int(light_rgb[0] * 255),
        int(light_rgb[1] * 255),
        int(light_rgb[2] * 255),
    )


def parse_link_record(parts: list[str]) -> tuple[str, str, int] | None:
    if len(parts) < 6:
        return None
    if parts[0].startswith("link"):
        parts = parts[1:]
    if len(parts) < 6:
        return None

    chr1, start1, end1, chr2, start2, end2 = parts[:6]
    try:
        start1_i = int(float(start1))
        end1_i = int(float(end1))
        start2_i = int(float(start2))
        end2_i = int(float(end2))
    except ValueError:
        return None

    length = min(abs(end1_i - start1_i), abs(end2_i - start2_i))
    return chr1, chr2, length


def accumulate_synteny_stats(links_path: Path | str) -> dict[tuple[str, str], int]:
    stats: dict[tuple[str, str], int] = defaultdict(int)
    with open(links_path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            record = parse_link_record(stripped.split())
            if record is None:
                continue
            chr1, chr2, length = record
            stats[(chr1, chr2)] += length
    return stats


def best_strain2_to_strain1_partners(
    stats: dict[tuple[str, str], int],
    strain1_name: str,
    strain2_name: str,
) -> dict[str, str]:
    """Map each strain2 chromosome ID to its best strain1 synteny partner."""
    candidates: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for (chr1, chr2), length in stats.items():
        if chr1.startswith(strain1_name) and chr2.startswith(strain2_name):
            candidates[chr2].append((chr1, length))
        elif chr2.startswith(strain1_name) and chr1.startswith(strain2_name):
            candidates[chr1].append((chr2, length))

    return {
        strain2_chr: max(partners, key=lambda item: item[1])[0]
        for strain2_chr, partners in candidates.items()
    }


def anchor_index_for_chr(
    chr_id: str,
    label: str,
    strain: str | None,
    strain1_name: str,
    strain2_name: str,
    strain2_to_strain1: dict[str, str],
    labels: dict[str, str],
) -> int | None:
    if strain == strain1_name:
        return parse_chromosome_number(label)

    partner_id = strain2_to_strain1.get(chr_id)
    if partner_id:
        return parse_chromosome_number(labels.get(partner_id, ""))
    return parse_chromosome_number(label)


def build_synteny_mirrored_chr_colors(
    links_path: Path | str,
    chr_ids: list[str],
    labels: dict[str, str],
    strains: dict[str, str | None],
    strain1_name: str,
    strain2_name: str,
    guide: dict[str, str] | None = None,
    path: Path | str | None = None,
) -> dict[str, str]:
    """Assign pop1 blue / pop2 orange shades by shared synteny anchor index."""
    colors = guide if guide is not None else load_color_guide(path)
    stats = accumulate_synteny_stats(links_path)
    strain2_to_strain1 = best_strain2_to_strain1_partners(stats, strain1_name, strain2_name)
    mirrored: dict[str, str] = {}

    for chr_id in chr_ids:
        strain = strains.get(chr_id)
        anchor_index = anchor_index_for_chr(
            chr_id,
            labels.get(chr_id, ""),
            strain,
            strain1_name,
            strain2_name,
            strain2_to_strain1,
            labels,
        )
        mirrored[chr_id] = population_chromosome_color(
            strain,
            anchor_index or 0,
            strain1_name,
            strain2_name,
            colors,
        )

    return mirrored


def darken_hex(hex_color: str, factor: float = 0.75) -> str:
    import colorsys

    hex_color = normalize_hex(hex_color).lstrip("#")
    rgb = tuple(int(hex_color[i : i + 2], 16) / 255 for i in (0, 2, 4))
    hue, saturation, value = colorsys.rgb_to_hsv(*rgb)
    value = max(0.0, min(1.0, value * factor))
    saturation = max(0.0, min(1.0, saturation * 1.05))
    dark_rgb = colorsys.hsv_to_rgb(hue, saturation, value)
    return "#{:02X}{:02X}{:02X}".format(
        int(dark_rgb[0] * 255),
        int(dark_rgb[1] * 255),
        int(dark_rgb[2] * 255),
    )


def parse_chromosome_number(label: str) -> int | None:
    match = re.search(r"(\d+)$", str(label).strip())
    if not match:
        return None
    return int(match.group(1))

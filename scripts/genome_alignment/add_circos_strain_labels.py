#!/usr/bin/env python3
"""Add left/right strain labels to a Circos SVG."""

import argparse
import html
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Add side labels to a Circos SVG")
    parser.add_argument("--input-svg", required=True)
    parser.add_argument("--output-svg", required=True)
    parser.add_argument("--left-label", required=True)
    parser.add_argument("--right-label", required=True)
    parser.add_argument("--font-size", type=float, default=36)
    parser.add_argument("--x-padding", type=float, default=25)
    return parser.parse_args()


def numeric_dimension(value):
    match = re.match(r"\s*([0-9.]+)", value or "")
    return float(match.group(1)) if match else None


def svg_dimensions(svg_text):
    tag_match = re.search(r"<svg\b[^>]*>", svg_text)
    if not tag_match:
        sys.exit("Input does not look like an SVG file")

    svg_tag = tag_match.group(0)
    width_match = re.search(r'\bwidth="([^"]+)"', svg_tag)
    height_match = re.search(r'\bheight="([^"]+)"', svg_tag)
    viewbox_match = re.search(r'\bviewBox="([^"]+)"', svg_tag)

    width = numeric_dimension(width_match.group(1)) if width_match else None
    height = numeric_dimension(height_match.group(1)) if height_match else None

    if (width is None or height is None) and viewbox_match:
        fields = viewbox_match.group(1).split()
        if len(fields) == 4:
            width = width if width is not None else float(fields[2])
            height = height if height is not None else float(fields[3])

    if width is None or height is None:
        sys.exit("Could not determine SVG width and height")

    return width, height


def label_elements(width, height, left_label, right_label, font_size, x_padding):
    y = height / 2
    right_x = width - x_padding
    common = (
        f' y="{y:.3f}" dominant-baseline="middle" '
        f'font-family="DejaVu Sans, Arial, sans-serif" '
        f'font-size="{font_size:g}" font-weight="700" fill="black"'
    )
    left = (
        f'<text x="{x_padding:.3f}"{common} text-anchor="start">'
        f'{html.escape(left_label)}</text>'
    )
    right = (
        f'<text x="{right_x:.3f}"{common} text-anchor="end">'
        f'{html.escape(right_label)}</text>'
    )
    return f"\n<g id=\"strain-side-labels\">\n  {left}\n  {right}\n</g>\n"


def main():
    args = parse_args()
    with open(args.input_svg) as handle:
        svg_text = handle.read()

    width, height = svg_dimensions(svg_text)
    labels = label_elements(
        width,
        height,
        args.left_label,
        args.right_label,
        args.font_size,
        args.x_padding,
    )

    if "</svg>" not in svg_text:
        sys.exit("Input SVG is missing closing </svg> tag")

    svg_text = svg_text.replace("</svg>", labels + "</svg>", 1)
    with open(args.output_svg, "w") as handle:
        handle.write(svg_text)


if __name__ == "__main__":
    main()

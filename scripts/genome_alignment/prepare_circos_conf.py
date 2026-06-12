#!/usr/bin/env python3
"""Copy the Circos config and ensure generated palette colors are included."""

import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare Circos config for rendering")
    parser.add_argument("--input", required=True, help="Source circos.conf")
    parser.add_argument("--output", required=True, help="Prepared circos.conf")
    parser.add_argument(
        "--palette-include",
        default="circos_palette_colors.conf",
        help="Color include file name relative to the render directory",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.input) as handle:
        conf = handle.read()

    include_line = f"<<include {args.palette_include}>>"
    if include_line not in conf:
        colors_close = "</colors>"
        if colors_close not in conf:
            raise SystemExit("Could not find </colors> block in Circos config")
        conf = conf.replace(colors_close, f"    {include_line}\n\n{colors_close}", 1)

    with open(args.output, "w") as handle:
        handle.write(conf)


if __name__ == "__main__":
    main()

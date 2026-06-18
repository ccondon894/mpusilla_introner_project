#!/usr/bin/env python3
"""Write present-introner body BED records for one sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genotype_matrix", type=Path)
    parser.add_argument("sample")
    parser.add_argument("output_bed", type=Path)
    parser.add_argument("--flank-length", type=int, default=100)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.genotype_matrix, sep="\t")
    for col in ("start", "end", "presence"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    rows = df[(df["sample"] == args.sample) & (df["presence"] == 1)].copy()
    rows = rows.sort_values(["contig", "start", "end", "ortholog_id"])

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.output_bed.open("w") as out:
        for _, row in rows.iterrows():
            if pd.isna(row["contig"]) or pd.isna(row["start"]) or pd.isna(row["end"]):
                continue
            start = int(row["start"]) + args.flank_length
            end = int(row["end"]) - args.flank_length
            if end <= start:
                continue
            name = f"{row['ortholog_id']}|{args.sample}"
            out.write(
                "\t".join(
                    [
                        str(row["contig"]),
                        str(start),
                        str(end),
                        name,
                        "0",
                        "+",
                    ]
                )
                + "\n"
            )
            written += 1

    print(f"Wrote {written} introner-body BED records for {args.sample}")


if __name__ == "__main__":
    main()

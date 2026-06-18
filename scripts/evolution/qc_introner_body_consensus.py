#!/usr/bin/env python3
"""Strict QC for introner-body consensus sequences from present calls."""

from __future__ import annotations

import argparse
import bisect
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = {"M", "D", "N", "=", "X"}
QUERY_OPS = {"M", "I", "S", "=", "X"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--consensus-fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bam", type=Path)
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument("--reference-sample", default="CCMP1545")
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--min-callable-fraction", type=float, default=0.90)
    parser.add_argument("--max-n-fraction", type=float, default=0.10)
    parser.add_argument("--clip-read-fraction", type=float, default=0.20)
    parser.add_argument("--max-highclip-read-fraction", type=float, default=0.15)
    parser.add_argument("--min-major-allele-fraction", type=float, default=0.80)
    parser.add_argument("--max-low-concordance-fraction", type=float, default=0.10)
    parser.add_argument("--concordance-merge-gap", type=int, default=2000)
    parser.add_argument("--max-depth-ratio", type=float, default=3.0)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    return parser.parse_args()


def read_bed(path: Path, sample: str) -> pd.DataFrame:
    cols = ["contig", "body_start", "body_end", "name", "score", "strand"]
    bed = pd.read_csv(
        path,
        sep="\t",
        names=cols,
        dtype={"contig": str, "name": str, "strand": str},
    )
    if bed.empty:
        return bed
    bed["name"] = bed["name"].astype(str)
    bed["ortholog_id"] = bed["name"].str.split("|", regex=False).str[0]
    bed["sample"] = sample
    bed["body_len"] = bed["body_end"].astype(int) - bed["body_start"].astype(int)
    return bed


def read_consensus(path: Path) -> dict[str, str]:
    seqs = {}
    for record in SeqIO.parse(path, "fasta"):
        name = record.id.split("::", 1)[0]
        seqs[name] = str(record.seq).upper()
    return seqs


def interval_index(bed: pd.DataFrame):
    by_contig = {}
    for contig, sub in bed.groupby("contig", sort=False):
        ordered = sub.sort_values("body_start")
        intervals = list(
            zip(
                ordered["body_start"].astype(int),
                ordered["body_end"].astype(int),
                ordered["name"],
            )
        )
        starts = [item[0] for item in intervals]
        by_contig[contig] = (starts, intervals)
    return by_contig


def overlapping_names(index, contig: str, start: int, end: int) -> list[str]:
    if contig not in index:
        return []
    starts, intervals = index[contig]
    stop = bisect.bisect_left(starts, end)
    names = []
    for idx in range(stop - 1, -1, -1):
        iv_start, iv_end, name = intervals[idx]
        if iv_end <= start:
            break
        if iv_end > start and iv_start < end:
            names.append(name)
    return names


def depth_stats(args: argparse.Namespace, bed: pd.DataFrame) -> dict[str, dict[str, float]]:
    stats = {
        row["name"]: {"depth_sum": 0.0, "depth_positions": 0, "callable_positions": 0}
        for _, row in bed.iterrows()
    }
    if bed.empty or args.sample == args.reference_sample:
        for row in stats.values():
            row["mean_depth"] = math.nan
            row["callable_fraction"] = 1.0
        return stats

    cmd = [
        "samtools",
        "depth",
        "-aa",
        "-q",
        str(args.min_baseq),
        "-Q",
        str(args.min_mapq),
        "-b",
        str(args.bed),
        str(args.bam),
    ]
    index = interval_index(bed)
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            if not line:
                continue
            contig, pos, depth = line.split("\t")[:3]
            pos0 = int(pos) - 1
            depth_value = int(depth)
            for name in overlapping_names(index, contig, pos0, pos0 + 1):
                stats[name]["depth_positions"] += 1
                stats[name]["depth_sum"] += depth_value
                if depth_value >= args.min_depth:
                    stats[name]["callable_positions"] += 1
        if proc.wait() != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    body_lens = dict(zip(bed["name"], bed["body_len"].astype(int)))
    for name, row in stats.items():
        body_len = body_lens[name]
        row["mean_depth"] = row["depth_sum"] / body_len if body_len > 0 else math.nan
        row["callable_fraction"] = (
            row["callable_positions"] / body_len if body_len > 0 else 0.0
        )
    return stats


def cigar_lengths(cigar: str) -> tuple[int, int, int]:
    ref_len = 0
    query_len = 0
    soft = 0
    for length_s, op in CIGAR_RE.findall(cigar):
        length = int(length_s)
        if op in REF_OPS:
            ref_len += length
        if op in QUERY_OPS:
            query_len += length
        if op == "S":
            soft += length
    return ref_len, query_len, soft


def clip_stats(args: argparse.Namespace, bed: pd.DataFrame) -> dict[str, dict[str, float]]:
    stats = {
        row["name"]: {"overlapping_reads": 0, "highclip_reads": 0}
        for _, row in bed.iterrows()
    }
    if bed.empty or args.sample == args.reference_sample:
        for row in stats.values():
            row["softclip_read_fraction"] = 0.0
        return stats

    cmd = [
        "samtools",
        "view",
        "-F",
        "3340",
        "-q",
        str(args.min_mapq),
        "-L",
        str(args.bed),
        str(args.bam),
    ]
    index = interval_index(bed)
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            contig = fields[2]
            if contig == "*":
                continue
            start = int(fields[3]) - 1
            cigar = fields[5]
            if cigar == "*":
                continue
            ref_len, query_len, soft = cigar_lengths(cigar)
            if ref_len <= 0 or query_len <= 0:
                continue
            end = start + ref_len
            soft_fraction = soft / query_len
            for name in overlapping_names(index, contig, start, end):
                stats[name]["overlapping_reads"] += 1
                if soft_fraction > args.clip_read_fraction:
                    stats[name]["highclip_reads"] += 1
        if proc.wait() != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    for row in stats.values():
        total = row["overlapping_reads"]
        row["softclip_read_fraction"] = row["highclip_reads"] / total if total else 0.0
    return stats


def parse_pileup_bases(ref_base: str, bases: str) -> dict[str, int]:
    counts = {base: 0 for base in "ACGT"}
    i = 0
    ref_base = ref_base.upper()
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            length_start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            if length_start == i:
                continue
            indel_len = int(bases[length_start:i])
            i += indel_len
            continue
        if char in ".,":
            if ref_base in counts:
                counts[ref_base] += 1
        else:
            base = char.upper()
            if base in counts:
                counts[base] += 1
        i += 1
    return counts


def concordance_stats(args: argparse.Namespace, bed: pd.DataFrame) -> dict[str, dict[str, float]]:
    stats = {
        row["name"]: {
            "concordance_positions": 0,
            "low_concordance_positions": 0,
            "major_allele_fraction_sum": 0.0,
        }
        for _, row in bed.iterrows()
    }
    if bed.empty or args.sample == args.reference_sample:
        for _, bed_row in bed.iterrows():
            stats[bed_row["name"]].update(
                {
                    "mean_major_allele_fraction": 1.0,
                    "low_concordance_fraction": 0.0,
                    "concordance_positions": int(bed_row["body_len"]),
                }
            )
        return stats

    if args.reference_fasta is None:
        raise ValueError("--reference-fasta is required for non-reference samples")

    index = interval_index(bed)
    for contig, sub in bed.groupby("contig", sort=False):
        ordered = sub.sort_values("body_start")
        chunks = []
        chunk_start = None
        chunk_end = None
        for start, end in ordered[["body_start", "body_end"]].astype(int).itertuples(index=False):
            if chunk_start is None or start > chunk_end + args.concordance_merge_gap:
                if chunk_start is not None:
                    chunks.append((chunk_start, chunk_end))
                chunk_start, chunk_end = start, end
            else:
                chunk_end = max(chunk_end, end)
        if chunk_start is not None:
            chunks.append((chunk_start, chunk_end))

        for start, end in chunks:
            region = f"{contig}:{start + 1}-{end}"
            cmd = [
                "samtools",
                "mpileup",
                "-aa",
                "-q",
                str(args.min_mapq),
                "-Q",
                str(args.min_baseq),
                "-f",
                str(args.reference_fasta),
                "-r",
                region,
                str(args.bam),
            ]
            with subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
            ) as proc:
                assert proc.stdout is not None
                for line in proc.stdout:
                    if not line:
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        continue
                    pileup_contig, pos, ref_base, depth_s, bases = fields[:5]
                    pos0 = int(pos) - 1
                    depth = int(depth_s)
                    names = overlapping_names(index, pileup_contig, pos0, pos0 + 1)
                    if not names or depth < args.min_depth:
                        continue
                    counts = parse_pileup_bases(ref_base, bases)
                    allele_depth = sum(counts.values())
                    if allele_depth < args.min_depth:
                        continue
                    major_fraction = max(counts.values()) / allele_depth
                    for name in names:
                        stats[name]["concordance_positions"] += 1
                        stats[name]["major_allele_fraction_sum"] += major_fraction
                        if major_fraction < args.min_major_allele_fraction:
                            stats[name]["low_concordance_positions"] += 1
                if proc.wait() != 0:
                    raise subprocess.CalledProcessError(proc.returncode, cmd)

    body_lens = dict(zip(bed["name"], bed["body_len"].astype(int)))
    for name, row in stats.items():
        body_len = body_lens[name]
        n_positions = row["concordance_positions"]
        row["mean_major_allele_fraction"] = (
            row["major_allele_fraction_sum"] / n_positions if n_positions else math.nan
        )
        row["low_concordance_fraction"] = (
            row["low_concordance_positions"] / body_len if body_len > 0 else 1.0
        )
    return stats


def fail_reasons(row: dict[str, object]) -> str:
    reasons = []
    if not row["length_ok"]:
        reasons.append("length_mismatch")
    if row["body_len"] <= 0:
        reasons.append("nonpositive_body_len")
    if row["n_fraction"] > row["max_n_fraction"]:
        reasons.append("high_n_fraction")
    if row["callable_fraction"] < row["min_callable_fraction"]:
        reasons.append("low_callable_fraction")
    if row["low_concordance_fraction"] > row["max_low_concordance_fraction"]:
        reasons.append("low_consensus_concordance")
    if row["depth_ratio"] > row["max_depth_ratio"]:
        reasons.append("high_depth_ratio")
    if row["consensus_missing"]:
        reasons.append("consensus_missing")
    return ",".join(reasons) if reasons else "PASS"


def main() -> None:
    args = parse_args()
    if args.sample != args.reference_sample and args.bam is None:
        raise ValueError("--bam is required for non-reference samples")

    bed = read_bed(args.bed, args.sample)
    consensus = read_consensus(args.consensus_fasta)
    depth = depth_stats(args, bed)
    clips = clip_stats(args, bed)
    concordance = concordance_stats(args, bed)

    mean_depths = [
        values["mean_depth"]
        for values in depth.values()
        if not math.isnan(values["mean_depth"]) and values["mean_depth"] > 0
    ]
    median_depth = float(np.median(mean_depths)) if mean_depths else math.nan

    rows = []
    for _, bed_row in bed.iterrows():
        name = bed_row["name"]
        seq = consensus.get(name, "")
        body_len = int(bed_row["body_len"])
        n_count = seq.count("N")
        n_fraction = n_count / len(seq) if seq else 1.0
        mean_depth = depth[name]["mean_depth"]
        depth_ratio = (
            mean_depth / median_depth
            if median_depth and not math.isnan(median_depth) and median_depth > 0
            else 1.0
        )
        row = {
            "sample": args.sample,
            "ortholog_id": bed_row["ortholog_id"],
            "contig": bed_row["contig"],
            "body_start": int(bed_row["body_start"]),
            "body_end": int(bed_row["body_end"]),
            "body_len": body_len,
            "consensus_source": (
                "reference_assembly"
                if args.sample == args.reference_sample
                else "read_consensus"
            ),
            "consensus_len": len(seq),
            "consensus_missing": name not in consensus,
            "length_ok": len(seq) == body_len and body_len > 0,
            "n_fraction": n_fraction,
            "callable_fraction": depth[name]["callable_fraction"],
            "mean_depth": mean_depth,
            "sample_median_body_depth": median_depth,
            "depth_ratio": depth_ratio,
            "overlapping_reads": clips[name]["overlapping_reads"],
            "highclip_reads": clips[name]["highclip_reads"],
            "softclip_read_fraction": clips[name]["softclip_read_fraction"],
            "concordance_positions": concordance[name]["concordance_positions"],
            "low_concordance_positions": concordance[name]["low_concordance_positions"],
            "mean_major_allele_fraction": concordance[name]["mean_major_allele_fraction"],
            "low_concordance_fraction": concordance[name]["low_concordance_fraction"],
            "min_callable_fraction": args.min_callable_fraction,
            "max_n_fraction": args.max_n_fraction,
            "max_highclip_read_fraction": args.max_highclip_read_fraction,
            "min_major_allele_fraction": args.min_major_allele_fraction,
            "max_low_concordance_fraction": args.max_low_concordance_fraction,
            "max_depth_ratio": args.max_depth_ratio,
        }
        row["qc_pass"] = fail_reasons(row) == "PASS"
        row["qc_fail_reason"] = fail_reasons(row)
        rows.append(row)

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out)} QC rows for {args.sample}: {int(out['qc_pass'].sum())} pass")


if __name__ == "__main__":
    main()

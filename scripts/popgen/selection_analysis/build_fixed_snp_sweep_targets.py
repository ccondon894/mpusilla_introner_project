#!/usr/bin/env python3
"""Build 4D SNP targets fixed between Group 1 and Group 2."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="All-sample 4D VCF.")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--group1-samples", required=True)
    parser.add_argument("--group2-samples", required=True)
    parser.add_argument(
        "--min-spacing",
        type=int,
        default=0,
        help="Minimum bp spacing between retained fixed SNP targets on the same contig. Use 0 to disable thinning.",
    )
    parser.add_argument("--mating-contig", required=True)
    parser.add_argument("--mating-start", type=int, required=True)
    parser.add_argument("--mating-end", type=int, required=True)
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def open_text(path: str | Path):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_field: str) -> int | None:
    gt = sample_field if len(sample_field) == 1 else sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|.", ""}:
        return None
    alt_count = 0
    for token in gt.replace("|", "/").split("/"):
        if token in {".", ""}:
            return None
        if token != "0":
            alt_count += 1
    return 1 if alt_count > 0 else 0


def overlaps_interval(
    contig: str, start: int, end: int, query_contig: str, query_start: int, query_end: int
) -> bool:
    return contig == query_contig and start <= query_end and end >= query_start


def allele_label(ref: str, alt: str, genotype: int) -> str:
    return alt if genotype else ref


def fixed_category(group1_gt: int, group2_gt: int) -> str:
    if group1_gt == 0 and group2_gt == 1:
        return "fixed_snp_group1_ref_group2_alt"
    if group1_gt == 1 and group2_gt == 0:
        return "fixed_snp_group1_alt_group2_ref"
    raise ValueError("Fixed SNP category requires different fixed alleles")


def build_fixed_snp_targets(
    vcf: str | Path,
    group1_samples: list[str],
    group2_samples: list[str],
    min_spacing: int,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    rows: list[dict[str, object]] = []
    stats = {
        "records_seen": 0,
        "skipped_missing": 0,
        "skipped_multiallelic": 0,
        "skipped_non_snp": 0,
        "skipped_not_fixed_within_groups": 0,
        "skipped_same_fixed_allele": 0,
        "skipped_mating_type": 0,
        "fixed_snp_candidates_before_spacing": 0,
        "skipped_min_spacing": 0,
        "targets_retained": 0,
    }
    sample_cols: dict[str, int] | None = None

    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [sample for sample in group1_samples + group2_samples if sample not in vcf_samples]
                if missing:
                    raise ValueError(f"VCF missing requested samples: {missing}")
                sample_cols = {sample: 9 + vcf_samples.index(sample) for sample in group1_samples + group2_samples}
                continue
            if not line or line.startswith("#"):
                continue
            if sample_cols is None:
                raise ValueError("VCF sample header was not found before data records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sample_cols.values()):
                stats["skipped_missing"] += 1
                continue
            stats["records_seen"] += 1

            contig = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            if "," in alt:
                stats["skipped_multiallelic"] += 1
                continue
            if alt in {".", ""} or len(ref) != 1 or len(alt) != 1:
                stats["skipped_non_snp"] += 1
                continue

            g1 = [parse_gt(fields[sample_cols[sample]]) for sample in group1_samples]
            g2 = [parse_gt(fields[sample_cols[sample]]) for sample in group2_samples]
            if any(gt is None for gt in g1 + g2):
                stats["skipped_missing"] += 1
                continue
            if len(set(g1)) != 1 or len(set(g2)) != 1:
                stats["skipped_not_fixed_within_groups"] += 1
                continue
            group1_gt = int(g1[0])
            group2_gt = int(g2[0])
            if group1_gt == group2_gt:
                stats["skipped_same_fixed_allele"] += 1
                continue
            if overlaps_interval(contig, pos, pos, mating_contig, mating_start, mating_end):
                stats["skipped_mating_type"] += 1
                continue

            rows.append(
                {
                    "ortholog_id": "",
                    "variant_id": f"{contig}:{pos}:{ref}>{alt}",
                    "contig": contig,
                    "start": pos,
                    "end": pos,
                    "midpoint": pos,
                    "length_bp": 1,
                    "ref_allele": ref,
                    "alt_allele": alt,
                    "group1_fixed_allele": allele_label(ref, alt, group1_gt),
                    "group2_fixed_allele": allele_label(ref, alt, group2_gt),
                    "group1_alt_count": int(sum(g1)),
                    "group2_alt_count": int(sum(g2)),
                    "group1_present_count": int(sum(g1)),
                    "group1_absent_count": len(group1_samples) - int(sum(g1)),
                    "group1_n_samples": len(group1_samples),
                    "group2_n_samples": len(group2_samples),
                    "frequency_category": fixed_category(group1_gt, group2_gt),
                    "primary_set": "fixed_snp_between_group1_group2",
                }
            )

    targets = pd.DataFrame(rows)
    if not targets.empty:
        targets = targets.sort_values(["contig", "midpoint", "variant_id"]).reset_index(drop=True)
    stats["fixed_snp_candidates_before_spacing"] = int(len(targets))

    if min_spacing > 0 and not targets.empty:
        keep_indices = []
        last_kept_by_contig: dict[str, int] = {}
        for idx, row in targets.iterrows():
            contig = str(row["contig"])
            pos = int(row["midpoint"])
            last_pos = last_kept_by_contig.get(contig)
            if last_pos is None or pos - last_pos >= min_spacing:
                keep_indices.append(idx)
                last_kept_by_contig[contig] = pos
        stats["skipped_min_spacing"] = int(len(targets) - len(keep_indices))
        targets = targets.loc[keep_indices].reset_index(drop=True)

    if not targets.empty:
        targets["ortholog_id"] = [f"fixed_snp_{idx:08d}" for idx in range(1, len(targets) + 1)]
    stats["targets_retained"] = int(len(targets))
    return targets, stats


def write_bed(targets: pd.DataFrame, output_bed: str | Path) -> None:
    columns = ["contig", "bed_start", "end", "ortholog_id", "frequency_category", "midpoint"]
    if targets.empty:
        Path(output_bed).write_text("")
        return
    bed = targets.copy()
    bed["bed_start"] = (bed["start"].astype(int) - 1).clip(lower=0)
    bed[columns].to_csv(output_bed, sep="\t", header=False, index=False)


def main() -> None:
    args = parse_args()
    output_tsv = Path(args.output_tsv)
    output_bed = Path(args.output_bed)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    targets, stats = build_fixed_snp_targets(
        vcf=args.vcf,
        group1_samples=parse_csv(args.group1_samples),
        group2_samples=parse_csv(args.group2_samples),
        min_spacing=args.min_spacing,
        mating_contig=args.mating_contig,
        mating_start=args.mating_start,
        mating_end=args.mating_end,
    )
    targets.to_csv(output_tsv, sep="\t", index=False)
    write_bed(targets, output_bed)

    for key, value in stats.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()

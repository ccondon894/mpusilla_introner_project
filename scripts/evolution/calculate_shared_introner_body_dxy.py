#!/usr/bin/env python3
"""
Calculate introner-body Dxy for fixed shared loci from the active all-samples
classification.

The active all-samples Dxy workflow already computes flanking Dxy from
consensus alignments. This script adds the corresponding introner-body Dxy for
fixed-present loci shared by Group 1 and Group 2, splitting them into broad
ancestral and independent origin classes using cross_group_status from the
genotype matrix.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


ANCESTRAL_CROSS_GROUP = {
    "ancestral",
    "likely_ancestral",
    "ancestral_low_identity",
}
INDEPENDENT_CROSS_GROUP = {
    "independent",
    "likely_independent",
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate introner-body Dxy for active all-samples shared loci"
    )
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--consensus-dir", required=True)
    parser.add_argument("--qc-table", required=True)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--group1-samples", nargs="+", required=True)
    parser.add_argument("--group2-samples", nargs="+", required=True)
    parser.add_argument(
        "--matrix-flank-length",
        type=int,
        default=100,
        help="Flank length included in genotype-matrix start/end coordinates",
    )
    return parser.parse_args()


def is_valid_base(base):
    return base.upper() not in {"-", "N"}


def calculate_pairwise_counts(seq1, seq2):
    differences = 0
    sites = 0
    for a, b in zip(seq1, seq2):
        if is_valid_base(a) and is_valid_base(b):
            sites += 1
            if a.upper() != b.upper():
                differences += 1
    return differences, sites


def calculate_dxy_ros(group1_seqs, group2_seqs):
    if not group1_seqs or not group2_seqs:
        return None

    total_differences = 0
    total_sites = 0
    for seq1 in group1_seqs:
        for seq2 in group2_seqs:
            differences, sites = calculate_pairwise_counts(seq1, seq2)
            total_differences += differences
            total_sites += sites

    return total_differences / total_sites if total_sites > 0 else None


def broad_origin(cross_group_status):
    if cross_group_status in ANCESTRAL_CROSS_GROUP:
        return "ancestral"
    if cross_group_status in INDEPENDENT_CROSS_GROUP:
        return "independent"
    return "unclassified"


def load_body_consensus(consensus_dir, samples):
    sequences = {}
    for sample in samples:
        fasta = Path(consensus_dir) / f"{sample}.introner_body.consensus.fa"
        if not fasta.exists():
            raise FileNotFoundError(f"Missing introner-body consensus FASTA: {fasta}")
        for record in SeqIO.parse(fasta, "fasta"):
            name = record.id.split("::", 1)[0]
            ortholog_id = name.split("|", 1)[0]
            sequences[(sample, ortholog_id)] = str(record.seq).upper()
    return sequences


def load_qc_table(path):
    qc = pd.read_csv(path, sep="\t")
    required = {
        "sample",
        "ortholog_id",
        "qc_pass",
        "qc_fail_reason",
        "callable_fraction",
        "softclip_read_fraction",
        "low_concordance_fraction",
        "mean_major_allele_fraction",
        "depth_ratio",
        "consensus_source",
    }
    missing = required - set(qc.columns)
    if missing:
        raise ValueError(f"QC table missing columns: {', '.join(sorted(missing))}")
    qc["qc_pass"] = qc["qc_pass"].astype(str).str.lower().isin({"true", "1"})
    return qc.set_index(["sample", "ortholog_id"], drop=False)


def consensus_body_sequence(row, consensus):
    key = (row["sample"], row["ortholog_id"])
    seq = consensus.get(key)
    if seq is None:
        return None
    if row.get("orientation", "") == "reverse":
        seq = str(Seq(seq).reverse_complement()).upper()
    return seq


def run_mafft(ortholog_id, sequences, alignment_dir):
    alignment_dir = Path(alignment_dir)
    alignment_dir.mkdir(parents=True, exist_ok=True)
    unaligned = alignment_dir / f"{ortholog_id}.introner_body.fa"
    aligned = alignment_dir / f"{ortholog_id}.introner_body.mafft.fa"

    with open(unaligned, "w") as handle:
        for sample, seq in sequences:
            handle.write(f">{sample}\n{seq}\n")

    result = subprocess.run(
        [
            "mafft",
            "--adjustdirection",
            "--maxiterate",
            "1000",
            "--globalpair",
            "--quiet",
            str(unaligned),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode != 0:
        print(
            f"MAFFT failed for {ortholog_id}: {result.stderr[:200]}",
            file=sys.stderr,
        )
        return None

    with open(aligned, "w") as handle:
        handle.write(result.stdout)
    unaligned.unlink(missing_ok=True)
    return aligned


def split_aligned_sequences(aligned_path, group1_samples, group2_samples):
    group1_set = set(group1_samples)
    group2_set = set(group2_samples)
    group1_seqs = []
    group2_seqs = []

    for record in SeqIO.parse(aligned_path, "fasta"):
        sample = record.id.replace("_R_", "")
        if sample in group1_set:
            group1_seqs.append(str(record.seq))
        elif sample in group2_set:
            group2_seqs.append(str(record.seq))

    return group1_seqs, group2_seqs


def main():
    args = parse_arguments()

    with open(args.classification) as handle:
        classification = json.load(handle)

    df = pd.read_csv(args.genotype_matrix, sep="\t")
    all_samples = args.group1_samples + args.group2_samples
    consensus = load_body_consensus(args.consensus_dir, all_samples)
    qc = load_qc_table(args.qc_table)

    rows = []
    skipped = 0

    for ortholog_id, info in sorted(classification.items()):
        if info["category"] != "group1_fixed_group2_fixed":
            continue

        cross_status = info.get("cross_group_status", "NA")
        origin = broad_origin(cross_status)
        if origin == "unclassified":
            skipped += 1
            continue

        locus = df[
            (df["ortholog_id"] == ortholog_id)
            & (df["sample"].isin(all_samples))
            & (df["presence"] == 1)
        ]

        sequences = []
        qc_reasons = []
        callable_fractions = []
        softclip_fractions = []
        low_concordance_fractions = []
        major_allele_fractions = []
        depth_ratios = []
        consensus_sources = []
        for sample in all_samples:
            sample_rows = locus[locus["sample"] == sample]
            if len(sample_rows) != 1:
                qc_reasons.append(f"{sample}:missing_genotype_row")
                continue
            row = sample_rows.iloc[0]
            qc_key = (sample, ortholog_id)
            if qc_key not in qc.index:
                qc_reasons.append(f"{sample}:missing_qc")
                continue
            qc_row = qc.loc[qc_key]
            callable_fractions.append(float(qc_row["callable_fraction"]))
            softclip_fractions.append(float(qc_row["softclip_read_fraction"]))
            low_concordance_fractions.append(float(qc_row["low_concordance_fraction"]))
            major_allele_fractions.append(float(qc_row["mean_major_allele_fraction"]))
            depth_ratios.append(float(qc_row["depth_ratio"]))
            consensus_sources.append(str(qc_row["consensus_source"]))
            if not bool(qc_row["qc_pass"]):
                qc_reasons.append(f"{sample}:{qc_row['qc_fail_reason']}")
                continue
            seq = consensus_body_sequence(row, consensus)
            expected_len = int(row["end"]) - int(row["start"]) - 2 * args.matrix_flank_length
            if seq is None:
                qc_reasons.append(f"{sample}:missing_consensus")
                continue
            if len(seq) != expected_len:
                qc_reasons.append(f"{sample}:length_mismatch_after_orientation")
                continue
            sequences.append((sample, seq))

        present_samples = {sample for sample, _ in sequences}
        if not set(args.group1_samples).issubset(present_samples):
            rows.append(
                {
                    "ortholog_id": ortholog_id,
                    "category": info["category"],
                    "cross_group_status": cross_status,
                    "ancestry_class": origin,
                    "within_group_status": info.get("within_group_status", ""),
                    "group1_present_count": info["group1_present_count"],
                    "group2_present_count": info["group2_present_count"],
                    "n_group1_seqs": 0,
                    "n_group2_seqs": 0,
                    "consensus_source": ",".join(sorted(set(consensus_sources))),
                    "n_consensus_sequences": len(sequences),
                    "n_qc_failed_present_carriers": len(all_samples) - len(sequences),
                    "qc_drop_reason": ";".join(qc_reasons) if qc_reasons else "missing_group1_consensus",
                    "mean_callable_fraction": (
                        float(np.mean(callable_fractions)) if callable_fractions else np.nan
                    ),
                    "max_softclip_read_fraction": (
                        max(softclip_fractions) if softclip_fractions else np.nan
                    ),
                    "max_low_concordance_fraction": (
                        max(low_concordance_fractions)
                        if low_concordance_fractions
                        else np.nan
                    ),
                    "min_mean_major_allele_fraction": (
                        min(major_allele_fractions) if major_allele_fractions else np.nan
                    ),
                    "max_depth_ratio": max(depth_ratios) if depth_ratios else np.nan,
                    "dxy_introner": None,
                }
            )
            skipped += 1
            continue
        if not set(args.group2_samples).issubset(present_samples):
            rows.append(
                {
                    "ortholog_id": ortholog_id,
                    "category": info["category"],
                    "cross_group_status": cross_status,
                    "ancestry_class": origin,
                    "within_group_status": info.get("within_group_status", ""),
                    "group1_present_count": info["group1_present_count"],
                    "group2_present_count": info["group2_present_count"],
                    "n_group1_seqs": 0,
                    "n_group2_seqs": 0,
                    "consensus_source": ",".join(sorted(set(consensus_sources))),
                    "n_consensus_sequences": len(sequences),
                    "n_qc_failed_present_carriers": len(all_samples) - len(sequences),
                    "qc_drop_reason": ";".join(qc_reasons) if qc_reasons else "missing_group2_consensus",
                    "mean_callable_fraction": (
                        float(np.mean(callable_fractions)) if callable_fractions else np.nan
                    ),
                    "max_softclip_read_fraction": (
                        max(softclip_fractions) if softclip_fractions else np.nan
                    ),
                    "max_low_concordance_fraction": (
                        max(low_concordance_fractions)
                        if low_concordance_fractions
                        else np.nan
                    ),
                    "min_mean_major_allele_fraction": (
                        min(major_allele_fractions) if major_allele_fractions else np.nan
                    ),
                    "max_depth_ratio": max(depth_ratios) if depth_ratios else np.nan,
                    "dxy_introner": None,
                }
            )
            skipped += 1
            continue

        aligned_path = run_mafft(ortholog_id, sequences, args.alignment_dir)
        if aligned_path is None:
            skipped += 1
            continue

        group1_seqs, group2_seqs = split_aligned_sequences(
            aligned_path,
            args.group1_samples,
            args.group2_samples,
        )
        dxy_introner = calculate_dxy_ros(group1_seqs, group2_seqs)

        rows.append(
            {
                "ortholog_id": ortholog_id,
                "category": info["category"],
                "cross_group_status": cross_status,
                "ancestry_class": origin,
                "within_group_status": info.get("within_group_status", ""),
                "group1_present_count": info["group1_present_count"],
                "group2_present_count": info["group2_present_count"],
                "n_group1_seqs": len(group1_seqs),
                "n_group2_seqs": len(group2_seqs),
                "consensus_source": ",".join(sorted(set(consensus_sources))),
                "n_consensus_sequences": len(sequences),
                "n_qc_failed_present_carriers": len(all_samples) - len(sequences),
                "qc_drop_reason": "PASS" if not qc_reasons else ";".join(qc_reasons),
                "mean_callable_fraction": (
                    float(np.mean(callable_fractions)) if callable_fractions else np.nan
                ),
                "max_softclip_read_fraction": (
                    max(softclip_fractions) if softclip_fractions else np.nan
                ),
                "max_low_concordance_fraction": (
                    max(low_concordance_fractions) if low_concordance_fractions else np.nan
                ),
                "min_mean_major_allele_fraction": (
                    min(major_allele_fractions) if major_allele_fractions else np.nan
                ),
                "max_depth_ratio": max(depth_ratios) if depth_ratios else np.nan,
                "dxy_introner": dxy_introner,
            }
        )

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    print(f"Processed {len(out)} fixed shared introner-body loci")
    print(f"Skipped {skipped} loci")
    if not out.empty:
        print(out["ancestry_class"].value_counts().to_string())


if __name__ == "__main__":
    main()

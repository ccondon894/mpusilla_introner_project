#!/usr/bin/env python3
"""Compare introner-body sequence decay in fixed and polymorphic introners."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats


MATING_CONTIG = "CCMP1545#0#scaffold_2"
MATING_START = 49808
MATING_END = 1730591
PRIMARY_CLASSES = ("fixed_present", "polymorphic")
LOW_IDENTITY_STATUS = "low_identity"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotype-matrix", required=True, type=Path)
    parser.add_argument("--consensus-dir", required=True, type=Path)
    parser.add_argument("--qc-table", required=True, type=Path)
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--per-locus", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--tests", required=True, type=Path)
    parser.add_argument("--plot", required=True, type=Path)
    parser.add_argument("--group1-samples", required=True, nargs="+")
    parser.add_argument("--group2-samples", required=True, nargs="+")
    parser.add_argument("--top-families", type=int, default=3)
    parser.add_argument("--matrix-flank-length", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=1545)
    parser.add_argument(
        "--max-loci",
        type=int,
        default=0,
        help="Optional smoke-test cap on candidate orthologs before alignment.",
    )
    return parser.parse_args()


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "within_group_status" not in df.columns:
        raise ValueError("Genotype matrix is missing required within_group_status column.")
    numeric = [
        "start",
        "end",
        "presence",
        "family",
        "group1_present_count",
        "group1_callable_count",
        "group1_n_samples",
        "group2_present_count",
        "group2_callable_count",
        "group2_n_samples",
    ]
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def overlaps(start: int, end: int, query_start: int, query_end: int) -> bool:
    return start < query_end and query_start < end


def mating_type_orthologs(df: pd.DataFrame) -> set[str]:
    ccmp = df[df["sample"] == "CCMP1545"]
    mask = (
        (ccmp["contig"] == MATING_CONTIG)
        & (ccmp["start"] < MATING_END)
        & (MATING_START < ccmp["end"])
    )
    return set(ccmp.loc[mask, "ortholog_id"])


def modal_present_family(present_rows: pd.DataFrame) -> pd.Series:
    usable = present_rows[present_rows["family"].notna() & (present_rows["family"] >= 0)]
    def mode_int(values: pd.Series) -> int | float:
        counts = Counter(int(v) for v in values)
        if not counts:
            return math.nan
        return counts.most_common(1)[0][0]

    return usable.groupby("ortholog_id")["family"].agg(mode_int)


def top_present_families(df: pd.DataFrame, top_n: int) -> list[int]:
    present = df[df["presence"] == 1].copy()
    family_by_ortholog = modal_present_family(present)
    counts = family_by_ortholog.value_counts().sort_values(ascending=False)
    return [int(v) for v in counts.head(top_n).index]


def classify_group(row: pd.Series, group: str) -> str:
    pattern = row[f"{group}_pattern"]
    present = int(row[f"{group}_present_count"])
    callable_count = int(row[f"{group}_callable_count"])
    n_samples = int(row[f"{group}_n_samples"])
    if pattern == "fixed_present":
        return "fixed_present"
    if callable_count == n_samples and 0 < present < n_samples:
        return "polymorphic"
    return "excluded"


def candidate_loci(
    df: pd.DataFrame,
    group_name: str,
    group_samples: list[str],
    families: list[int],
    mt_orthologs: set[str],
) -> pd.DataFrame:
    present = df[df["presence"] == 1].copy()
    family_by_ortholog = modal_present_family(present)
    one_per_ortholog = df.drop_duplicates("ortholog_id").copy()
    one_per_ortholog["analysis_family"] = one_per_ortholog["ortholog_id"].map(
        family_by_ortholog
    )
    one_per_ortholog = one_per_ortholog[
        one_per_ortholog["analysis_family"].isin(families)
        & ~one_per_ortholog["ortholog_id"].isin(mt_orthologs)
        & (
            one_per_ortholog["within_group_status"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            != LOW_IDENTITY_STATUS
        )
    ].copy()
    one_per_ortholog["analysis_class"] = one_per_ortholog.apply(
        lambda row: classify_group(row, group_name), axis=1
    )
    one_per_ortholog = one_per_ortholog[
        one_per_ortholog["analysis_class"].isin(PRIMARY_CLASSES)
    ].copy()
    one_per_ortholog["analysis_scope"] = (
        "group1_primary" if group_name == "group1" else "group2_secondary"
    )
    one_per_ortholog["analysis_group"] = group_name
    one_per_ortholog["analysis_samples"] = ",".join(group_samples)
    return one_per_ortholog


def load_body_consensus(consensus_dir: Path, samples: list[str]) -> dict[tuple[str, str], str]:
    sequences = {}
    for sample in samples:
        fasta = consensus_dir / f"{sample}.introner_body.consensus.fa"
        if not fasta.exists():
            raise FileNotFoundError(f"Missing introner-body consensus FASTA: {fasta}")
        for record in SeqIO.parse(fasta, "fasta"):
            name = record.id.split("::", 1)[0]
            ortholog_id = name.split("|", 1)[0]
            sequences[(sample, ortholog_id)] = str(record.seq).upper()
    return sequences


def load_qc_table(path: Path) -> pd.DataFrame:
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


def consensus_body_sequence(
    row: pd.Series,
    consensus: dict[tuple[str, str], str],
) -> str | None:
    key = (row["sample"], row["ortholog_id"])
    seq = consensus.get(key)
    if seq is None:
        return None
    if row.get("orientation", "") == "reverse":
        seq = str(Seq(seq).reverse_complement()).upper()
    return seq


def sanitize_name(value: str) -> str:
    return value.replace("/", "_").replace(" ", "_")


def write_singleton_alignment(path: Path, sequences: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for sample, seq in sequences:
            handle.write(f">{sample}\n{seq}\n")
    return path


def run_mafft(ortholog_id: str, scope: str, sequences: list[tuple[str, str]], outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = f"{sanitize_name(scope)}.{sanitize_name(ortholog_id)}"
    unaligned = outdir / f"{prefix}.introner_body.fa"
    aligned = outdir / f"{prefix}.introner_body.mafft.fa"
    with unaligned.open("w") as handle:
        for sample, seq in sequences:
            handle.write(f">{sample}\n{seq}\n")

    if len(sequences) == 1:
        unaligned.replace(aligned)
        return aligned, ""

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
        return None, result.stderr[:500]
    aligned.write_text(result.stdout)
    unaligned.unlink(missing_ok=True)
    return aligned, ""


def is_valid_base(base: str) -> bool:
    return base.upper() not in {"-", "N"}


def pairwise_counts(seq1: str, seq2: str) -> tuple[int, int]:
    differences = 0
    sites = 0
    for a, b in zip(seq1, seq2):
        if is_valid_base(a) and is_valid_base(b):
            sites += 1
            if a.upper() != b.upper():
                differences += 1
    return differences, sites


def calculate_pi(sequences: list[str]) -> tuple[float | None, int, int, int]:
    if len(sequences) < 2:
        return None, 0, 0, 0
    total_diff = 0
    total_sites = 0
    comparisons = 0
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            diff, sites = pairwise_counts(sequences[i], sequences[j])
            total_diff += diff
            total_sites += sites
            comparisons += 1
    pi = total_diff / total_sites if total_sites else None
    return pi, total_diff, total_sites, comparisons


def alignment_column_metrics(sequences: list[str]) -> tuple[int, int, int, int, float]:
    if not sequences:
        return 0, 0, 0, 0, math.nan
    aln_len = len(sequences[0])
    segregating = 0
    valid_columns = 0
    gap_chars = 0
    total_chars = aln_len * len(sequences)
    for col_idx in range(aln_len):
        column = [seq[col_idx].upper() for seq in sequences]
        gap_chars += sum(base == "-" for base in column)
        valid = [base for base in column if is_valid_base(base)]
        if valid:
            valid_columns += 1
        if len(set(valid)) > 1:
            segregating += 1
    gap_fraction = gap_chars / total_chars if total_chars else math.nan
    return aln_len, valid_columns, segregating, gap_chars, gap_fraction


def alignment_sample_names_and_sequences(aligned_path: Path) -> tuple[list[str], list[str]]:
    names = []
    seqs = []
    for record in SeqIO.parse(aligned_path, "fasta"):
        names.append(record.id.replace("_R_", ""))
        seqs.append(str(record.seq).upper())
    return names, seqs


def analyze_locus(
    locus: pd.Series,
    df: pd.DataFrame,
    group_samples: list[str],
    consensus: dict[tuple[str, str], str],
    qc: pd.DataFrame,
    args: argparse.Namespace,
) -> dict[str, object]:
    ortholog_id = locus["ortholog_id"]
    present_rows = df[
        (df["ortholog_id"] == ortholog_id)
        & (df["sample"].isin(group_samples))
        & (df["presence"] == 1)
    ].copy()

    sequences = []
    skipped_samples = []
    failed_samples = []
    qc_reasons = []
    callable_fractions = []
    softclip_fractions = []
    low_concordance_fractions = []
    major_allele_fractions = []
    depth_ratios = []
    consensus_sources = []
    expected_lens = []
    for _, row in present_rows.iterrows():
        sample = row["sample"]
        expected_len = int(row["end"]) - int(row["start"]) - 2 * args.matrix_flank_length
        expected_lens.append(expected_len)
        qc_key = (sample, ortholog_id)
        if qc_key not in qc.index:
            skipped_samples.append(sample)
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
            failed_samples.append(sample)
            qc_reasons.append(f"{sample}:{qc_row['qc_fail_reason']}")
            continue
        seq = consensus_body_sequence(row, consensus)
        if seq is None:
            skipped_samples.append(sample)
            qc_reasons.append(f"{sample}:missing_consensus")
            continue
        if len(seq) != expected_len:
            failed_samples.append(sample)
            qc_reasons.append(f"{sample}:length_mismatch_after_orientation")
            continue
        sequences.append((sample, seq))

    qc_complete = len(sequences) == len(present_rows) and not failed_samples and not skipped_samples
    if not qc_complete:
        sequences = []

    aligned_path = None
    mafft_error = ""
    aligned_names = []
    aligned_seqs = []
    if sequences:
        aligned_path, mafft_error = run_mafft(
            ortholog_id,
            locus["analysis_scope"],
            sequences,
            args.alignment_dir,
        )
        if aligned_path is not None:
            aligned_names, aligned_seqs = alignment_sample_names_and_sequences(aligned_path)

    pi, pairwise_differences, pairwise_sites, comparisons = calculate_pi(aligned_seqs)
    aln_len, valid_cols, seg_sites, gap_chars, gap_fraction = alignment_column_metrics(
        aligned_seqs
    )
    segregating_per_valid_bp = seg_sites / valid_cols if valid_cols else math.nan

    return {
        "analysis_scope": locus["analysis_scope"],
        "analysis_group": locus["analysis_group"],
        "analysis_class": locus["analysis_class"],
        "ortholog_id": ortholog_id,
        "family": int(locus["analysis_family"]),
        "group_present_count": int(locus[f"{locus['analysis_group']}_present_count"]),
        "group_callable_count": int(locus[f"{locus['analysis_group']}_callable_count"]),
        "group_n_samples": int(locus[f"{locus['analysis_group']}_n_samples"]),
        "n_present_rows": len(present_rows),
        "n_sequences": len(sequences),
        "n_consensus_sequences": len(sequences),
        "n_qc_failed_present_carriers": len(failed_samples) + len(skipped_samples),
        "qc_drop_reason": "PASS" if qc_complete else ";".join(qc_reasons),
        "consensus_source": ",".join(sorted(set(consensus_sources))) if consensus_sources else "",
        "mean_callable_fraction": np.mean(callable_fractions) if callable_fractions else math.nan,
        "max_softclip_read_fraction": max(softclip_fractions) if softclip_fractions else math.nan,
        "max_low_concordance_fraction": (
            max(low_concordance_fractions) if low_concordance_fractions else math.nan
        ),
        "min_mean_major_allele_fraction": (
            min(major_allele_fractions) if major_allele_fractions else math.nan
        ),
        "max_depth_ratio": max(depth_ratios) if depth_ratios else math.nan,
        "n_aligned_sequences": len(aligned_seqs),
        "present_samples": ",".join(sample for sample, _ in sequences),
        "skipped_samples": ",".join(skipped_samples),
        "body_len_min": min(expected_lens) if expected_lens else math.nan,
        "body_len_max": max(expected_lens) if expected_lens else math.nan,
        "alignment_length": aln_len,
        "valid_alignment_columns": valid_cols,
        "pairwise_comparisons": comparisons,
        "pairwise_differences": pairwise_differences,
        "pairwise_valid_sites": pairwise_sites,
        "pi_introner_body": pi,
        "segregating_body_sites": seg_sites,
        "segregating_sites_per_valid_bp": segregating_per_valid_bp,
        "gap_chars": gap_chars,
        "gap_fraction": gap_fraction,
        "alignment_path": "" if aligned_path is None else str(aligned_path),
        "mafft_error": mafft_error,
    }


def summarize(per_locus: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [
        ("all_top3_families", per_locus),
        *[
            (f"family_{family}", sub)
            for family, sub in per_locus.groupby("family", sort=True)
        ],
    ]
    for family_label, fam_df in groups:
        for keys, sub in fam_df.groupby(["analysis_scope", "analysis_class"], dropna=False):
            scope, analysis_class = keys
            pi_values = sub["pi_introner_body"].dropna()
            rows.append(
                {
                    "family_group": family_label,
                    "analysis_scope": scope,
                    "analysis_class": analysis_class,
                    "n_loci": len(sub),
                    "n_loci_with_pi": len(pi_values),
                    "median_pi_introner_body": pi_values.median(),
                    "mean_pi_introner_body": pi_values.mean(),
                    "median_segregating_sites_per_valid_bp": sub[
                        "segregating_sites_per_valid_bp"
                    ].dropna().median(),
                    "median_gap_fraction": sub["gap_fraction"].dropna().median(),
                    "median_n_aligned_sequences": sub["n_aligned_sequences"].median(),
                }
            )
    return pd.DataFrame(rows)


def mean_difference(df: pd.DataFrame, column: str) -> float:
    fixed = df.loc[df["analysis_class"] == "fixed_present", column].dropna()
    poly = df.loc[df["analysis_class"] == "polymorphic", column].dropna()
    return float(poly.mean() - fixed.mean())


def permutation_test(
    df: pd.DataFrame,
    column: str,
    rng: np.random.Generator,
    n_perm: int,
) -> dict[str, float | int]:
    data = df[df["analysis_class"].isin(PRIMARY_CLASSES) & df[column].notna()].copy()
    data = data.reset_index(drop=True)
    if data["analysis_class"].nunique() < 2:
        return {
            "observed_polymorphic_minus_fixed": math.nan,
            "permutation_p_two_sided": math.nan,
            "permutation_p_polymorphic_lower": math.nan,
            "permutations": 0,
            "null_mean": math.nan,
            "null_sd": math.nan,
        }

    observed = mean_difference(data, column)
    labels = data["analysis_class"].to_numpy().copy()
    values = data[column].to_numpy()
    strata = [idx.to_numpy() for _, idx in data.groupby("family").groups.items()]
    null = np.empty(n_perm)
    for idx_perm in range(n_perm):
        perm_labels = labels.copy()
        for idx in strata:
            perm_labels[idx] = rng.permutation(perm_labels[idx])
        fixed = values[perm_labels == "fixed_present"]
        poly = values[perm_labels == "polymorphic"]
        null[idx_perm] = poly.mean() - fixed.mean()

    p_two = (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    p_lower = (np.sum(null <= observed) + 1) / (len(null) + 1)
    return {
        "observed_polymorphic_minus_fixed": observed,
        "permutation_p_two_sided": p_two,
        "permutation_p_polymorphic_lower": p_lower,
        "permutations": len(null),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
    }


def mann_whitney(df: pd.DataFrame, column: str) -> dict[str, float | int]:
    fixed = df.loc[df["analysis_class"] == "fixed_present", column].dropna()
    poly = df.loc[df["analysis_class"] == "polymorphic", column].dropna()
    if len(fixed) == 0 or len(poly) == 0:
        return {
            "fixed_n": len(fixed),
            "polymorphic_n": len(poly),
            "mannwhitney_u": math.nan,
            "mannwhitney_p_two_sided": math.nan,
            "mannwhitney_p_polymorphic_lower": math.nan,
        }
    two = stats.mannwhitneyu(poly, fixed, alternative="two-sided")
    lower = stats.mannwhitneyu(poly, fixed, alternative="less")
    return {
        "fixed_n": len(fixed),
        "polymorphic_n": len(poly),
        "mannwhitney_u": float(two.statistic),
        "mannwhitney_p_two_sided": float(two.pvalue),
        "mannwhitney_p_polymorphic_lower": float(lower.pvalue),
    }


def tests(per_locus: pd.DataFrame, rng: np.random.Generator, n_perm: int) -> pd.DataFrame:
    rows = []
    metrics = [
        "pi_introner_body",
        "segregating_sites_per_valid_bp",
        "gap_fraction",
    ]
    for scope, scope_df in per_locus.groupby("analysis_scope", sort=True):
        groups = [("all_top3_families", scope_df)]
        groups.extend((f"family_{family}", sub) for family, sub in scope_df.groupby("family"))
        for family_group, sub in groups:
            for metric in metrics:
                perm = (
                    permutation_test(sub, metric, rng, n_perm)
                    if family_group == "all_top3_families"
                    else {
                        "observed_polymorphic_minus_fixed": mean_difference(
                            sub[sub[metric].notna()], metric
                        )
                        if sub["analysis_class"].nunique() == 2
                        else math.nan,
                        "permutation_p_two_sided": math.nan,
                        "permutation_p_polymorphic_lower": math.nan,
                        "permutations": 0,
                        "null_mean": math.nan,
                        "null_sd": math.nan,
                    }
                )
                mw = mann_whitney(sub, metric)
                rows.append(
                    {
                        "analysis_scope": scope,
                        "family_group": family_group,
                        "metric": metric,
                        **perm,
                        **mw,
                    }
                )
    return pd.DataFrame(rows)


def plot_results(per_locus: pd.DataFrame, path: Path) -> None:
    plot_df = per_locus[
        (per_locus["analysis_scope"] == "group1_primary")
        & per_locus["analysis_class"].isin(PRIMARY_CLASSES)
    ].copy()
    families = sorted(plot_df["family"].dropna().unique())
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4), sharey=False)
    metrics = [
        ("pi_introner_body", "Body pi"),
        ("segregating_sites_per_valid_bp", "Segregating sites / valid bp"),
        ("gap_fraction", "Gap fraction"),
    ]
    colors = {"fixed_present": "#4C78A8", "polymorphic": "#F58518"}
    for ax, (metric, title) in zip(axes, metrics):
        positions = []
        labels = []
        data = []
        facecolors = []
        pos = 0
        for family in families:
            for cls in PRIMARY_CLASSES:
                values = plot_df.loc[
                    (plot_df["family"] == family) & (plot_df["analysis_class"] == cls),
                    metric,
                ].dropna()
                positions.append(pos)
                labels.append(f"F{int(family)}\n{cls.replace('_', ' ')}")
                data.append(values.to_numpy())
                facecolors.append(colors[cls])
                pos += 1
            pos += 0.6
        box = ax.boxplot(data, positions=positions, widths=0.45, patch_artist=True, showfliers=False)
        for patch, color in zip(box["boxes"], facecolors):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
        ax.set_title(title)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    handles = [
        plt.Line2D([0], [0], color=colors[cls], lw=6, alpha=0.5, label=cls.replace("_", " "))
        for cls in PRIMARY_CLASSES
    ]
    axes[0].legend(handles=handles, frameon=False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def validate(
    per_locus: pd.DataFrame,
    top_families: list[int],
    low_identity_orthologs: set[str],
) -> None:
    if per_locus.empty:
        raise ValueError("No per-locus rows were produced.")
    if not set(per_locus["family"].unique()).issubset(set(top_families)):
        raise ValueError("Output contains families outside the selected top families.")
    if (per_locus["body_len_min"] <= 0).any() or (per_locus["body_len_max"] <= 0).any():
        raise ValueError("Output contains non-positive body lengths.")
    retained_low_identity = set(per_locus["ortholog_id"]) & low_identity_orthologs
    if retained_low_identity:
        raise ValueError(
            "Output contains low-identity orthologs: "
            + ", ".join(sorted(retained_low_identity))
        )
    group1 = per_locus[per_locus["analysis_scope"] == "group1_primary"]
    fixed = group1[group1["analysis_class"] == "fixed_present"]
    poly = group1[group1["analysis_class"] == "polymorphic"]
    if not (fixed["group_present_count"] == 11).all():
        raise ValueError("Group 1 fixed-present rows do not all have count 11.")
    if not (
        (poly["group_callable_count"] == 11)
        & (poly["group_present_count"] > 0)
        & (poly["group_present_count"] < 11)
    ).all():
        raise ValueError("Group 1 polymorphic rows do not match complete-call definition.")


def main() -> None:
    args = parse_args()
    if args.permutations < 1:
        raise ValueError("--permutations must be at least 1.")

    df = read_matrix(args.genotype_matrix)
    low_identity_orthologs = set(
        df.loc[
            df["within_group_status"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .eq(LOW_IDENTITY_STATUS),
            "ortholog_id",
        ]
    )
    mt_orthologs = mating_type_orthologs(df)
    top_families = top_present_families(df, args.top_families)
    all_samples = sorted(set(args.group1_samples + args.group2_samples))
    consensus = load_body_consensus(args.consensus_dir, all_samples)
    qc = load_qc_table(args.qc_table)

    group1_loci = candidate_loci(
        df, "group1", args.group1_samples, top_families, mt_orthologs
    )
    group2_loci = candidate_loci(
        df, "group2", args.group2_samples, top_families, mt_orthologs
    )
    candidates = pd.concat([group1_loci, group2_loci], ignore_index=True)
    if args.max_loci:
        candidates = (
            candidates.sort_values(["analysis_scope", "family", "analysis_class", "ortholog_id"])
            .groupby(["analysis_scope", "family", "analysis_class"], group_keys=False)
            .head(args.max_loci)
        )

    group_samples_by_scope = {
        "group1_primary": args.group1_samples,
        "group2_secondary": args.group2_samples,
    }
    rows = []
    for _, locus in candidates.sort_values(["analysis_scope", "family", "ortholog_id"]).iterrows():
        rows.append(
            analyze_locus(
                locus,
                df,
                group_samples_by_scope[locus["analysis_scope"]],
                consensus,
                qc,
                args,
            )
        )

    per_locus = pd.DataFrame(rows)
    validate(per_locus, top_families, low_identity_orthologs)
    summary = summarize(per_locus)
    rng = np.random.default_rng(args.seed)
    test_df = tests(per_locus, rng, args.permutations)

    args.per_locus.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.tests.parent.mkdir(parents=True, exist_ok=True)
    per_locus.to_csv(args.per_locus, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)
    test_df.to_csv(args.tests, sep="\t", index=False)
    plot_results(per_locus, args.plot)

    print(f"Top families: {','.join(map(str, top_families))}")
    print(
        "Excluded low-identity orthologs from body pi analysis: "
        f"{len(low_identity_orthologs)} matrix-wide"
    )
    print(f"Wrote {len(per_locus)} per-locus rows to {args.per_locus}")
    print(f"Wrote summary to {args.summary}")
    print(f"Wrote tests to {args.tests}")


if __name__ == "__main__":
    main()

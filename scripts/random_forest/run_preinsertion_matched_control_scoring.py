#!/usr/bin/env python3
"""Score RCC1749 pre-insertion sites against GC/same-gene matched controls."""

import argparse
import random
import shlex
import sys
from pathlib import Path

import pandas as pd

from build_preinsertion_features import (
    add_window_features,
    build_exclusion_index,
    overlaps_interval,
    parse_windows,
    read_exons,
    read_fai,
    read_reverse_oriented_contigs,
)
from introner_rf.config import ModelConfig
from introner_rf.data import prepare_dataset
from introner_rf.feature_selection import get_feature_cols
from introner_rf.matching import build_matched_control_dataset
from introner_rf.modeling import train_model
from introner_rf.validation import evaluate_predictions


DEFAULT_FEATURE_SETS = ["composition", "kmer_left_right", "all"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--positive-matrix",
        default=(
            "rf_results/ccmp1545_train_rcc1749_preinsertion_scoring_mask3_orientnorm_exclude28/"
            "rcc1749_preinsertion_sequence_features.mask3.orientnorm.exclude28.tsv"
        ),
        help="Pre-insertion feature matrix; label=1 rows are used as cases.",
    )
    parser.add_argument(
        "--train-matrix",
        default="rf_results/ccmp1545_present_absent_direct.mask3.w25_50_100.codon_microc_tpm.tsv",
    )
    parser.add_argument("--genotype-matrix", default="../../results/genotyping/genotype_matrix.final.tsv")
    parser.add_argument("--gtf", default="../../results/annotations/RCC1749.gtf")
    parser.add_argument("--fa", default="../../results/assemblies/RCC1749.vg_paths.fa")
    parser.add_argument(
        "--mummer-coords",
        default="../../results/genome_alignment/mummer/CCMP1545_vs_RCC1749.coords",
    )
    parser.add_argument("--orientation-threshold", type=float, default=0.8)
    parser.add_argument("--exclude-contig", action="append", default=["RCC1749#0#intronerless_contig_28"])
    parser.add_argument("--sample", default="RCC1749")
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--junction-mask-bp", type=int, default=3)
    parser.add_argument("--feature-sets", default=",".join(DEFAULT_FEATURE_SETS))
    parser.add_argument("--gc-col", default="gc_combined_100")
    parser.add_argument("--gc-tolerance", type=float, default=0.05)
    parser.add_argument("--same-gene-candidates-per-positive", type=int, default=20)
    parser.add_argument("--global-candidates-per-positive", type=int, default=20)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        default="rf_results/ccmp1545_train_rcc1749_preinsertion_matched_controls_mask3_orientnorm_exclude28",
    )
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def write_config(args, output_dir, windows, feature_sets):
    lines = [
        "RCC1749 pre-insertion matched-control scoring",
        "",
        "Purpose:",
        "Train the CCMP1545 present/absent RF and score RCC1749 empty orthologs against exonic controls matched by flank GC, with an additional same-gene+GC matched test.",
        "",
        "Command:",
        " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "",
        "Configuration:",
        f"positive_matrix: {args.positive_matrix}",
        f"train_matrix: {args.train_matrix}",
        f"genotype_matrix: {args.genotype_matrix}",
        f"gtf: {args.gtf}",
        f"fa: {args.fa}",
        f"mummer_coords: {args.mummer_coords}",
        f"exclude_contig: {','.join(args.exclude_contig)}",
        f"windows: {','.join(map(str, windows))}",
        f"junction_mask_bp: {args.junction_mask_bp}",
        f"feature_sets: {','.join(feature_sets)}",
        f"gc_col: {args.gc_col}",
        f"gc_tolerance: {args.gc_tolerance}",
        f"same_gene_candidates_per_positive: {args.same_gene_candidates_per_positive}",
        f"global_candidates_per_positive: {args.global_candidates_per_positive}",
        f"n_estimators: {args.n_estimators}",
        "class_weight: balanced",
        f"random_state: {args.random_state}",
        "",
        "Matched test sets:",
        "gc_only: controls selected from the global exonic candidate pool by nearest GC within tolerance.",
        "same_gene_gc: controls selected from the same GTF gene as the empty ortholog, also constrained by GC tolerance.",
    ]
    (output_dir / "config.txt").write_text("\n".join(lines) + "\n")


def build_gene_index(exons_by_contig):
    gene_to_exons = {}
    weighted_exons = []
    weights = []
    for contig, exons in exons_by_contig.items():
        for exon in exons:
            gene_to_exons.setdefault(exon["gene"], []).append(exon)
            width = max(0, exon["end"] - exon["start"])
            if width:
                weighted_exons.append(exon)
                weights.append(width)
    return gene_to_exons, weighted_exons, weights


def site_bounds(exon, max_window, fasta_sizes):
    low = max(int(exon["start"]) + max_window, 0)
    high = min(int(exon["end"]) - max_window, fasta_sizes.get(exon["contig"], 0))
    return low, high


def sample_site_from_exons(rng, exons, fasta_sizes, max_window):
    eligible = []
    weights = []
    for exon in exons:
        low, high = site_bounds(exon, max_window, fasta_sizes)
        if high <= low:
            continue
        eligible.append((exon, low, high))
        weights.append(high - low)
    if not eligible:
        return None, None
    exon, low, high = rng.choices(eligible, weights=weights, k=1)[0]
    return exon, rng.randrange(low, high)


def make_control_row(exon, site, sample, source, control_idx, reverse_contigs):
    return {
        "sequence_id": f"{exon['contig']}:{site}-{site}:{source}:control_{control_idx:08d}",
        "ortholog_id": "",
        "sample": sample,
        "source": source,
        "contig": exon["contig"],
        "start": site,
        "end": site,
        "site": site,
        "feature_start": site,
        "feature_end": site,
        "gene": exon["gene"],
        "transcript_id": exon["transcript_id"],
        "strand": exon["strand"],
        "family": "exonic_control",
        "absent_span": pd.NA,
        "exon_start": exon["start"],
        "exon_end": exon["end"],
        "reverse_oriented_contig": int(exon["contig"] in reverse_contigs),
        "label": 0,
    }


def generate_control_candidates(
    positives,
    exons_by_contig,
    fasta_sizes,
    exclusion_index,
    reverse_contigs,
    max_window,
    sample,
    same_gene_per_positive,
    global_per_positive,
    random_state,
):
    rng = random.Random(random_state)
    gene_to_exons, weighted_exons, exon_weights = build_gene_index(exons_by_contig)
    used_sites = set()
    controls = []
    skipped_same_gene = 0
    skipped_global = 0
    control_idx = 0

    for case in positives.to_dict("records"):
        gene = str(case.get("gene", ""))
        gene_exons = gene_to_exons.get(gene, [])

        made = 0
        attempts = max(same_gene_per_positive * 100, 100)
        for _ in range(attempts):
            if made >= same_gene_per_positive:
                break
            exon, site = sample_site_from_exons(rng, gene_exons, fasta_sizes, max_window)
            if exon is None:
                break
            site_key = (exon["contig"], site)
            if site_key in used_sites:
                continue
            if overlaps_interval(exclusion_index, exon["contig"], site - max_window, site + max_window):
                continue
            used_sites.add(site_key)
            control_idx += 1
            controls.append(
                make_control_row(
                    exon,
                    site,
                    sample,
                    "rcc1749_same_gene_exonic_control_candidate",
                    control_idx,
                    reverse_contigs,
                )
            )
            made += 1
        skipped_same_gene += same_gene_per_positive - made

        made = 0
        attempts = max(global_per_positive * 100, 100)
        for _ in range(attempts):
            if made >= global_per_positive:
                break
            exon = rng.choices(weighted_exons, weights=exon_weights, k=1)[0]
            low, high = site_bounds(exon, max_window, fasta_sizes)
            if high <= low:
                continue
            site = rng.randrange(low, high)
            site_key = (exon["contig"], site)
            if site_key in used_sites:
                continue
            if overlaps_interval(exclusion_index, exon["contig"], site - max_window, site + max_window):
                continue
            used_sites.add(site_key)
            control_idx += 1
            controls.append(
                make_control_row(
                    exon,
                    site,
                    sample,
                    "rcc1749_global_exonic_control_candidate",
                    control_idx,
                    reverse_contigs,
                )
            )
            made += 1
        skipped_global += global_per_positive - made

    stats = {
        "n_positive": len(positives),
        "same_gene_candidates_requested": len(positives) * same_gene_per_positive,
        "same_gene_candidates_made": sum(
            row["source"] == "rcc1749_same_gene_exonic_control_candidate"
            for row in controls
        ),
        "same_gene_candidates_skipped": skipped_same_gene,
        "global_candidates_requested": len(positives) * global_per_positive,
        "global_candidates_made": sum(
            row["source"] == "rcc1749_global_exonic_control_candidate"
            for row in controls
        ),
        "global_candidates_skipped": skipped_global,
    }
    return pd.DataFrame(controls), stats


def metadata_columns(df):
    preferred = [
        "sequence_id",
        "ortholog_id",
        "sample",
        "source",
        "gene",
        "family",
        "contig",
        "start",
        "end",
        "site",
        "feature_start",
        "feature_end",
        "absent_span",
        "label",
        "group_id",
        "reverse_oriented_contig",
        "matched_pair_id",
        "match_case_sequence_id",
        "match_case_gene",
        "match_gc_col",
        "match_gc_delta",
        "match_abs_gc_delta",
        "match_rank",
    ]
    return [col for col in preferred if col in df.columns]


def score_matched_dataset(train_df, test_df, feature_sets, windows, model_config):
    rows = []
    prediction_frames = []
    feature_rows = []
    for feature_set in feature_sets:
        train_cols = get_feature_cols(train_df, windows, feature_set=feature_set)
        test_numeric = {
            col for col in test_df.columns
            if pd.api.types.is_numeric_dtype(test_df[col])
        }
        feature_cols = [col for col in train_cols if col in test_numeric]
        if not feature_cols:
            raise ValueError(f"No shared numeric features selected for {feature_set!r}")

        clf = train_model(train_df[feature_cols], train_df["label"], model_config=model_config)
        y_true = test_df["label"]
        y_prob = clf.predict_proba(test_df[feature_cols])[:, 1]
        y_pred = clf.predict(test_df[feature_cols])
        metrics = evaluate_predictions(y_true, y_pred, y_prob)
        positive_probs = y_prob[y_true.to_numpy() == 1]
        control_probs = y_prob[y_true.to_numpy() == 0]

        row = {
            "feature_set": feature_set,
            "n_train": len(train_df),
            "n_train_pos": int(train_df["label"].eq(1).sum()),
            "n_train_neg": int(train_df["label"].eq(0).sum()),
            "n_test": len(test_df),
            "n_test_pos": int(test_df["label"].eq(1).sum()),
            "n_test_neg": int(test_df["label"].eq(0).sum()),
            "n_train_selected_features": len(train_cols),
            "n_shared_features": len(feature_cols),
            "mean_prob_present_for_empty_orthologs": float(positive_probs.mean()),
            "mean_prob_present_for_matched_controls": float(control_probs.mean()),
            "median_prob_present_for_empty_orthologs": float(pd.Series(positive_probs).median()),
            "median_prob_present_for_matched_controls": float(pd.Series(control_probs).median()),
            **metrics,
        }
        rows.append(row)

        predictions = test_df[metadata_columns(test_df)].copy()
        predictions["feature_set"] = feature_set
        predictions["y_pred"] = y_pred
        predictions["prob_present"] = y_prob
        predictions["correct"] = y_true.to_numpy() == y_pred
        prediction_frames.append(predictions)
        feature_rows.extend(
            {"feature_set": feature_set, "feature": feature}
            for feature in feature_cols
        )
        print(
            f"  {feature_set}: ROC AUC={row['roc_auc']:.3f}, "
            f"mean positive/control prob={row['mean_prob_present_for_empty_orthologs']:.3f}/"
            f"{row['mean_prob_present_for_matched_controls']:.3f}"
        )

    return pd.DataFrame(rows), pd.concat(prediction_frames, ignore_index=True), pd.DataFrame(feature_rows)


def main():
    args = parse_args()
    windows = parse_windows(args.windows)
    window_labels = parse_csv(args.windows)
    feature_sets = parse_csv(args.feature_sets)
    max_window = max(windows)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_config(args, output_dir, windows, feature_sets)

    print("Loading positive pre-insertion matrix...")
    base_df = pd.read_csv(args.positive_matrix, sep="\t", low_memory=False)
    positives = base_df.loc[base_df["label"].eq(1)].copy()

    print("Loading RCC1749 annotation/FASTA context...")
    exons_by_contig = read_exons(args.gtf)
    excluded_contigs = set(args.exclude_contig)
    for contig in excluded_contigs:
        exons_by_contig.pop(contig, None)
    fasta_sizes = read_fai(args.fa)
    fasta_index = pd.NA
    from Bio import SeqIO

    fasta_index = SeqIO.index(args.fa, "fasta")
    reverse_contigs = read_reverse_oriented_contigs(
        args.mummer_coords,
        args.orientation_threshold,
    )
    genotype_df = pd.read_csv(args.genotype_matrix, sep="\t", low_memory=False)
    genotype_df["start"] = genotype_df["start"].astype(int)
    genotype_df["end"] = genotype_df["end"].astype(int)
    exclusion_index = build_exclusion_index(genotype_df, args.sample)

    print("Generating exonic control candidate pool...")
    controls, candidate_stats = generate_control_candidates(
        positives,
        exons_by_contig,
        fasta_sizes,
        exclusion_index,
        reverse_contigs,
        max_window,
        args.sample,
        args.same_gene_candidates_per_positive,
        args.global_candidates_per_positive,
        args.random_state,
    )

    print(f"Computing features for {len(controls)} control candidates...")
    control_features = [
        add_window_features(
            row,
            fasta_index,
            windows,
            junction_mask_bp=args.junction_mask_bp,
        )
        for row in controls.to_dict("records")
    ]
    controls = pd.concat([controls.reset_index(drop=True), pd.DataFrame(control_features)], axis=1)
    controls["group_id"] = controls["gene"].where(controls["gene"].astype(str).ne(""), controls["contig"])

    candidate_stats.update({
        "n_control_candidates": len(controls),
        "n_control_candidate_genes": controls["gene"].nunique(),
        "n_positive_genes": positives["gene"].nunique(),
        "gc_col": args.gc_col,
        "gc_tolerance": args.gc_tolerance,
    })
    pd.DataFrame([candidate_stats]).to_csv(
        output_dir / "control_candidate_summary.tsv",
        sep="\t",
        index=False,
    )

    combined = pd.concat([positives, controls], ignore_index=True, sort=False)
    combined = prepare_dataset(combined)
    train_df = prepare_dataset(pd.read_csv(args.train_matrix, sep="\t", low_memory=False))
    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=-1,
    )

    match_configs = [
        {
            "name": "gc_only",
            "same_gene": False,
            "same_contig": False,
            "candidate_sources": None,
        },
        {
            "name": "same_gene_gc",
            "same_gene": True,
            "same_contig": False,
            "candidate_sources": {"rcc1749_same_gene_exonic_control_candidate"},
        },
    ]

    all_summary = []
    for match_config in match_configs:
        condition = match_config["name"]
        print(f"Matching controls for {condition}...")
        match_input = combined
        if match_config["candidate_sources"] is not None:
            control_mask = (
                match_input["label"].eq(0)
                & match_input["source"].isin(match_config["candidate_sources"])
            )
            match_input = pd.concat(
                [match_input.loc[match_input["label"].eq(1)], match_input.loc[control_mask]],
                ignore_index=True,
                sort=False,
            )

        matched_df, match_stats = build_matched_control_dataset(
            match_input,
            gc_col=args.gc_col,
            controls_per_positive=1,
            gc_tolerance=args.gc_tolerance,
            same_contig=match_config["same_contig"],
            same_gene=match_config["same_gene"],
            drop_unannotated_positives=True,
            replace_controls=False,
            random_state=args.random_state,
        )
        matched_df.to_csv(output_dir / f"{condition}.matched_dataset.tsv", sep="\t", index=False)
        pd.DataFrame([match_stats]).to_csv(
            output_dir / f"{condition}.match_summary.tsv",
            sep="\t",
            index=False,
        )

        print(f"Scoring {condition} matched test set...")
        summary, predictions, features = score_matched_dataset(
            train_df,
            prepare_dataset(matched_df),
            feature_sets,
            window_labels,
            model_config,
        )
        summary.insert(0, "match_condition", condition)
        for key, value in match_stats.items():
            summary[f"match_{key}"] = value
        summary.to_csv(output_dir / f"{condition}.score_summary.tsv", sep="\t", index=False)
        predictions.insert(0, "match_condition", condition)
        predictions.to_csv(output_dir / f"{condition}.predictions.tsv", sep="\t", index=False)
        features.insert(0, "match_condition", condition)
        features.to_csv(output_dir / f"{condition}.features.tsv", sep="\t", index=False)
        all_summary.append(summary)

    pd.concat(all_summary, ignore_index=True).to_csv(
        output_dir / "summary.tsv",
        sep="\t",
        index=False,
    )
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

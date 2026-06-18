#!/usr/bin/env python3
"""Leave-one-introner-family-out present/absent RF tests in one sample."""

import argparse
import shlex
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import confusion_matrix

from introner_rf.config import ModelConfig
from introner_rf.data import prepare_dataset
from introner_rf.feature_selection import get_feature_cols
from introner_rf.modeling import train_model
from introner_rf.validation import evaluate_predictions


DEFAULT_FEATURE_SETS = [
    "composition",
    "kmer_left_right",
    "all",
    "microc",
    "kmer_left_right_microc_expression",
    "all_microc_codon_expression",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        default="rf_results/ccmp1545_present_absent_direct.mask3.w25_50_100.codon_microc_tpm.tsv",
    )
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--feature-sets", default=",".join(DEFAULT_FEATURE_SETS))
    parser.add_argument(
        "--negative-mode",
        choices=["same_gene", "random"],
        default="same_gene",
        help="How to choose test negatives for each held-out family.",
    )
    parser.add_argument("--min-test-pos", type=int, default=20)
    parser.add_argument("--min-test-neg", type=int, default=20)
    parser.add_argument("--random-test-neg-fraction", type=float, default=0.30)
    parser.add_argument("--max-random-test-neg-ratio", type=float, default=1.0)
    parser.add_argument(
        "--exclude-heldout-genes-from-train",
        action="store_true",
        help="Also remove loci in held-out positive genes from the training set.",
    )
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        default="rf_results/leave_one_family_out_ccmp1545_present_absent_mask3",
    )
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def metadata_columns(df):
    preferred = [
        "sequence_id",
        "ortholog_id",
        "sample",
        "gene",
        "gene_from_genotype_matrix",
        "family",
        "contig",
        "start",
        "end",
        "feature_start",
        "feature_end",
        "label",
        "sample_presence_status",
        "class_source",
        "group_id",
    ]
    return [col for col in preferred if col in df.columns]


def write_config(args, output_dir, feature_sets, windows):
    lines = [
        "Leave-one-family-out CCMP1545 present/absent RF test",
        "",
        "Purpose:",
        "Train on present introners from all but one introner family plus absent loci, then test on the held-out family's present loci and an independent negative set.",
        "",
        "Command:",
        " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "",
        "Configuration:",
        f"matrix: {args.matrix}",
        f"windows: {','.join(windows)}",
        f"feature_sets: {','.join(feature_sets)}",
        f"negative_mode: {args.negative_mode}",
        f"min_test_pos: {args.min_test_pos}",
        f"min_test_neg: {args.min_test_neg}",
        f"exclude_heldout_genes_from_train: {args.exclude_heldout_genes_from_train}",
        f"n_estimators: {args.n_estimators}",
        "class_weight: balanced",
        f"random_state: {args.random_state}",
        "",
        "Notes:",
        "Absences do not carry introner-family labels in this matrix. In same_gene mode, test negatives are absent loci in genes that contain the held-out family's positives; the remaining absences are used for training.",
    ]
    (output_dir / "config.txt").write_text("\n".join(lines) + "\n")


def choose_test_negatives(df, heldout_pos, args, family):
    negatives = df.loc[df["label"].eq(0)]
    if args.negative_mode == "same_gene":
        heldout_genes = set(heldout_pos["gene"].dropna().astype(str))
        return negatives.loc[negatives["gene"].astype(str).isin(heldout_genes)]

    max_by_fraction = int(len(negatives) * args.random_test_neg_fraction)
    max_by_ratio = int(len(heldout_pos) * args.max_random_test_neg_ratio)
    n_test_neg = min(len(negatives), max_by_fraction, max_by_ratio)
    if n_test_neg <= 0:
        return negatives.iloc[0:0]
    return negatives.sample(n=n_test_neg, random_state=args.random_state + int(family))


def evaluate_family_feature_set(train_df, test_df, windows, feature_set, model_config):
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
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    row = {
        "feature_set": feature_set,
        "n_features": len(feature_cols),
        "tn": int(cm[0, 0]),
        "fp": int(cm[0, 1]),
        "fn": int(cm[1, 0]),
        "tp": int(cm[1, 1]),
        **metrics,
    }

    predictions = test_df[metadata_columns(test_df)].copy()
    predictions["feature_set"] = feature_set
    predictions["y_pred"] = y_pred
    predictions["prob_present"] = y_prob
    predictions["correct"] = y_true.to_numpy() == y_pred
    return row, predictions, feature_cols


def main():
    args = parse_args()
    windows = parse_csv(args.windows)
    feature_sets = parse_csv(args.feature_sets)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_config(args, output_dir, feature_sets, windows)

    df = prepare_dataset(pd.read_csv(args.matrix, sep="\t", low_memory=False))
    df["family"] = df["family"].astype("string")
    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=-1,
    )

    family_counts = (
        df.groupby(["label", "family"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["label", "n"], ascending=[False, False])
    )
    family_counts.to_csv(output_dir / "family_counts.tsv", sep="\t", index=False)

    rows = []
    skipped = []
    prediction_frames = []
    feature_rows = []
    families = sorted(df.loc[df["label"].eq(1), "family"].dropna().unique(), key=lambda x: int(x))
    for family in families:
        heldout_pos = df.loc[df["label"].eq(1) & df["family"].eq(family)].copy()
        test_neg = choose_test_negatives(df, heldout_pos, args, family).copy()
        if len(heldout_pos) < args.min_test_pos or len(test_neg) < args.min_test_neg:
            skipped.append({
                "family": family,
                "n_test_pos": len(heldout_pos),
                "n_test_neg": len(test_neg),
                "reason": "below_minimum_test_counts",
            })
            continue

        test_idx = set(heldout_pos.index) | set(test_neg.index)
        train_df = df.drop(index=list(test_idx)).copy()
        train_df = train_df.loc[
            ~(train_df["label"].eq(1) & train_df["family"].eq(family))
        ].copy()
        if args.exclude_heldout_genes_from_train:
            heldout_genes = set(heldout_pos["gene"].dropna().astype(str))
            train_df = train_df.loc[~train_df["gene"].astype(str).isin(heldout_genes)].copy()

        test_df = pd.concat([heldout_pos, test_neg], ignore_index=False)
        for feature_set in feature_sets:
            row, predictions, feature_cols = evaluate_family_feature_set(
                train_df,
                test_df,
                windows,
                feature_set,
                model_config,
            )
            row.update({
                "heldout_family": family,
                "negative_mode": args.negative_mode,
                "n_train": len(train_df),
                "n_train_pos": int(train_df["label"].eq(1).sum()),
                "n_train_neg": int(train_df["label"].eq(0).sum()),
                "n_test": len(test_df),
                "n_test_pos": int(test_df["label"].eq(1).sum()),
                "n_test_neg": int(test_df["label"].eq(0).sum()),
            })
            rows.append(row)
            predictions["heldout_family"] = family
            prediction_frames.append(predictions)
            feature_rows.extend(
                {"heldout_family": family, "feature_set": feature_set, "feature": feature}
                for feature in feature_cols
            )
            print(
                f"family {family} {feature_set}: "
                f"ROC AUC={row['roc_auc']:.3f}, balanced accuracy={row['balanced_accuracy']:.3f}"
            )

    pd.DataFrame(rows).to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    if prediction_frames:
        pd.concat(prediction_frames, ignore_index=True).to_csv(
            output_dir / "predictions.tsv",
            sep="\t",
            index=False,
        )
    pd.DataFrame(feature_rows).to_csv(output_dir / "features.tsv", sep="\t", index=False)
    pd.DataFrame(skipped).to_csv(output_dir / "skipped_families.tsv", sep="\t", index=False)
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

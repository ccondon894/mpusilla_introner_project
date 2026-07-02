#!/usr/bin/env python3
"""Train a present/absent RF in one sample and evaluate in another sample."""

import argparse
import re
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
    "codon",
    "microc",
    "expression",
    "microc_expression",
    "kmer_left_right_expression",
    "kmer_left_right_microc",
    "kmer_left_right_microc_expression",
    "all_microc",
    "all_microc_codon",
    "all_microc_expression",
    "all_microc_codon_expression",
]

DELTA_COMPARISONS = [
    ("kmer_left_right", "kmer_left_right_microc"),
    ("kmer_left_right", "kmer_left_right_expression"),
    ("kmer_left_right_microc", "kmer_left_right_microc_expression"),
    ("all", "all_microc"),
    ("all_microc", "all_microc_codon"),
    ("all_microc", "all_microc_expression"),
    ("all_microc_codon", "all_microc_codon_expression"),
]

METRIC_COLS = ["roc_auc", "pr_auc", "balanced_accuracy", "accuracy", "f1"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-matrix", required=True)
    parser.add_argument("--test-matrix", required=True)
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument(
        "--feature-sets",
        default=",".join(DEFAULT_FEATURE_SETS),
        help="Comma-separated feature sets to train/evaluate.",
    )
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--output-prefix",
        default="rf_results/ccmp1545_train_rcc1749_test_present_absent",
    )
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_feature_names(df):
    rename = {}
    for col in df.columns:
        new_col = col
        new_col = re.sub(
            r"^microc_[^_]+_(\d+bp_\d+bp_.+)$",
            r"microc_sample_\1",
            new_col,
        )
        new_col = re.sub(
            r"^microc_[^_]+_(anchors|coverage)_cpm_(.+)$",
            r"microc_sample_\1_cpm_\2",
            new_col,
        )
        new_col = re.sub(
            r"^expr_tpm_\d+_(.+)$",
            r"expr_tpm_rep_\1",
            new_col,
        )
        rename[col] = new_col
    return df.rename(columns=rename)


def metadata_columns(df):
    preferred = [
        "sequence_id",
        "ortholog_id",
        "sample",
        "gene",
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


def evaluate_feature_set(train_df, test_df, windows, feature_set, model_config):
    train_cols = get_feature_cols(train_df, windows, feature_set=feature_set)
    test_numeric = {
        col for col in test_df.columns
        if pd.api.types.is_numeric_dtype(test_df[col])
    }
    feature_cols = [col for col in train_cols if col in test_numeric]
    if not feature_cols:
        raise ValueError(f"No shared numeric features selected for {feature_set!r}")

    clf = train_model(
        train_df[feature_cols],
        train_df["label"],
        model_config=model_config,
    )

    y_true = test_df["label"]
    y_prob = clf.predict_proba(test_df[feature_cols])[:, 1]
    y_pred = clf.predict(test_df[feature_cols])
    metrics = evaluate_predictions(y_true, y_pred, y_prob)
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])

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


def summarize_deltas(summary):
    rows = []
    indexed = summary.set_index("feature_set")
    for baseline, expanded in DELTA_COMPARISONS:
        if baseline not in indexed.index or expanded not in indexed.index:
            continue
        row = {
            "baseline_feature_set": baseline,
            "expanded_feature_set": expanded,
            "baseline_n_shared_features": int(indexed.loc[baseline, "n_shared_features"]),
            "expanded_n_shared_features": int(indexed.loc[expanded, "n_shared_features"]),
            "delta_n_shared_features": int(
                indexed.loc[expanded, "n_shared_features"]
                - indexed.loc[baseline, "n_shared_features"]
            ),
        }
        for metric in METRIC_COLS:
            row[f"delta_{metric}"] = indexed.loc[expanded, metric] - indexed.loc[baseline, metric]
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    windows = parse_csv(args.windows)
    feature_sets = parse_csv(args.feature_sets)
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("Loading matrices...")
    train_df = normalize_feature_names(pd.read_csv(args.train_matrix, sep="\t", low_memory=False))
    test_df = normalize_feature_names(pd.read_csv(args.test_matrix, sep="\t", low_memory=False))
    train_df = prepare_dataset(train_df)
    test_df = prepare_dataset(test_df)

    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=-1,
    )

    rows = []
    prediction_frames = []
    feature_rows = []
    for feature_set in feature_sets:
        print(f"Training/evaluating {feature_set}...")
        row, predictions, feature_cols = evaluate_feature_set(
            train_df,
            test_df,
            windows,
            feature_set,
            model_config,
        )
        rows.append(row)
        prediction_frames.append(predictions)
        feature_rows.extend(
            {"feature_set": feature_set, "feature": feature}
            for feature in feature_cols
        )
        print(
            f"  ROC AUC={row['roc_auc']:.3f}, "
            f"balanced accuracy={row['balanced_accuracy']:.3f}, "
            f"shared features={row['n_shared_features']}"
        )

    summary = pd.DataFrame(rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    features = pd.DataFrame(feature_rows)
    deltas = summarize_deltas(summary)

    summary_path = output_prefix.with_suffix(".summary.tsv")
    predictions_path = output_prefix.with_suffix(".predictions.tsv")
    features_path = output_prefix.with_suffix(".features.tsv")
    deltas_path = output_prefix.with_suffix(".feature_set_deltas.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    predictions.to_csv(predictions_path, sep="\t", index=False)
    features.to_csv(features_path, sep="\t", index=False)
    deltas.to_csv(deltas_path, sep="\t", index=False)
    print(summary)
    print(f"Wrote {summary_path}")
    print(f"Wrote {predictions_path}")
    print(f"Wrote {features_path}")
    print(f"Wrote {deltas_path}")


if __name__ == "__main__":
    main()

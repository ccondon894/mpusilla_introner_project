#!/usr/bin/env python3
"""Train a CCMP1545 present/absent RF and score RCC1749 pre-insertion sites."""

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


DEFAULT_FEATURE_SETS = ["composition", "kmer_left_right", "all"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--train-matrix",
        default="rf_results/ccmp1545_present_absent_direct.mask3.w25_50_100.codon_microc_tpm.tsv",
    )
    parser.add_argument(
        "--test-matrix",
        default="rf_results/rcc1749_preinsertion_sequence_features.tsv",
    )
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--feature-sets", default=",".join(DEFAULT_FEATURE_SETS))
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        default="rf_results/ccmp1545_train_rcc1749_preinsertion_scoring",
    )
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


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
    ]
    return [col for col in preferred if col in df.columns]


def write_config(args, output_dir, feature_sets, windows):
    lines = [
        "CCMP1545-trained scoring of RCC1749 pre-insertion candidates",
        "",
        "Purpose:",
        "Train on CCMP1545 present-vs-absent loci, then score RCC1749 loci that are empty orthologs of CCMP1545 introners against sampled RCC1749 exonic controls.",
        "",
        "Command:",
        " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "",
        "Configuration:",
        f"train_matrix: {args.train_matrix}",
        f"test_matrix: {args.test_matrix}",
        f"windows: {','.join(windows)}",
        f"feature_sets: {','.join(feature_sets)}",
        f"n_estimators: {args.n_estimators}",
        "class_weight: balanced",
        f"random_state: {args.random_state}",
        "",
        "Notes:",
        "This runner is most interpretable with sequence-only feature sets. If RCC1749 Micro-C or expression features are added to the pre-insertion matrix later, those feature sets can be supplied explicitly.",
    ]
    (output_dir / "config.txt").write_text("\n".join(lines) + "\n")


def evaluate_feature_set(train_df, test_df, windows, feature_set, model_config):
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
        "mean_prob_present_for_exonic_controls": float(control_probs.mean()),
        "median_prob_present_for_empty_orthologs": float(pd.Series(positive_probs).median()),
        "median_prob_present_for_exonic_controls": float(pd.Series(control_probs).median()),
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

    train_df = prepare_dataset(pd.read_csv(args.train_matrix, sep="\t", low_memory=False))
    test_df = prepare_dataset(pd.read_csv(args.test_matrix, sep="\t", low_memory=False))
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
            f"{feature_set}: ROC AUC={row['roc_auc']:.3f}, "
            f"balanced accuracy={row['balanced_accuracy']:.3f}, "
            f"mean positive/control prob={row['mean_prob_present_for_empty_orthologs']:.3f}/"
            f"{row['mean_prob_present_for_exonic_controls']:.3f}"
        )

    pd.DataFrame(rows).to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    pd.concat(prediction_frames, ignore_index=True).to_csv(
        output_dir / "predictions.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(feature_rows).to_csv(output_dir / "features.tsv", sep="\t", index=False)
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

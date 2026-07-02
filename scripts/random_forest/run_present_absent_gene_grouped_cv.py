#!/usr/bin/env python3
"""Run gene-grouped CV panels for present/absent introner RF models."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedGroupKFold

from introner_rf.config import ModelConfig
from introner_rf.data import load_prepared_dataset
from introner_rf.feature_selection import get_feature_cols
from introner_rf.modeling import train_model
from introner_rf.validation import evaluate_predictions, summarize_cv


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

DEFAULT_DELTA_COMPARISONS = [
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
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--feature-sets", default=",".join(DEFAULT_FEATURE_SETS))
    parser.add_argument("--group-col", default="group_id")
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


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
        "family",
    ]
    return [col for col in preferred if col in df.columns]


def balance_labels(df, random_state):
    counts = df["label"].value_counts()
    n = int(counts.min())
    return (
        df.groupby("label", group_keys=False)
        .sample(n=n, random_state=random_state)
        .sample(frac=1.0, random_state=random_state)
        .reset_index(drop=True)
    )


def cv_splits(df, group_col, n_splits, random_state):
    cv = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_state,
    )
    X_placeholder = np.zeros((len(df), 1))
    return list(cv.split(X_placeholder, df["label"], df[group_col]))


def run_feature_set_cv(df, panel, feature_set, feature_cols, splits, group_col, model_config):
    rows = []
    prediction_frames = []
    meta_cols = metadata_columns(df)

    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train_df = df.iloc[train_idx]
        test_df = df.iloc[test_idx]
        train_groups = set(train_df[group_col].astype(str))
        test_groups = set(test_df[group_col].astype(str))
        group_overlap = train_groups & test_groups
        if group_overlap:
            raise RuntimeError(
                f"{panel}/{feature_set}/fold {fold} has {len(group_overlap)} "
                f"overlapping {group_col} groups"
            )

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

        rows.append({
            "panel": panel,
            "feature_set": feature_set,
            "fold": fold,
            "n_features": len(feature_cols),
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_train_pos": int(train_df["label"].eq(1).sum()),
            "n_train_neg": int(train_df["label"].eq(0).sum()),
            "n_test_pos": int(y_true.eq(1).sum()),
            "n_test_neg": int(y_true.eq(0).sum()),
            "n_train_groups": len(train_groups),
            "n_test_groups": len(test_groups),
            "n_group_overlap": len(group_overlap),
            "tn": int(cm[0, 0]),
            "fp": int(cm[0, 1]),
            "fn": int(cm[1, 0]),
            "tp": int(cm[1, 1]),
            **metrics,
        })

        predictions = test_df[meta_cols].copy()
        predictions["panel"] = panel
        predictions["feature_set"] = feature_set
        predictions["fold"] = fold
        predictions["y_pred"] = y_pred
        predictions["prob_present"] = y_prob
        predictions["correct"] = y_true.to_numpy() == y_pred
        prediction_frames.append(predictions)

    return pd.DataFrame(rows), pd.concat(prediction_frames, ignore_index=True)


def comparison_row(panel, feature_set, df, feature_cols, summary):
    row = {
        "panel": panel,
        "feature_set": feature_set,
        "n_rows": len(df),
        "n_pos": int(df["label"].eq(1).sum()),
        "n_neg": int(df["label"].eq(0).sum()),
        "n_groups": df["group_id"].nunique(dropna=False),
        "n_features": len(feature_cols),
    }
    for metric in METRIC_COLS:
        row[metric] = summary.loc[metric, "mean"]
        row[f"{metric}_std"] = summary.loc[metric, "std"]
        row[f"{metric}_sem"] = summary.loc[metric, "sem"]
        row[f"{metric}_ci95"] = summary.loc[metric, "ci95"]
    return row


def summarize_deltas(metrics_df):
    rows = []
    for panel, panel_df in metrics_df.groupby("panel", sort=False):
        indexed = panel_df.set_index(["feature_set", "fold"])
        for baseline, expanded in DEFAULT_DELTA_COMPARISONS:
            if baseline not in indexed.index.get_level_values("feature_set"):
                continue
            if expanded not in indexed.index.get_level_values("feature_set"):
                continue
            base = indexed.loc[baseline]
            exp = indexed.loc[expanded]
            shared_folds = sorted(set(base.index) & set(exp.index))
            row = {
                "panel": panel,
                "baseline_feature_set": baseline,
                "expanded_feature_set": expanded,
                "n_folds": len(shared_folds),
            }
            for metric in METRIC_COLS:
                deltas = exp.loc[shared_folds, metric] - base.loc[shared_folds, metric]
                row[f"delta_{metric}_mean"] = deltas.mean()
                row[f"delta_{metric}_std"] = deltas.std()
                row[f"delta_{metric}_sem"] = deltas.std() / np.sqrt(len(deltas))
                row[f"delta_{metric}_ci95"] = 1.96 * row[f"delta_{metric}_sem"]
            rows.append(row)
    return pd.DataFrame(rows)


def windows_tag(windows):
    return "w" + "_".join(str(window) for window in windows)


def write_panel_outputs(out_dir, panel, feature_set, windows, metrics_df, predictions_df, summary):
    stem = f"{panel}_{feature_set}.{windows_tag(windows)}.gene_grouped_cv"
    metrics_df.to_csv(out_dir / f"{stem}_metrics.tsv", sep="\t", index=False)
    predictions_df.to_csv(out_dir / f"{stem}_predictions.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / f"{stem}_summary.tsv", sep="\t")


def run_panel(df, panel, feature_sets, windows, args, output_dir):
    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=-1,
    )
    splits = cv_splits(df, args.group_col, args.n_splits, args.random_state)
    comparison_rows = []
    metrics_frames = []
    prediction_frames = []
    feature_rows = []

    for feature_set in feature_sets:
        feature_cols = get_feature_cols(df, windows, feature_set=feature_set)
        if not feature_cols:
            raise ValueError(f"No features selected for {feature_set!r}")
        print(f"{panel}: {feature_set}: {len(feature_cols)} features")

        metrics_df, predictions_df = run_feature_set_cv(
            df,
            panel,
            feature_set,
            feature_cols,
            splits,
            args.group_col,
            model_config,
        )
        summary = summarize_cv(metrics_df)
        write_panel_outputs(
            output_dir,
            panel,
            feature_set,
            windows,
            metrics_df,
            predictions_df,
            summary,
        )
        comparison_rows.append(comparison_row(panel, feature_set, df, feature_cols, summary))
        metrics_frames.append(metrics_df)
        prediction_frames.append(predictions_df)
        feature_rows.extend(
            {"panel": panel, "feature_set": feature_set, "feature": feature}
            for feature in feature_cols
        )

    return (
        pd.DataFrame(comparison_rows),
        pd.concat(metrics_frames, ignore_index=True),
        pd.concat(prediction_frames, ignore_index=True),
        pd.DataFrame(feature_rows),
    )


def main():
    args = parse_args()
    windows = parse_csv(args.windows)
    feature_sets = parse_csv(args.feature_sets)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_prepared_dataset(args.matrix)
    balanced_df = balance_labels(df, args.random_state)
    balanced_df.to_csv(output_dir / "balanced_dataset.tsv", sep="\t", index=False)

    outputs = []
    outputs.append(run_panel(df, "full", feature_sets, windows, args, output_dir))
    outputs.append(run_panel(balanced_df, "balanced", feature_sets, windows, args, output_dir))

    comparison = pd.concat([item[0] for item in outputs], ignore_index=True)
    metrics = pd.concat([item[1] for item in outputs], ignore_index=True)
    predictions = pd.concat([item[2] for item in outputs], ignore_index=True)
    features = pd.concat([item[3] for item in outputs], ignore_index=True)
    deltas = summarize_deltas(metrics)

    comparison.to_csv(output_dir / "comparison.tsv", sep="\t", index=False)
    metrics.to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    predictions.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    features.to_csv(output_dir / "features.tsv", sep="\t", index=False)
    deltas.to_csv(output_dir / "feature_set_deltas.tsv", sep="\t", index=False)
    print(comparison)
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

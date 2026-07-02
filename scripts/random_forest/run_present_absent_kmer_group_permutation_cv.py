#!/usr/bin/env python3
"""Gene-grouped CV k-mer permutation importance for present/absent RF models."""

import argparse
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from introner_rf.config import ModelConfig
from introner_rf.data import load_prepared_dataset
from introner_rf.feature_selection import get_feature_cols
from introner_rf.modeling import train_model


KMER_RE = re.compile(r"^kmer_([ACGT]{3})_(left|right)_(\d+)$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--feature-set", default="kmer_left_right")
    parser.add_argument("--group-col", default="group_id")
    parser.add_argument("--panel", choices=("balanced", "full"), default="balanced")
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-repeats", type=int, default=100)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel permutation jobs. Training still uses this value, but prediction is forced to one thread inside each permutation task.",
    )
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


def balance_labels(df, random_state):
    counts = df["label"].value_counts()
    n = int(counts.min())
    return (
        df.groupby("label", group_keys=False)
        .sample(n=n, random_state=random_state)
        .sample(frac=1.0, random_state=random_state)
        .reset_index(drop=True)
    )


def kmer_group(feature):
    match = KMER_RE.match(feature)
    if not match:
        return None
    return match.group(1)


def build_feature_groups(feature_cols):
    groups = {}
    for col in feature_cols:
        group = kmer_group(col)
        if group is not None:
            groups.setdefault(group, []).append(col)
    if not groups:
        raise ValueError("No k-mer feature groups found")
    return groups


def cv_splits(df, group_col, n_splits, random_state):
    cv = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_state,
    )
    X_placeholder = np.zeros((len(df), 1))
    return list(cv.split(X_placeholder, df["label"], df[group_col]))


def write_config(args, output_dir, windows):
    lines = [
        "Within-sample present/absent k-mer permutation importance",
        "",
        "Purpose:",
        "Train kmer_left_right random forests under gene-grouped CV, then jointly permute each 3-mer group on held-out genes.",
        "",
        "Command:",
        " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "",
        "Configuration:",
        f"matrix: {args.matrix}",
        f"windows: {','.join(windows)}",
        f"feature_set: {args.feature_set}",
        f"group_col: {args.group_col}",
        f"panel: {args.panel}",
        f"n_splits: {args.n_splits}",
        f"n_repeats: {args.n_repeats}",
        f"n_estimators: {args.n_estimators}",
        f"n_jobs: {args.n_jobs}",
        "class_weight: balanced",
        f"random_state: {args.random_state}",
        "",
        "Importance definition:",
        "baseline ROC-AUC minus permuted ROC-AUC after jointly permuting all left/right/window features for each 3-mer group in held-out CV rows.",
    ]
    (output_dir / "config.txt").write_text("\n".join(lines) + "\n")


def summarize_importance(by_repeat):
    summary = (
        by_repeat
        .groupby("kmer", as_index=False)
        .agg(
            n_features=("n_features", "first"),
            features=("features", "first"),
            roc_auc_importance_mean=("roc_auc_importance", "mean"),
            roc_auc_importance_std=("roc_auc_importance", "std"),
            roc_auc_importance_sem=("roc_auc_importance", lambda s: s.std() / np.sqrt(len(s))),
            pr_auc_importance_mean=("pr_auc_importance", "mean"),
            pr_auc_importance_std=("pr_auc_importance", "std"),
            pr_auc_importance_sem=("pr_auc_importance", lambda s: s.std() / np.sqrt(len(s))),
            baseline_roc_auc_mean=("baseline_roc_auc", "mean"),
            permuted_roc_auc_mean=("permuted_roc_auc", "mean"),
            baseline_pr_auc_mean=("baseline_pr_auc", "mean"),
            permuted_pr_auc_mean=("permuted_pr_auc", "mean"),
            n_folds=("fold", "nunique"),
            n_repeats=("repeat", "nunique"),
        )
    )
    summary["roc_auc_importance_ci95"] = 1.96 * summary["roc_auc_importance_sem"]
    summary["pr_auc_importance_ci95"] = 1.96 * summary["pr_auc_importance_sem"]
    summary["roc_auc_importance_rank"] = (
        summary["roc_auc_importance_mean"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    summary["pr_auc_importance_rank"] = (
        summary["pr_auc_importance_mean"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    return summary.sort_values("roc_auc_importance_rank")


def permutation_group_rows(
    clf,
    X_test,
    y_test,
    baseline_roc_auc,
    baseline_pr_auc,
    panel,
    feature_set,
    fold,
    group,
    cols,
    n_repeats,
    random_state,
):
    rows = []
    for repeat in range(1, n_repeats + 1):
        rng = np.random.default_rng(
            random_state + fold * 1_000_000 + repeat * 1_000 + sum(map(ord, group))
        )
        permuted = X_test.copy()
        order = rng.permutation(len(permuted))
        permuted.loc[:, cols] = permuted[cols].to_numpy()[order, :]
        permuted_prob = clf.predict_proba(permuted)[:, 1]
        permuted_roc_auc = roc_auc_score(y_test, permuted_prob)
        permuted_pr_auc = average_precision_score(y_test, permuted_prob)
        rows.append({
            "panel": panel,
            "feature_set": feature_set,
            "fold": fold,
            "kmer": group,
            "repeat": repeat,
            "n_features": len(cols),
            "features": ",".join(cols),
            "baseline_roc_auc": baseline_roc_auc,
            "permuted_roc_auc": permuted_roc_auc,
            "roc_auc_importance": baseline_roc_auc - permuted_roc_auc,
            "baseline_pr_auc": baseline_pr_auc,
            "permuted_pr_auc": permuted_pr_auc,
            "pr_auc_importance": baseline_pr_auc - permuted_pr_auc,
        })
    return rows


def main():
    args = parse_args()
    windows = parse_csv(args.windows)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_config(args, output_dir, windows)

    df = load_prepared_dataset(args.matrix)
    if args.panel == "balanced":
        df = balance_labels(df, args.random_state)
        df.to_csv(output_dir / "balanced_dataset.tsv", sep="\t", index=False)

    feature_cols = get_feature_cols(df, windows, feature_set=args.feature_set)
    if not feature_cols:
        raise ValueError(f"No features selected for {args.feature_set!r}")
    feature_groups = build_feature_groups(feature_cols)
    print(f"Using {len(feature_cols)} features in {len(feature_groups)} k-mer groups")
    print(f"Using {args.n_jobs} parallel worker processes for permutation tasks", flush=True)

    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=args.n_jobs,
    )
    splits = cv_splits(df, args.group_col, args.n_splits, args.random_state)
    baseline_rows = []
    prediction_frames = []
    repeat_rows = []
    meta_cols = metadata_columns(df)

    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train_df = df.iloc[train_idx]
        test_df = df.iloc[test_idx]
        train_groups = set(train_df[args.group_col].astype(str))
        test_groups = set(test_df[args.group_col].astype(str))
        group_overlap = train_groups & test_groups
        if group_overlap:
            raise RuntimeError(
                f"fold {fold} has {len(group_overlap)} overlapping {args.group_col} groups"
            )

        clf = train_model(
            train_df[feature_cols],
            train_df["label"],
            model_config=model_config,
        )
        if args.n_jobs != 1:
            clf.set_params(n_jobs=1)
        X_test = test_df[feature_cols].copy()
        y_test = test_df["label"]
        baseline_prob = clf.predict_proba(X_test)[:, 1]
        baseline_roc_auc = roc_auc_score(y_test, baseline_prob)
        baseline_pr_auc = average_precision_score(y_test, baseline_prob)
        baseline_rows.append({
            "panel": args.panel,
            "feature_set": args.feature_set,
            "fold": fold,
            "n_train": len(train_idx),
            "n_train_pos": int(train_df["label"].eq(1).sum()),
            "n_train_neg": int(train_df["label"].eq(0).sum()),
            "n_test": len(test_idx),
            "n_test_pos": int(y_test.eq(1).sum()),
            "n_test_neg": int(y_test.eq(0).sum()),
            "n_train_groups": len(train_groups),
            "n_test_groups": len(test_groups),
            "n_group_overlap": len(group_overlap),
            "n_features": len(feature_cols),
            "n_kmer_groups": len(feature_groups),
            "baseline_roc_auc": baseline_roc_auc,
            "baseline_pr_auc": baseline_pr_auc,
        })

        predictions = test_df[meta_cols].copy()
        predictions["panel"] = args.panel
        predictions["feature_set"] = args.feature_set
        predictions["fold"] = fold
        predictions["prob_present"] = baseline_prob
        prediction_frames.append(predictions)

        print(f"fold {fold}: running {len(feature_groups)} k-mer permutation groups", flush=True)
        tasks = [
            delayed(permutation_group_rows)(
                clf,
                X_test,
                y_test,
                baseline_roc_auc,
                baseline_pr_auc,
                args.panel,
                args.feature_set,
                fold,
                group,
                cols,
                args.n_repeats,
                args.random_state,
            )
            for group, cols in sorted(feature_groups.items())
        ]
        group_rows = Parallel(n_jobs=args.n_jobs, backend="loky")(tasks)
        repeat_rows.extend(row for rows in group_rows for row in rows)

    baseline = pd.DataFrame(baseline_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    by_repeat = pd.DataFrame(repeat_rows)
    summary = summarize_importance(by_repeat)
    feature_groups_df = pd.DataFrame(
        [
            {"kmer": group, "feature": feature}
            for group, cols in feature_groups.items()
            for feature in cols
        ]
    )

    baseline.to_csv(output_dir / "baseline.tsv", sep="\t", index=False)
    predictions.to_csv(output_dir / "baseline_predictions.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "kmer_permutation_importance.tsv", sep="\t", index=False)
    by_repeat.to_csv(output_dir / "kmer_permutation_importance_by_fold_repeat.tsv", sep="\t", index=False)
    feature_groups_df.to_csv(output_dir / "feature_groups.tsv", sep="\t", index=False)
    print(summary.head(25).to_string(index=False))
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

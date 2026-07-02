#!/usr/bin/env python3
"""Grouped k-mer permutation for a model trained in one matrix and tested in another."""

import argparse
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score, roc_auc_score

from introner_rf.config import ModelConfig
from introner_rf.data import prepare_dataset
from introner_rf.feature_selection import get_feature_cols
from introner_rf.modeling import train_model


KMER_RE = re.compile(r"^kmer_([ACGT]{3})_(left|right)_(\d+)$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--train-matrix",
        default="rf_results/ccmp1545_present_absent_direct.mask3.w25_50_100.codon_microc_tpm.tsv",
    )
    parser.add_argument(
        "--test-matrix",
        default=(
            "rf_results/ccmp1545_train_rcc1749_preinsertion_matched_controls_mask3_orientnorm_exclude28/"
            "same_gene_gc.matched_dataset.tsv"
        ),
    )
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--feature-set", default="kmer_left_right")
    parser.add_argument(
        "--grouping",
        choices=("kmer", "kmer_side", "side_window", "window"),
        default="kmer",
    )
    parser.add_argument("--n-repeats", type=int, default=20)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel permutation jobs. Training still uses this value, but prediction is forced to one thread inside each permutation task.",
    )
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        default="rf_results/ccmp1545_train_rcc1749_preinsertion_same_gene_gc_kmer_group_permutation",
    )
    return parser.parse_args()


def parse_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def feature_group(feature, grouping):
    match = KMER_RE.match(feature)
    if not match:
        return None
    kmer, side, window = match.groups()
    if grouping == "kmer":
        return kmer
    if grouping == "kmer_side":
        return f"{kmer}_{side}"
    if grouping == "side_window":
        return f"{side}_{window}"
    if grouping == "window":
        return window
    raise ValueError(f"Unknown grouping: {grouping}")


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
        "label",
        "matched_pair_id",
        "match_case_sequence_id",
        "match_case_gene",
        "match_gc_col",
        "match_abs_gc_delta",
        "match_rank",
    ]
    return [col for col in preferred if col in df.columns]


def write_config(args, output_dir, windows):
    lines = [
        "Cross-sample grouped k-mer permutation importance",
        "",
        "Purpose:",
        "Train the CCMP1545 present/absent random forest with kmer_left_right features, then measure which 3-mer groups matter when scoring RCC1749 same-gene+GC matched empty-ortholog controls.",
        "",
        "Command:",
        " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "",
        "Configuration:",
        f"train_matrix: {args.train_matrix}",
        f"test_matrix: {args.test_matrix}",
        f"windows: {','.join(windows)}",
        f"feature_set: {args.feature_set}",
        f"grouping: {args.grouping}",
        f"n_repeats: {args.n_repeats}",
        f"n_estimators: {args.n_estimators}",
        f"n_jobs: {args.n_jobs}",
        "class_weight: balanced",
        f"random_state: {args.random_state}",
        "",
        "Importance definition:",
        "baseline ROC-AUC minus permuted ROC-AUC after jointly permuting every feature in a k-mer group across RCC1749 test rows.",
    ]
    (output_dir / "config.txt").write_text("\n".join(lines) + "\n")


def build_feature_groups(feature_cols, grouping):
    feature_groups = {}
    for col in feature_cols:
        group = feature_group(col, grouping)
        if group is not None:
            feature_groups.setdefault(group, []).append(col)
    if not feature_groups:
        raise ValueError(f"No feature groups found with grouping {grouping!r}")
    return feature_groups


def permutation_group_rows(
    clf,
    X_test,
    y_test,
    baseline_roc_auc,
    baseline_pr_auc,
    group,
    cols,
    n_repeats,
    random_state,
):
    rows = []
    for repeat in range(1, n_repeats + 1):
        rng = np.random.default_rng(random_state + repeat * 1_000 + sum(map(ord, group)))
        permuted = X_test.copy()
        order = rng.permutation(len(permuted))
        permuted.loc[:, cols] = permuted[cols].to_numpy()[order, :]
        permuted_prob = clf.predict_proba(permuted)[:, 1]
        permuted_roc_auc = roc_auc_score(y_test, permuted_prob)
        permuted_pr_auc = average_precision_score(y_test, permuted_prob)
        rows.append({
            "group": group,
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

    print("Loading matrices...")
    train_df = prepare_dataset(pd.read_csv(args.train_matrix, sep="\t", low_memory=False))
    test_df = prepare_dataset(pd.read_csv(args.test_matrix, sep="\t", low_memory=False))

    train_cols = get_feature_cols(train_df, windows, feature_set=args.feature_set)
    test_numeric = {
        col for col in test_df.columns
        if pd.api.types.is_numeric_dtype(test_df[col])
    }
    feature_cols = [col for col in train_cols if col in test_numeric]
    if not feature_cols:
        raise ValueError(f"No shared numeric features selected for {args.feature_set!r}")
    feature_groups = build_feature_groups(feature_cols, args.grouping)
    print(
        f"Using {len(feature_cols)} features in {len(feature_groups)} "
        f"{args.grouping} groups"
    )
    print(f"Using {args.n_jobs} parallel worker processes for permutation tasks", flush=True)

    model_config = ModelConfig(
        n_estimators=args.n_estimators,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=args.n_jobs,
    )
    clf = train_model(train_df[feature_cols], train_df["label"], model_config=model_config)
    if args.n_jobs != 1:
        clf.set_params(n_jobs=1)

    X_test = test_df[feature_cols].copy()
    y_test = test_df["label"]
    baseline_prob = clf.predict_proba(X_test)[:, 1]
    baseline_roc_auc = roc_auc_score(y_test, baseline_prob)
    baseline_pr_auc = average_precision_score(y_test, baseline_prob)
    pd.DataFrame([{
        "feature_set": args.feature_set,
        "grouping": args.grouping,
        "n_train": len(train_df),
        "n_train_pos": int(train_df["label"].eq(1).sum()),
        "n_train_neg": int(train_df["label"].eq(0).sum()),
        "n_test": len(test_df),
        "n_test_pos": int(test_df["label"].eq(1).sum()),
        "n_test_neg": int(test_df["label"].eq(0).sum()),
        "n_features": len(feature_cols),
        "n_feature_groups": len(feature_groups),
        "baseline_roc_auc": baseline_roc_auc,
        "baseline_pr_auc": baseline_pr_auc,
    }]).to_csv(output_dir / "baseline.tsv", sep="\t", index=False)

    predictions = test_df[metadata_columns(test_df)].copy()
    predictions["prob_present"] = baseline_prob
    predictions.to_csv(output_dir / "baseline_predictions.tsv", sep="\t", index=False)

    tasks = [
        delayed(permutation_group_rows)(
            clf,
            X_test,
            y_test,
            baseline_roc_auc,
            baseline_pr_auc,
            group,
            cols,
            args.n_repeats,
            args.random_state,
        )
        for group, cols in sorted(feature_groups.items())
    ]
    group_rows = Parallel(n_jobs=args.n_jobs, backend="loky")(tasks)
    rows = [row for group_result in group_rows for row in group_result]

    by_repeat = pd.DataFrame(rows)
    summary = (
        by_repeat
        .groupby("group", as_index=False)
        .agg(
            n_features=("n_features", "first"),
            features=("features", "first"),
            roc_auc_importance_mean=("roc_auc_importance", "mean"),
            roc_auc_importance_std=("roc_auc_importance", "std"),
            roc_auc_importance_sem=("roc_auc_importance", lambda s: s.std() / np.sqrt(len(s))),
            pr_auc_importance_mean=("pr_auc_importance", "mean"),
            pr_auc_importance_std=("pr_auc_importance", "std"),
            baseline_roc_auc_mean=("baseline_roc_auc", "mean"),
            permuted_roc_auc_mean=("permuted_roc_auc", "mean"),
            baseline_pr_auc_mean=("baseline_pr_auc", "mean"),
            permuted_pr_auc_mean=("permuted_pr_auc", "mean"),
            n_repeats=("repeat", "nunique"),
        )
    )
    summary["roc_auc_importance_ci95"] = 1.96 * summary["roc_auc_importance_sem"]
    summary = summary.sort_values("roc_auc_importance_mean", ascending=False)

    summary.to_csv(output_dir / f"{args.grouping}_permutation_importance.tsv", sep="\t", index=False)
    by_repeat.to_csv(output_dir / f"{args.grouping}_permutation_importance_by_repeat.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"group": group, "feature": feature} for group, cols in feature_groups.items() for feature in cols]
    ).to_csv(output_dir / "feature_groups.tsv", sep="\t", index=False)

    print(summary.head(25).to_string(index=False))
    print(f"Wrote outputs under {output_dir}")


if __name__ == "__main__":
    main()

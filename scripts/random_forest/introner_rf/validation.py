from numpy.typing import NDArray
import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold, train_test_split

from introner_rf.config import ModelConfig
from introner_rf.modeling import train_model


def split_dataset(df, feature_cols, label_col="label"):
    X = df[feature_cols]
    y = df[label_col]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.2,
        stratify=y,
        random_state=42,
    )
    return X_train, X_test, y_train, y_test


def evaluate_model(model, X_test: NDArray, y_test: NDArray):
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    metrics = {
        "accuracy": accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred),
        "recall": recall_score(y_test, y_pred),
        "f1": f1_score(y_test, y_pred),
        "roc_auc": roc_auc_score(y_test, y_prob),
        "pr_auc": average_precision_score(y_test, y_prob),
    }

    cm = confusion_matrix(y_test, y_pred)
    print(classification_report(y_test, y_pred))
    print(cm)

    return metrics, cm


def eval_permutation_importance(clf, X_test, y_test, feature_cols):
    result = permutation_importance(
        clf,
        X_test,
        y_test,
        n_repeats=10,
        random_state=42,
        n_jobs=-1,
        scoring="roc_auc",
    )

    imp = pd.DataFrame({
        "feature": feature_cols,
        "importance_mean": result.importances_mean,
        "importance_std": result.importances_std,
    }).sort_values("importance_mean", ascending=False)
    return imp


def eval_grouped_permutation_importance(
    df,
    feature_cols,
    group_col,
    label_col="label",
    n_splits=5,
    model_config: ModelConfig | None = None,
    n_repeats=10,
    scoring="roc_auc",
    random_state=42,
):
    X = df[feature_cols]
    y = df[label_col]
    groups = df[group_col]

    cv = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_state,
    )

    fold_rows = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        clf = train_model(X_train, y_train, model_config=model_config)
        result = permutation_importance(
            clf,
            X_test,
            y_test,
            n_repeats=n_repeats,
            random_state=random_state + fold,
            n_jobs=-1,
            scoring=scoring,
        )

        fold_rows.append(pd.DataFrame({
            "fold": fold,
            "feature": feature_cols,
            "importance_mean": result.importances_mean,
            "importance_std": result.importances_std,
            "n_test": len(test_idx),
            "n_test_pos": int(y_test.sum()),
            "n_test_neg": int((y_test == 0).sum()),
        }))

    fold_importance = pd.concat(fold_rows, ignore_index=True)
    summary = (
        fold_importance
        .groupby("feature", as_index=False)
        .agg(
            importance_mean=("importance_mean", "mean"),
            importance_std_across_folds=("importance_mean", "std"),
            importance_sem=("importance_mean", lambda s: s.std() / np.sqrt(len(s))),
            mean_within_fold_std=("importance_std", "mean"),
            n_folds=("fold", "nunique"),
        )
    )
    summary["importance_ci95"] = 1.96 * summary["importance_sem"]
    summary = summary.sort_values("importance_mean", ascending=False)
    return summary, fold_importance


def evaluate_predictions(y_true, y_pred, y_prob):
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred),
        "recall": recall_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_prob),
        "pr_auc": average_precision_score(y_true, y_prob),
    }


def run_grouped_cv(
    df,
    feature_cols,
    group_col,
    label_col="label",
    n_splits=5,
    model_config: ModelConfig | None = None,
):
    X = df[feature_cols]
    y = df[label_col]
    groups = df[group_col]

    cv = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=42,
    )

    rows = []
    prediction_rows = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        prediction_cols = [
            "sequence_id",
            "gene",
            "contig",
            "start",
            "end",
            "label",
            "group_id",
            "matched_pair_id",
            "match_case_sequence_id",
            "match_gc_delta",
            "match_abs_gc_delta",
            "match_rank",
        ]
        prediction_cols = [c for c in prediction_cols if c in df.columns]
        pred_rows = df.iloc[test_idx][prediction_cols].copy()

        clf = train_model(X_train, y_train, model_config=model_config)

        y_pred = clf.predict(X_test)
        y_prob = clf.predict_proba(X_test)[:, 1]

        metrics = evaluate_predictions(y_test, y_pred, y_prob)
        metrics["fold"] = fold
        metrics["n_train"] = len(train_idx)
        metrics["n_test"] = len(test_idx)
        metrics["n_test_pos"] = int(y_test.sum())
        metrics["n_test_neg"] = int((y_test == 0).sum())
        rows.append(metrics)

        pred_rows["fold"] = fold
        pred_rows["y_pred"] = y_pred
        pred_rows["prob_introner"] = y_prob
        pred_rows[label_col] = y_test.values
        pred_rows["correct"] = pred_rows[label_col] == pred_rows["y_pred"]
        prediction_rows.append(pred_rows)

    metrics_df = pd.DataFrame(rows)
    predictions_df = pd.concat(prediction_rows, ignore_index=True)

    return metrics_df, predictions_df


def summarize_cv(metrics_df):
    metric_cols = [
        "accuracy",
        "balanced_accuracy",
        "precision",
        "recall",
        "f1",
        "roc_auc",
        "pr_auc",
    ]

    summary = metrics_df[metric_cols].agg(["mean", "std"]).T
    summary["sem"] = summary["std"] / np.sqrt(len(metrics_df))
    summary["ci95"] = 1.96 * summary["sem"]
    return summary

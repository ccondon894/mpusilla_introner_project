from pathlib import Path

import pandas as pd

from introner_rf.config import FeatureConfig, MatchConfig, ModelConfig
from introner_rf.data import load_prepared_dataset
from introner_rf.feature_selection import DEFAULT_FEATURE_SET, get_feature_cols
from introner_rf.matching import build_matched_control_dataset, subset_family_vs_introns
from introner_rf.modeling import train_model
from introner_rf.validation import (
    eval_permutation_importance,
    eval_grouped_permutation_importance,
    evaluate_model,
    run_grouped_cv,
    split_dataset,
    summarize_cv,
)


def output_path(output_dir, filename):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / filename


def prefixed_output_path(output_dir, output_prefix, suffix):
    output_prefix = Path(output_prefix)
    if output_prefix.is_absolute() or output_prefix.parent != Path("."):
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        return output_prefix.with_name(f"{output_prefix.name}.{suffix}")
    return output_path(output_dir, f"{output_prefix.name}.{suffix}")


def run_matched_control_cv(
    df,
    feature_config: FeatureConfig,
    match_config: MatchConfig | None = None,
    model_config: ModelConfig | None = None,
):
    match_config = match_config or MatchConfig()
    matched_df, stats = build_matched_control_dataset(
        df,
        gc_col=match_config.gc_col,
        controls_per_positive=match_config.controls_per_positive,
        gc_tolerance=match_config.gc_tolerance,
        same_contig=match_config.same_contig,
        same_gene=match_config.same_gene,
        drop_unannotated_positives=match_config.drop_unannotated_positives,
        replace_controls=match_config.replace_controls,
        random_state=match_config.random_state,
    )

    print("Matched-control dataset:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    matched_df.to_csv(
        prefixed_output_path(
            match_config.output_dir,
            match_config.output_prefix,
            "dataset.tsv",
        ),
        sep="\t",
        index=False,
    )
    pd.DataFrame([stats]).to_csv(
        prefixed_output_path(
            match_config.output_dir,
            match_config.output_prefix,
            "match_summary.tsv",
        ),
        sep="\t",
        index=False,
    )

    feature_cols = get_feature_cols(
        matched_df,
        feature_config.windows,
        feature_set=feature_config.feature_set,
    )
    print(f"Using feature set {feature_config.feature_set!r}: {len(feature_cols)} features")
    metrics_df, predictions_df = run_grouped_cv(
        matched_df,
        feature_cols,
        group_col=match_config.group_col,
        n_splits=5,
        model_config=model_config,
    )
    metrics_df.to_csv(
        prefixed_output_path(
            match_config.output_dir,
            match_config.output_prefix,
            "grouped_cv_metrics.tsv",
        ),
        sep="\t",
        index=False,
    )
    predictions_df.to_csv(
        prefixed_output_path(
            match_config.output_dir,
            match_config.output_prefix,
            "grouped_cv_predictions.tsv",
        ),
        sep="\t",
        index=False,
    )

    summary = summarize_cv(metrics_df)
    summary.to_csv(
        prefixed_output_path(
            match_config.output_dir,
            match_config.output_prefix,
            "grouped_cv_summary.tsv",
        ),
        sep="\t",
    )
    print(summary)
    return matched_df, metrics_df, predictions_df, summary


def run_random_forest(
    matrix,
    windows,
    feature_set=DEFAULT_FEATURE_SET,
    permutation=False,
    run_cv=False,
    run_matched_controls=False,
    match_gc_window=250,
    match_controls_per_positive=1,
    match_gc_tolerance=None,
    match_across_contigs=False,
    match_same_gene=False,
    match_keep_unannotated_positives=False,
    match_replace_controls=False,
    matched_group_col="matched_pair_id",
    matched_output_prefix="matched_controls",
    output_dir="rf_results",
    matched_permutation=False,
    permutation_repeats=10,
    model_config: ModelConfig | None = None,
):
    feature_config = FeatureConfig(windows=windows.split(","), feature_set=feature_set)
    model_config = model_config or ModelConfig()

    print("loading dataset...")
    df = load_prepared_dataset(matrix)
    feature_cols = get_feature_cols(
        df,
        feature_config.windows,
        feature_set=feature_config.feature_set,
    )
    print(f"Using feature set {feature_config.feature_set!r}: {len(feature_cols)} features")
    X_train, X_test, y_train, y_test = split_dataset(df, feature_cols)

    print("training model...")
    clf = train_model(X_train, y_train, model_config=model_config)

    print("evaluating model...")
    metrics, cm = evaluate_model(clf, X_test, y_test)

    print(metrics)
    print(cm)

    if permutation:
        print("running permutation importance...")
        imp = eval_permutation_importance(clf, X_test, y_test, feature_cols)
        imp.to_csv(output_path(output_dir, "importances.tsv"), sep="\t", index=False)

    if run_cv:
        for family in sorted(df.loc[df["label"] == 1, "family"].dropna().unique()):
            print("Running CV on family:", family)
            family_df = subset_family_vs_introns(
                df,
                family,
                neg_pos_ratio=3,
                random_state=42,
            )
            family_feature_cols = get_feature_cols(
                family_df,
                feature_config.windows,
                feature_set=feature_config.feature_set,
            )
            metrics_df, predictions_df = run_grouped_cv(
                family_df,
                family_feature_cols,
                group_col="group_id",
                n_splits=5,
                model_config=model_config,
            )
            metrics_df.to_csv(
                output_path(output_dir, f"grouped_cv_metrics_family_{family}.tsv"),
                sep="\t",
                index=False,
            )
            predictions_df.to_csv(
                output_path(output_dir, f"grouped_cv_predictions_family_{family}.tsv"),
                sep="\t",
                index=False,
            )

            summary = summarize_cv(metrics_df)
            summary.to_csv(
                output_path(output_dir, f"grouped_cv_summary_family_{family}.tsv"),
                sep="\t",
            )

    if run_matched_controls:
        print("running matched-control grouped CV...")
        match_config = MatchConfig(
            gc_window=match_gc_window,
            controls_per_positive=match_controls_per_positive,
            gc_tolerance=match_gc_tolerance,
            same_contig=not match_across_contigs,
            same_gene=match_same_gene,
            drop_unannotated_positives=not match_keep_unannotated_positives,
            replace_controls=match_replace_controls,
            group_col=matched_group_col,
            output_prefix=matched_output_prefix,
            output_dir=output_dir,
            random_state=42,
        )
        matched_df, _, _, _ = run_matched_control_cv(
            df,
            feature_config=feature_config,
            match_config=match_config,
            model_config=model_config,
        )
        if matched_permutation:
            print("running matched-control grouped permutation importance...")
            matched_feature_cols = get_feature_cols(
                matched_df,
                feature_config.windows,
                feature_set=feature_config.feature_set,
            )
            imp_summary, imp_folds = eval_grouped_permutation_importance(
                matched_df,
                matched_feature_cols,
                group_col=match_config.group_col,
                n_splits=5,
                model_config=model_config,
                n_repeats=permutation_repeats,
            )
            imp_summary.to_csv(
                prefixed_output_path(
                    match_config.output_dir,
                    match_config.output_prefix,
                    "grouped_permutation_importance.tsv",
                ),
                sep="\t",
                index=False,
            )
            imp_folds.to_csv(
                prefixed_output_path(
                    match_config.output_dir,
                    match_config.output_prefix,
                    "grouped_permutation_importance_by_fold.tsv",
                ),
                sep="\t",
                index=False,
            )
            print(imp_summary.head(25).to_string(index=False))

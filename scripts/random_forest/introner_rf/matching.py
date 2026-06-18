import numpy as np
import pandas as pd


def subset_family_vs_introns(df, family, neg_pos_ratio=3, random_state=42):
    family_df = df[(df["label"] == 1) & (df["family"] == family)].copy()
    intron_df = df[(df["label"] == 0)].copy()

    n_pos = len(family_df)
    n_neg = min(len(intron_df), n_pos * neg_pos_ratio)

    intron_df = intron_df.sample(
        n=n_neg,
        random_state=random_state,
        replace=False
    )

    out = pd.concat([family_df, intron_df], ignore_index=True)
    out["family_label"] = out["label"]
    return out


def build_matched_control_dataset(
    df,
    gc_col="gc_combined_250",
    controls_per_positive=1,
    gc_tolerance=None,
    same_contig=True,
    same_gene=False,
    drop_unannotated_positives=True,
    replace_controls=False,
    random_state=42,
):
    """Build an introner/control set matched by genomic context and flank GC."""
    if gc_col not in df.columns:
        raise ValueError(f"Cannot match controls: missing column {gc_col!r}")

    if same_gene and "gene" not in df.columns:
        raise ValueError("Cannot match controls within gene: missing column 'gene'")

    if controls_per_positive < 1:
        raise ValueError("controls_per_positive must be at least 1")

    positives = df.loc[df["label"] == 1].copy()
    controls = df.loc[df["label"] == 0].copy()

    n_positive_initial = len(positives)
    if drop_unannotated_positives or same_gene:
        positives = positives.loc[positives["gene"].notna()].copy()

    positives = positives.loc[positives[gc_col].notna()].copy()
    controls = controls.loc[controls[gc_col].notna()].copy()
    if same_gene:
        controls = controls.loc[controls["gene"].notna()].copy()
    controls["_control_row_id"] = np.arange(len(controls))

    positives = positives.sample(frac=1.0, random_state=random_state)

    used_control_ids = set()
    matched_rows = []
    n_skipped = 0
    pair_idx = 0

    for _, case in positives.iterrows():
        case_gc = case[gc_col]
        candidates = controls

        if same_contig:
            candidates = candidates.loc[candidates["contig"] == case["contig"]]

        if same_gene:
            candidates = candidates.loc[candidates["gene"] == case["gene"]]

        if not replace_controls:
            candidates = candidates.loc[
                ~candidates["_control_row_id"].isin(used_control_ids)
            ]

        if candidates.empty:
            n_skipped += 1
            continue

        candidates = candidates.copy()
        candidates["_match_gc_delta"] = candidates[gc_col] - case_gc
        candidates["_match_abs_gc_delta"] = candidates["_match_gc_delta"].abs()

        if gc_tolerance is not None:
            candidates = candidates.loc[
                candidates["_match_abs_gc_delta"] <= gc_tolerance
            ]

        if len(candidates) < controls_per_positive:
            n_skipped += 1
            continue

        chosen = candidates.nsmallest(
            controls_per_positive,
            ["_match_abs_gc_delta", "_control_row_id"],
        )

        pair_idx += 1
        pair_id = f"match_{pair_idx:06d}"
        case_dict = case.to_dict()
        case_dict.update({
            "matched_pair_id": pair_id,
            "match_case_sequence_id": case["sequence_id"],
            "match_case_gene": case["gene"],
            "match_case_contig": case["contig"],
            "match_case_start": case["start"],
            "match_case_end": case["end"],
            "match_gc_col": gc_col,
            "match_gc_delta": 0.0,
            "match_abs_gc_delta": 0.0,
            "match_rank": 0,
        })
        matched_rows.append(case_dict)

        for rank, (_, control) in enumerate(chosen.iterrows(), start=1):
            control_dict = control.to_dict()
            control_id = control_dict.pop("_control_row_id")
            gc_delta = control_dict.pop("_match_gc_delta")
            abs_gc_delta = control_dict.pop("_match_abs_gc_delta")

            control_dict.update({
                "matched_pair_id": pair_id,
                "match_case_sequence_id": case["sequence_id"],
                "match_case_gene": case["gene"],
                "match_case_contig": case["contig"],
                "match_case_start": case["start"],
                "match_case_end": case["end"],
                "match_gc_col": gc_col,
                "match_gc_delta": gc_delta,
                "match_abs_gc_delta": abs_gc_delta,
                "match_rank": rank,
            })
            matched_rows.append(control_dict)
            used_control_ids.add(control_id)

    if not matched_rows:
        raise ValueError("No matched controls were found with the requested settings")

    matched_df = pd.DataFrame(matched_rows)
    matched_df = matched_df.drop(columns=["_control_row_id"], errors="ignore")

    stats = {
        "n_positive_initial": n_positive_initial,
        "n_positive_eligible": len(positives),
        "n_positive_matched": pair_idx,
        "n_positive_skipped": n_skipped,
        "n_controls_available": len(controls),
        "n_controls_used": int((matched_df["label"] == 0).sum()),
        "controls_per_positive": controls_per_positive,
        "gc_col": gc_col,
        "gc_tolerance": gc_tolerance,
        "same_contig": same_contig,
        "same_gene": same_gene,
        "drop_unannotated_positives": drop_unannotated_positives,
        "replace_controls": replace_controls,
        "mean_abs_gc_delta": matched_df.loc[
            matched_df["label"] == 0, "match_abs_gc_delta"
        ].mean(),
        "median_abs_gc_delta": matched_df.loc[
            matched_df["label"] == 0, "match_abs_gc_delta"
        ].median(),
        "max_abs_gc_delta": matched_df.loc[
            matched_df["label"] == 0, "match_abs_gc_delta"
        ].max(),
    }

    return matched_df, stats

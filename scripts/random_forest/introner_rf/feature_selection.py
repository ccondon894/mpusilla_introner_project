import pandas as pd


DEFAULT_FEATURE_SET = "kmer_left_right"

FEATURE_SET_CHOICES = (
    "all",
    "kmer",
    "kmer_left_right",
    "kmer_delta",
    "scalar",
    "composition",
    "codon",
    "microc",
    "all_microc",
    "kmer_microc",
    "kmer_left_right_microc",
    "kmer_left_right_codon",
    "kmer_left_right_microc_codon",
    "recomb",
    "kmer_left_right_recomb",
    "kmer_left_right_microc_recomb",
    "kmer_left_right_microc_codon_recomb",
    "expression",
    "microc_expression",
    "kmer_left_right_expression",
    "kmer_left_right_microc_expression",
)

METADATA_COLS = {
    "sequence_id",
    "ortholog_id",
    "sample",
    "source",
    "gene",
    "gene_from_genotype_matrix",
    "gtf_gene_source",
    "gtf_gene_overlap_bp",
    "gtf_gene_n_overlaps",
    "transcript_id",
    "strand",
    "group_id",
    "family",
    "family_label",
    "class_source",
    "contig",
    "start",
    "end",
    "site",
    "feature_start",
    "feature_end",
    "absent_span",
    "exon_start",
    "exon_end",
    "group1_introner_status",
    "group1_absent_vs_fixed_status",
    "group1_pattern",
    "group1_present_count",
    "group1_absent_count",
    "group1_missing_count",
    "group1_callable_count",
    "group1_n_samples",
    "within_group_status",
    "within_group_identity",
    "label",
    "boundary_support_ortholog_id",
    "boundary_support_best_junction_name",
    "boundary_support_best_junction_score",
    "boundary_support_start_delta",
    "boundary_support_end_delta",
    "boundary_support_max_abs_boundary_delta",
    "boundary_support_sum_abs_boundary_delta",
    "boundary_support_applied",
    "matched_pair_id",
    "match_case_sequence_id",
    "match_case_gene",
    "match_case_contig",
    "match_case_start",
    "match_case_end",
    "match_gc_col",
    "match_gc_delta",
    "match_abs_gc_delta",
    "match_rank",
    "microc_midpoint",
}

PYRHO_DIAGNOSTIC_COLS = {
    "pyrho_contig_available",
    "pyrho_any_overlap",
}

CODON_DIAGNOSTIC_COLS = {
    "codon_match_found",
    "codon_cds_intron_start",
    "codon_cds_intron_end",
    "codon_start_delta",
    "codon_end_delta",
    "codon_abs_start_delta",
    "codon_abs_end_delta",
    "codon_boundary_distance_sum",
    "codon_has_left_cds_boundary",
    "codon_has_right_cds_boundary",
    "codon_has_both_cds_boundaries",
}


def parse_feature_col(col):
    rest, window = col.rsplit("_", 1)
    type_part, view = rest.rsplit("_", 1)
    feature_type = "kmer" if type_part.startswith("kmer_") else type_part
    return feature_type, view, int(window)


def feature_set_filter(feature_set):
    feature_set = (
        feature_set
        .replace("_microc", "")
        .replace("_codon", "")
        .replace("_recomb", "")
        .replace("_expression", "")
    )
    if feature_set == "all":
        return lambda _col, _feature_type, _view, _window: True
    if feature_set == "kmer":
        return lambda _col, feature_type, _view, _window: feature_type == "kmer"
    if feature_set == "kmer_left_right":
        return (
            lambda _col, feature_type, view, _window:
            feature_type == "kmer" and view in {"left", "right"}
        )
    if feature_set == "kmer_delta":
        return (
            lambda _col, feature_type, view, _window:
            feature_type == "kmer" and view == "delta"
        )
    if feature_set == "scalar":
        return lambda _col, feature_type, _view, _window: feature_type != "kmer"
    if feature_set == "composition":
        composition_types = {"gc", "cpg_density", "gt_density", "ag_density"}
        return lambda _col, feature_type, _view, _window: feature_type in composition_types
    raise ValueError(f"Unknown feature_set: {feature_set}")


def microc_feature_cols(df):
    return [
        c for c in df.columns
        if c not in METADATA_COLS
        and c.startswith("microc_")
        and not c.endswith("_domain_name")
        and pd.api.types.is_numeric_dtype(df[c])
    ]


def codon_feature_cols(df):
    return [
        c for c in df.columns
        if c not in METADATA_COLS
        and c not in CODON_DIAGNOSTIC_COLS
        and c.startswith("codon_")
        and pd.api.types.is_numeric_dtype(df[c])
    ]


def recomb_feature_cols(df):
    return [
        c for c in df.columns
        if c not in METADATA_COLS
        and c not in PYRHO_DIAGNOSTIC_COLS
        and c.startswith("pyrho_")
        and pd.api.types.is_numeric_dtype(df[c])
    ]


def expression_feature_cols(df):
    return [
        c for c in df.columns
        if c not in METADATA_COLS
        and c.startswith("expr_")
        and pd.api.types.is_numeric_dtype(df[c])
    ]


def get_feature_cols(df, windows: list[int | str] | None, feature_set=DEFAULT_FEATURE_SET) -> list:
    candidate_cols = [c for c in df.columns if c not in METADATA_COLS]
    microc_cols = microc_feature_cols(df)
    codon_cols = codon_feature_cols(df)
    recomb_cols = recomb_feature_cols(df)
    expression_cols = expression_feature_cols(df)
    sequence_candidates = [
        c for c in candidate_cols
        if (
            not c.startswith("microc_")
            and not c.startswith("codon_")
            and not c.startswith("pyrho_")
            and not c.startswith("expr_")
        )
    ]

    if feature_set == "microc":
        return microc_cols
    if feature_set == "codon":
        return codon_cols
    if feature_set == "recomb":
        return recomb_cols
    if feature_set == "expression":
        return expression_cols
    if feature_set == "microc_expression":
        return microc_cols + expression_cols

    selector = feature_set_filter(feature_set)
    selected_cols = []
    for col in sequence_candidates:
        try:
            feature_type, view, window = parse_feature_col(col)
        except ValueError:
            continue
        if selector(col, feature_type, view, window):
            selected_cols.append(col)

    if not windows:
        sequence_cols = selected_cols
    else:
        suffixes = tuple(f"_{w}" for w in windows)
        sequence_cols = [c for c in selected_cols if c.endswith(suffixes)]

    out_cols = list(sequence_cols)
    if "_microc" in feature_set:
        out_cols += microc_cols
    if "_codon" in feature_set:
        out_cols += codon_cols
    if "_recomb" in feature_set:
        out_cols += recomb_cols
    if "_expression" in feature_set:
        out_cols += expression_cols
    return out_cols

#!/usr/bin/env python3
"""Test whether introner-containing genes have elevated NMD isoform rates.

The unit of observation is a common gene in one strain.  NMD-positive and
NMD-negative multi-exon isoforms form a binomial response, so the analysis
tests the *fraction* of NMD-prone isoforms rather than the raw number of NMD
isoforms.  Population-averaged binomial GEEs account for repeated observations
of the same common gene across strains.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
SQANTI_CATEGORIES = {
    "full-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog",
}
PRIMARY_FORMULA_SUFFIX = (
    "C(strain) + log1p_long_read_count + log_gene_span + "
    "background_exon_count + log1p_n_isoforms"
)
RAW_EXON_FORMULA_SUFFIX = (
    "C(strain) + log1p_long_read_count + log_gene_span + "
    "max_annotated_exons + log1p_n_isoforms"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    parser.add_argument("--minimum-gene-read-support", type=int, default=10)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def extract_mating_type_genes(
    gtf_path: Path,
    scaffold: str = "CCMP1545#0#scaffold_2",
    start: int = 49808,
    end: int = 1730591,
) -> set[str]:
    genes: set[str] = set()
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[0] != scaffold:
                continue
            if int(fields[3]) > end or int(fields[4]) < start:
                continue
            match = re.search(r'gene_id[=\s]+"?([^;"]+)', fields[8])
            if match:
                genes.add(match.group(1))
    return genes


def coerce_nmd(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    mapped = series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    return mapped


def prepare_sqanti_isoforms(sample: str, path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    frame["exons"] = pd.to_numeric(frame["exons"], errors="coerce")
    frame["predicted_NMD"] = coerce_nmd(frame["predicted_NMD"])
    frame = frame[
        frame["structural_category"].isin(SQANTI_CATEGORIES)
        & (frame["exons"] > 1)
        & frame["associated_gene"].notna()
        & frame["predicted_NMD"].isin([True, False])
    ].copy()
    frame = frame.drop_duplicates("isoform")
    count_columns = [
        column for column in frame.columns if re.fullmatch(r"\d+_[A-Z]\d+", column)
    ]
    if not count_columns:
        raise ValueError(f"No R2C2 count columns found in {path}")
    for column in count_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0)
    frame["isoform_read_count"] = frame[count_columns].sum(axis=1)
    frame["strain"] = sample
    frame["common_gene_id"] = frame["associated_gene"].astype(str)
    return frame[
        [
            "common_gene_id",
            "strain",
            "isoform",
            "predicted_NMD",
            "isoform_read_count",
        ]
    ]


def aggregate_nmd(isoforms: pd.DataFrame, minimum_isoform_reads: int = 0) -> pd.DataFrame:
    frame = isoforms[isoforms["isoform_read_count"] >= minimum_isoform_reads].copy()
    frame["nmd_positive_isoform"] = frame["predicted_NMD"].astype(int)
    frame["nmd_negative_isoform"] = 1 - frame["nmd_positive_isoform"]
    frame["nmd_positive_reads"] = (
        frame["isoform_read_count"] * frame["nmd_positive_isoform"]
    )
    frame["nmd_negative_reads"] = (
        frame["isoform_read_count"] * frame["nmd_negative_isoform"]
    )
    result = (
        frame.groupby(["common_gene_id", "strain"], sort=False)
        .agg(
            n_nmd_positive=("nmd_positive_isoform", "sum"),
            n_nmd_negative=("nmd_negative_isoform", "sum"),
            n_nmd_positive_reads=("nmd_positive_reads", "sum"),
            n_nmd_negative_reads=("nmd_negative_reads", "sum"),
        )
        .reset_index()
    )
    result["n_nmd_evaluable"] = result["n_nmd_positive"] + result["n_nmd_negative"]
    result["nmd_isoform_fraction"] = (
        result["n_nmd_positive"] / result["n_nmd_evaluable"]
    )
    result["n_nmd_evaluable_reads"] = (
        result["n_nmd_positive_reads"] + result["n_nmd_negative_reads"]
    )
    result["nmd_read_fraction"] = np.where(
        result["n_nmd_evaluable_reads"] > 0,
        result["n_nmd_positive_reads"] / result["n_nmd_evaluable_reads"],
        np.nan,
    )
    return result


def direct_introner_counts(turnover_dir: Path) -> pd.DataFrame:
    loci = pd.read_csv(turnover_dir / "locus_gene_map.tsv", sep="\t", low_memory=False)
    crosswalk = pd.read_csv(turnover_dir / "sample_gene_crosswalk.tsv", sep="\t")
    crosswalk = crosswalk.dropna(subset=["common_gene_id"]).drop_duplicates(
        ["sample", "native_gene_id"]
    )
    loci = loci[(loci["presence"] == 1) & loci["native_gene_id"].notna()].copy()
    # The locus map also contains an ortholog-consensus common-gene assignment.
    # Drop it here because the direct definition intentionally uses each
    # strain's native gene followed by the exact SQANTI/GTF crosswalk.
    loci = loci.drop(columns="common_gene_id", errors="ignore")
    loci = loci.merge(
        crosswalk[["sample", "native_gene_id", "common_gene_id"]],
        on=["sample", "native_gene_id"],
        how="inner",
    )
    return (
        loci.groupby(["common_gene_id", "sample"], sort=False)["ortholog_id"]
        .nunique()
        .rename("direct_current_introner_count")
        .reset_index()
        .rename(columns={"sample": "strain"})
    )


def build_model_data(root: Path, minimum_isoform_reads: int = 0) -> pd.DataFrame:
    turnover_dir = root / "results" / "expression" / "isoform_analysis" / "turnover"
    architecture = pd.read_csv(
        turnover_dir / "isoform_introner_model_data.tsv", sep="\t", low_memory=False
    )
    sqanti_paths = {
        "CCMP1545": root
        / "data"
        / "sqanti3_output_834"
        / "834_isoforms_classification.filtered.txt",
        "RCC1614": root
        / "data"
        / "sqanti3_output_1614"
        / "1614_isoforms_classification.filtered.txt",
        "RCC1749": root
        / "data"
        / "sqanti3_output_1749"
        / "1749_isoforms_classification.filtered.txt",
    }
    isoforms = pd.concat(
        [prepare_sqanti_isoforms(sample, path) for sample, path in sqanti_paths.items()],
        ignore_index=True,
    )
    mt_genes = extract_mating_type_genes(
        root / "results" / "annotations" / "CCMP1545.gtf"
    )
    isoforms = isoforms[~isoforms["common_gene_id"].isin(mt_genes)].copy()
    nmd = aggregate_nmd(isoforms, minimum_isoform_reads=minimum_isoform_reads)

    keep = [
        "common_gene_id",
        "strain",
        "n_isoforms",
        "total_long_read_count",
        "gene_span",
        "max_annotated_exons",
        "current_introner_count",
        "current_missing_locus_count",
    ]
    result = nmd.merge(architecture[keep], on=["common_gene_id", "strain"], how="inner")
    result = result.merge(
        direct_introner_counts(turnover_dir),
        on=["common_gene_id", "strain"],
        how="left",
    )
    result["direct_current_introner_count"] = (
        result["direct_current_introner_count"].fillna(0).astype(int)
    )
    result["has_introner"] = (result["direct_current_introner_count"] > 0).astype(int)
    result["clean_has_introner"] = (result["current_introner_count"] > 0).astype(int)
    result["additional_introner_count"] = np.maximum(
        result["direct_current_introner_count"] - 1, 0
    )
    # Each spliced introner can itself split an annotated exon. Subtracting
    # current introners better approximates background gene architecture than
    # controlling for the raw exon count, which partly contains the exposure.
    result["background_exon_count"] = np.maximum(
        result["max_annotated_exons"] - result["direct_current_introner_count"], 1
    )
    result["log1p_n_isoforms"] = np.log1p(result["n_isoforms"])
    result["log1p_long_read_count"] = np.log1p(result["total_long_read_count"])
    result["log_gene_span"] = np.log(result["gene_span"])
    result["minimum_isoform_reads"] = minimum_isoform_reads
    return result


def fit_gee(
    data: pd.DataFrame,
    model_name: str,
    predictor: str,
    formula_suffix: str,
    outcome: str = "isoform",
) -> tuple[pd.DataFrame, dict[str, object]]:
    if outcome == "isoform":
        fraction = "nmd_isoform_fraction"
        denominator = "n_nmd_evaluable"
    elif outcome == "read":
        fraction = "nmd_read_fraction"
        denominator = "n_nmd_evaluable_reads"
    else:
        raise ValueError(outcome)

    required = [
        fraction,
        denominator,
        predictor,
        "strain",
        "common_gene_id",
    ]
    for term in re.findall(r"[A-Za-z_]\w*", formula_suffix):
        if term in data.columns:
            required.append(term)
    frame = data.replace([np.inf, -np.inf], np.nan).dropna(subset=set(required)).copy()
    frame = frame[frame[denominator] > 0].copy()
    formula = f"{fraction} ~ {predictor} + {formula_suffix}"
    # Unique isoforms are the binomial trials in the primary analysis. For the
    # read-abundance fraction, use a fractional-logit GEE with one weight per
    # gene-strain row; R2C2 reads are not independent biological replicates.
    analysis_weights = (
        frame[denominator] if outcome == "isoform" else pd.Series(1.0, index=frame.index)
    )
    model = sm.GEE.from_formula(
        formula,
        groups="common_gene_id",
        data=frame,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
        weights=analysis_weights,
    )
    result = model.fit(cov_type="robust", maxiter=200)
    conf = result.conf_int()
    coefficient_rows = []
    for term in result.params.index:
        coefficient_rows.append(
            {
                "model": model_name,
                "outcome": outcome,
                "predictor": predictor,
                "term": term,
                "estimate": result.params[term],
                "std_error": result.bse[term],
                "p_value": result.pvalues[term],
                "ci_lower": conf.loc[term, 0],
                "ci_upper": conf.loc[term, 1],
                "odds_ratio": np.exp(result.params[term]),
                "or_ci_lower": np.exp(conf.loc[term, 0]),
                "or_ci_upper": np.exp(conf.loc[term, 1]),
                "n_gene_strain_rows": len(frame),
                "n_genes": frame["common_gene_id"].nunique(),
                "converged": bool(result.converged),
            }
        )

    fitted = np.clip(np.asarray(result.fittedvalues), 1e-8, 1 - 1e-8)
    observed = frame[fraction].to_numpy()
    weights = analysis_weights.to_numpy()
    pearson = np.sum(weights * (observed - fitted) ** 2 / (fitted * (1 - fitted)))
    dispersion = pearson / max(len(frame) - len(result.params), 1)
    diagnostic = {
        "model": model_name,
        "outcome": outcome,
        "predictor": predictor,
        "formula": formula,
        "n_gene_strain_rows": len(frame),
        "n_genes": frame["common_gene_id"].nunique(),
        "n_introner_rows": int((frame[predictor] > 0).sum()),
        "n_nmd_positive": int(frame["n_nmd_positive"].sum()),
        "n_nmd_evaluable": int(frame["n_nmd_evaluable"].sum()),
        "pearson_dispersion": dispersion,
        "working_correlation": float(np.asarray(result.cov_struct.dep_params).squeeze()),
        "converged": bool(result.converged),
    }
    return pd.DataFrame(coefficient_rows), diagnostic


def make_plot(coefficients: pd.DataFrame, output: Path) -> None:
    wanted = coefficients[
        coefficients["model"].isin(
            [
                "binary_unadjusted",
                "binary_detection_adjusted",
                "binary_architecture_adjusted",
                "binary_raw_exon_adjusted_sensitivity",
                "burden_architecture_adjusted",
                "binary_plus_additional_burden",
                "binary_read_fractional_sensitivity",
                "binary_clean_mapping_sensitivity",
            ]
        )
        & (
            coefficients["term"].isin(
                ["has_introner", "direct_current_introner_count", "clean_has_introner"]
            )
        )
    ].copy()
    labels = {
        "binary_unadjusted": "Any introner: strain adjusted",
        "binary_detection_adjusted": "Any introner: detection adjusted",
        "binary_architecture_adjusted": "Any introner: architecture adjusted",
        "binary_raw_exon_adjusted_sensitivity": "Any introner: raw-exon sensitivity",
        "burden_architecture_adjusted": "Per introner: architecture adjusted",
        "binary_plus_additional_burden": "Any introner: burden-separated",
        "binary_read_fractional_sensitivity": "Any introner: read-fraction sensitivity",
        "binary_clean_mapping_sensitivity": "Any introner: clean-map sensitivity",
    }
    wanted["label"] = wanted["model"].map(labels)
    wanted = wanted.drop_duplicates("model").set_index("model").loc[list(labels)].reset_index()
    positions = np.arange(len(wanted))[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.errorbar(
        wanted["odds_ratio"],
        positions,
        xerr=[
            wanted["odds_ratio"] - wanted["or_ci_lower"],
            wanted["or_ci_upper"] - wanted["odds_ratio"],
        ],
        fmt="o",
        color="#2B6F8A",
        ecolor="#2B6F8A",
        capsize=3,
    )
    ax.axvline(1, color="#555555", linestyle="--", linewidth=1)
    ax.set_yticks(positions, wanted["label"])
    ax.set_xlabel("Odds ratio for NMD-positive status")
    ax.set_title("NMD propensity in introner-containing genes")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output = (
        args.output_dir.resolve()
        if args.output_dir
        else root / "analysis" / "nmd_introner_adjusted" / "results"
    )
    output.mkdir(parents=True, exist_ok=True)

    all_data = build_model_data(root, minimum_isoform_reads=0)
    primary = all_data[
        (all_data["total_long_read_count"] >= args.minimum_gene_read_support)
        & all_data["gene_span"].notna()
        & all_data["max_annotated_exons"].notna()
    ].copy()
    primary.to_csv(output / "nmd_gene_strain_data.tsv", sep="\t", index=False)

    specifications = [
        ("binary_unadjusted", "has_introner", "C(strain)", "isoform"),
        (
            "binary_detection_adjusted",
            "has_introner",
            "C(strain) + log1p_long_read_count + log1p_n_isoforms",
            "isoform",
        ),
        (
            "binary_architecture_adjusted",
            "has_introner",
            PRIMARY_FORMULA_SUFFIX,
            "isoform",
        ),
        (
            "burden_architecture_adjusted",
            "direct_current_introner_count",
            PRIMARY_FORMULA_SUFFIX,
            "isoform",
        ),
        (
            "binary_raw_exon_adjusted_sensitivity",
            "has_introner",
            RAW_EXON_FORMULA_SUFFIX,
            "isoform",
        ),
        (
            "binary_plus_additional_burden",
            "has_introner",
            "additional_introner_count + " + PRIMARY_FORMULA_SUFFIX,
            "isoform",
        ),
        (
            "binary_read_fractional_sensitivity",
            "has_introner",
            PRIMARY_FORMULA_SUFFIX,
            "read",
        ),
        (
            "binary_clean_mapping_sensitivity",
            "clean_has_introner",
            PRIMARY_FORMULA_SUFFIX,
            "isoform",
        ),
    ]
    coefficients, diagnostics = [], []
    for model_name, predictor, suffix, outcome in specifications:
        coef, diagnostic = fit_gee(primary, model_name, predictor, suffix, outcome)
        coefficients.append(coef)
        diagnostics.append(diagnostic)
    coefficients = pd.concat(coefficients, ignore_index=True)
    diagnostics = pd.DataFrame(diagnostics)
    coefficients.to_csv(output / "model_coefficients.tsv", sep="\t", index=False)
    diagnostics.to_csv(output / "model_diagnostics.tsv", sep="\t", index=False)

    sensitivity_rows = []
    for isoform_threshold in (0, 2, 5):
        threshold_data = build_model_data(root, minimum_isoform_reads=isoform_threshold)
        for gene_threshold in (5, 10, 20):
            frame = threshold_data[
                (threshold_data["total_long_read_count"] >= gene_threshold)
                & threshold_data["gene_span"].notna()
                & threshold_data["max_annotated_exons"].notna()
            ].copy()
            model_name = f"binary_isoform{isoform_threshold}_gene{gene_threshold}"
            coef, diagnostic = fit_gee(
                frame, model_name, "has_introner", PRIMARY_FORMULA_SUFFIX, "isoform"
            )
            focal = coef[coef["term"] == "has_introner"].iloc[0].to_dict()
            focal.update(
                {
                    "minimum_isoform_reads": isoform_threshold,
                    "minimum_gene_read_support": gene_threshold,
                    "pearson_dispersion": diagnostic["pearson_dispersion"],
                }
            )
            sensitivity_rows.append(focal)
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(output / "read_support_sensitivity.tsv", sep="\t", index=False)

    strain_rows = []
    for strain, frame in primary.groupby("strain", sort=False):
        coef, _ = fit_gee(
            frame,
            f"architecture_adjusted_{strain}",
            "has_introner",
            (
                "log1p_long_read_count + log_gene_span + "
                "background_exon_count + log1p_n_isoforms"
            ),
            "isoform",
        )
        strain_rows.append(coef[coef["term"] == "has_introner"])
    pd.concat(strain_rows, ignore_index=True).to_csv(
        output / "strain_stratified_sensitivity.tsv", sep="\t", index=False
    )

    mt_genes = extract_mating_type_genes(
        root / "results" / "annotations" / "CCMP1545.gtf"
    )
    checks = [
        (
            "one row per common gene and strain",
            not primary.duplicated(["common_gene_id", "strain"]).any(),
        ),
        (
            "NMD denominator equals positive plus negative isoforms",
            (
                primary["n_nmd_evaluable"]
                == primary["n_nmd_positive"] + primary["n_nmd_negative"]
            ).all(),
        ),
        (
            "binary introner indicator matches direct introner count",
            (
                primary["has_introner"]
                == (primary["direct_current_introner_count"] > 0).astype(int)
            ).all(),
        ),
        (
            "mating-type genes are excluded",
            not primary["common_gene_id"].isin(mt_genes).any(),
        ),
        ("all fitted GEEs converged", diagnostics["converged"].all()),
        (
            "all coefficient estimates and confidence limits are finite",
            np.isfinite(
                coefficients[["estimate", "ci_lower", "ci_upper"]].to_numpy()
            ).all(),
        ),
    ]
    validation_lines = ["Adjusted NMD model validation", "=============================", ""]
    validation_lines.extend(
        f"{'PASS' if passed else 'FAIL'}\t{label}" for label, passed in checks
    )
    (output / "validation_summary.txt").write_text(
        "\n".join(validation_lines) + "\n"
    )
    failed = [label for label, passed in checks if not passed]
    if failed:
        raise RuntimeError("Validation failed: " + "; ".join(failed))

    make_plot(coefficients, output / "adjusted_nmd_introner_effects.png")
    primary_focal = coefficients[
        (coefficients["model"] == "binary_architecture_adjusted")
        & (coefficients["term"] == "has_introner")
    ].iloc[0]
    burden_focal = coefficients[
        (coefficients["model"] == "burden_architecture_adjusted")
        & (coefficients["term"] == "direct_current_introner_count")
    ].iloc[0]
    weighted_focal = coefficients[
        (coefficients["model"] == "binary_read_fractional_sensitivity")
        & (coefficients["term"] == "has_introner")
    ].iloc[0]
    summary = [
        "Adjusted NMD-introner association analysis",
        "==========================================",
        "",
        "Primary model: population-averaged binomial GEE clustered by common gene.",
        "Outcome: NMD-positive versus NMD-negative multi-exon SQANTI isoforms.",
        "Adjustment: strain, total R2C2 support, total detected isoform count,",
        "gene span, and background exon count after subtracting current introners.",
        "",
        f"Gene-strain rows: {len(primary)}",
        f"Common genes: {primary['common_gene_id'].nunique()}",
        f"NMD-evaluable isoforms: {int(primary['n_nmd_evaluable'].sum())}",
        f"NMD-positive isoforms: {int(primary['n_nmd_positive'].sum())}",
        f"Introner-containing gene-strain rows: {int(primary['has_introner'].sum())}",
        "",
        "Focal effects:",
        (
            "  Any current introner: "
            f"OR={primary_focal.odds_ratio:.4f} "
            f"({primary_focal.or_ci_lower:.4f}-{primary_focal.or_ci_upper:.4f}), "
            f"p={primary_focal.p_value:.4g}"
        ),
        (
            "  Per current introner: "
            f"OR={burden_focal.odds_ratio:.4f} "
            f"({burden_focal.or_ci_lower:.4f}-{burden_focal.or_ci_upper:.4f}), "
            f"p={burden_focal.p_value:.4g}"
        ),
        (
            "  Read-fraction any introner: "
            f"OR={weighted_focal.odds_ratio:.4f} "
            f"({weighted_focal.or_ci_lower:.4f}-{weighted_focal.or_ci_upper:.4f}), "
            f"p={weighted_focal.p_value:.4g}"
        ),
        "",
        "Interpretation:",
        "  The binomial denominator controls the number of NMD-evaluable isoforms.",
        "  log1p_n_isoforms additionally adjusts the mean model for total detected",
        "  isoform diversity. The association remains observational and cross-gene.",
        "  The raw maximum-exon model is retained only as a sensitivity because a",
        "  spliced introner can itself increase the annotated exon count.",
    ]
    (output / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()

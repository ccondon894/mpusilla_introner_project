"""Validate internal consistency of the isoform-turnover workflow outputs."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def check(condition, message, results):
    results.append((bool(condition), message))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    locus_map = pd.read_csv(args.data_dir / "locus_gene_map.tsv", sep="\t")
    events = pd.read_csv(args.data_dir / "locus_level_reference_comparisons.tsv", sep="\t")
    model_data = pd.read_csv(args.data_dir / "isoform_introner_model_data.tsv", sep="\t")
    coefficients = pd.read_csv(args.model_dir / "model_coefficients.tsv", sep="\t")
    sidecar = pd.read_csv(
        args.sidecar,
        sep="\t",
        usecols=["ortholog_id", "gene"],
    ).drop_duplicates("ortholog_id")

    results = []
    check(
        not locus_map.duplicated(["ortholog_id", "sample"]).any(),
        "one locus-map row per ortholog and sample",
        results,
    )
    check(set(locus_map["presence"].unique()) <= {1, 2, 3}, "genotype calls use 1/2/3", results)
    check(
        not events.loc[events["event_status"] == "uncallable", ["gain_count", "loss_count"]]
        .astype(bool)
        .any(axis=None),
        "missing comparisons never become gains or losses",
        results,
    )
    check(
        not model_data.duplicated(["common_gene_id", "strain"]).any(),
        "one model row per common gene and strain",
        results,
    )
    check(
        (model_data["n_extra_isoforms"] == model_data["n_isoforms"] - 1).all(),
        "extra-isoform outcome equals total minus one",
        results,
    )
    check(
        ((model_data["reference_relative_gain_count"] > 0)
         & (model_data["reference_relative_loss_count"] > 0)).any(),
        "locus-wise preparation permits simultaneous gain and loss",
        results,
    )
    check(
        np.isfinite(coefficients[["coefficient", "std_error", "p_value"]]).all(axis=None),
        "all saved model coefficients are finite",
        results,
    )
    within_gene_models = coefficients[
        coefficients["model"].isin(
            [
                "within_gene_conditional_poisson_all",
                "within_gene_conditional_logit_all",
            ]
        )
    ]
    check(
        set(within_gene_models["model"])
        == {
            "within_gene_conditional_poisson_all",
            "within_gene_conditional_logit_all",
        }
        and (within_gene_models["term"] == "current_introner_count").groupby(
            within_gene_models["model"]
        ).any().all(),
        "within-gene models include current introner count",
        results,
    )
    check(
        not within_gene_models["term"].isin({"gain", "loss"}).any(),
        "within-gene models do not use reference-relative gain/loss terms",
        results,
    )

    reference_map = locus_map.loc[
        locus_map["sample"] == "CCMP1545", ["ortholog_id", "native_gene_id"]
    ].merge(sidecar, on="ortholog_id", how="inner")
    # The overlap sidecar can retain historical loci that were subsequently
    # removed from the final genotype matrix (for example, low-identity groups).
    # Validate the mappings that remain in the current workflow instead of
    # requiring filtered loci to be resurrected solely for validation.
    sidecar_loci_in_workflow = sidecar["ortholog_id"].isin(locus_map["ortholog_id"])
    check(
        len(reference_map) == int(sidecar_loci_in_workflow.sum()),
        "all current clean sidecar loci have workflow mappings",
        results,
    )
    check(
        (reference_map["native_gene_id"] == reference_map["gene"]).all(),
        "workflow CCMP1545 mappings agree with the clean overlap sidecar",
        results,
    )

    lines = ["Isoform turnover workflow validation", "===================================="]
    for passed, message in results:
        lines.append(f"{'PASS' if passed else 'FAIL'}\t{message}")
    lines.extend(
        [
            "",
            f"Model rows: {len(model_data)}",
            f"Common genes: {model_data['common_gene_id'].nunique()}",
            "Rows with both gain and loss: "
            f"{int(((model_data['reference_relative_gain_count'] > 0) & (model_data['reference_relative_loss_count'] > 0)).sum())}",
            "Common genes with within-gene introner-count variation: "
            f"{int((model_data.groupby('common_gene_id')['current_introner_count'].nunique() > 1).sum())}",
            f"Clean sidecar mappings checked: {len(reference_map)}",
            f"Historical sidecar loci absent from current workflow: {int((~sidecar_loci_in_workflow).sum())}",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    if not all(passed for passed, _ in results):
        raise SystemExit("Isoform turnover validation failed; see validation_summary.txt")


if __name__ == "__main__":
    main()

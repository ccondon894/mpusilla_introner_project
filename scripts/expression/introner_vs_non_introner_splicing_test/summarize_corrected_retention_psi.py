#!/usr/bin/env python3
"""Create evidence-guarded summaries from corrected retention-PSI outputs."""

from pathlib import Path

import pandas as pd
from scipy import stats


import argparse

OUT = Path("results/expression/functional/introner_vs_non_introner_splicing_test/corrected_retention_psi")


def summarize_strata(loci: pd.DataFrame) -> pd.DataFrame:
    eligible = loci.primary_boundary_eligible.astype(bool)
    strata = {
        "any_informative_signal": eligible & (loci.total_psi_signal > 0),
        "total_signal_ge_10": eligible & (loci.total_psi_signal >= 10),
        "junction_signal_gt_0": eligible & (loci.junction_signal_exact > 0),
        "junction_signal_ge_5": eligible & (loci.junction_signal_exact >= 5),
        "junction_signal_ge_10": eligible & (loci.junction_signal_exact >= 10),
    }
    rows = []
    for stratum, mask in strata.items():
        for cls, x in loci[mask].groupby("analysis_class", sort=False):
            rows.append({
                "stratum": stratum,
                "analysis_class": cls,
                "n_loci": len(x),
                "mean_retention_psi": x.retention_psi.mean(),
                "median_retention_psi": x.retention_psi.median(),
                "q25_retention_psi": x.retention_psi.quantile(.25),
                "q75_retention_psi": x.retention_psi.quantile(.75),
                "low_retention_psi_fraction": (x.retention_psi <= .10).mean(),
                "mean_retained_support": x.retained_support.mean(),
                "median_junction_signal": x.junction_signal_exact.median(),
            })
    return pd.DataFrame(rows)


def summarize_replicates(rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, replicate, cls), x in rep.groupby(
            ["sample", "replicate", "analysis_class"], sort=False):
        confirmed = x[
            x.primary_boundary_eligible.astype(bool)
            & (x.S_junction_signal_exact >= 5)
        ]
        rows.append({
            "sample": sample, "replicate": replicate, "analysis_class": cls,
            "n_loci": len(x),
            "n_with_any_signal": int((x.informative_signal > 0).sum()),
            "n_junction_confirmed_ge_5": len(confirmed),
            "mean_U5": x.U5_donor_crossing.mean(),
            "mean_U3": x.U3_acceptor_crossing.mean(),
            "mean_S": x.S_junction_signal_exact.mean(),
            "mean_psi_junction_confirmed_ge_5": confirmed.retention_psi.mean(),
            "median_psi_junction_confirmed_ge_5": confirmed.retention_psi.median(),
            "low_retention_fraction_junction_confirmed_ge_5": (
                confirmed.retention_psi <= .10).mean(),
        })
    return pd.DataFrame(rows)


def main() -> None:
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    OUT = parser.parse_args().outdir
    loci = pd.read_csv(OUT / "per_locus_scores.tsv", sep="\t", low_memory=False)
    rep = pd.read_csv(OUT / "per_locus_replicate_U5_U3_S.tsv", sep="\t",
                      low_memory=False)
    strata = summarize_strata(loci)
    rep_summary = summarize_replicates(rep)
    boundary = []
    for cls, x in loci.groupby("analysis_class", sort=False):
        with_unspliced = x[(x.donor_crossing + x.acceptor_crossing) > 0]
        boundary.append({
            "analysis_class": cls,
            "n_loci": len(x),
            "n_with_unspliced_boundary_support": len(with_unspliced),
            "both_boundaries_observed_fraction": (
                (with_unspliced.donor_crossing > 0)
                & (with_unspliced.acceptor_crossing > 0)).mean(),
            "median_boundary_agreement_ratio": (
                with_unspliced.boundary_agreement_ratio.median()),
            "agreement_ratio_ge_0_5_fraction": (
                with_unspliced.boundary_agreement_ratio >= .5).mean(),
        })
    strata.to_csv(OUT / "evidence_guarded_class_summary.tsv", sep="\t", index=False)
    rep_summary.to_csv(OUT / "sample_replicate_class_summary.tsv", sep="\t", index=False)
    pd.DataFrame(boundary).to_csv(
        OUT / "boundary_agreement_audit.tsv", sep="\t", index=False)
    confirmed = loci[
        loci.primary_boundary_eligible.astype(bool)
        & (loci.junction_signal_exact >= 5)
    ]
    tests = []
    for focal, reference in [
        ("fixed_present", "conventional_intron"),
        ("polymorphic", "conventional_intron"),
        ("polymorphic", "fixed_present"),
    ]:
        a = confirmed.loc[confirmed.analysis_class == focal, "retention_psi"]
        b = confirmed.loc[confirmed.analysis_class == reference, "retention_psi"]
        test = stats.mannwhitneyu(a, b, alternative="two-sided")
        tests.append({
            "eligibility": "junction_signal_ge_5",
            "focal_class": focal, "reference_class": reference,
            "n_focal": len(a), "n_reference": len(b),
            "median_focal": a.median(), "median_reference": b.median(),
            "mean_focal": a.mean(), "mean_reference": b.mean(),
            "mann_whitney_u": test.statistic, "mann_whitney_p": test.pvalue,
        })
    pd.DataFrame(tests).to_csv(
        OUT / "evidence_guarded_class_tests.tsv", sep="\t", index=False)
    primary_summary = strata[strata.stratum == "junction_signal_ge_5"]
    primary_tests = pd.DataFrame(tests)
    resolver = pd.read_csv(OUT / "boundary_resolver_audit.tsv", sep="\t")
    with (OUT / "PRIMARY_RESULTS.md").open("w") as out:
        out.write("# Regtools-resolved boundary-crossing PSI\n\n")
        out.write(
            "Introner boundaries require an exact regtools junction pair observed "
            "in at least two replicates, with pooled support >=10 and both "
            "endpoints within 25 bp of a plausible predicted interval. Duplicate "
            "resolved genomic coordinates are excluded from the primary analysis. "
            "Conventional introns use exact GTF boundaries. The mating-type region "
            "is excluded.\n\n## Boundary-resolution audit\n\n"
        )
        out.write(resolver.to_markdown(index=False))
        out.write("\n\n## Primary class summary\n\n")
        out.write(primary_summary.to_markdown(index=False))
        out.write("\n\n## Descriptive locus-level tests\n\n")
        out.write(primary_tests.to_markdown(index=False))
        out.write("\n")


if __name__ == "__main__":
    main()

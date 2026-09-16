#!/usr/bin/env python3
"""Approximate the exon-normalized intron-depth PSI used by Gozashti et al.

The paper describes PSI as relative RNA-seq read depth across an intron versus
its adjacent exons, but does not provide a directly reusable implementation.
For each locus and replicate this script therefore sums aligned base depth over
the full intron and over equally sized windows immediately outside both splice
boundaries. Replicates are pooled as base counts before calculating

    depth_PSI = min(1, mean_intron_depth / mean_adjacent_exon_depth).

The 50-bp window is primary; 20- and 100-bp estimates are sensitivity checks.
This isolated analysis uses the same locus construction and mating-type-region
exclusion as corrected_retention_psi.py and does not modify workflow outputs.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from corrected_retention_psi import construct_loci


WINDOWS = (20, 50, 100)


def arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", type=Path,
                   default=Path("results/genotyping/genotype_matrix.final.tsv"))
    p.add_argument("--annotations", type=Path, default=Path("results/annotations"))
    p.add_argument("--splice-dir", type=Path,
                   default=Path("results/expression/splice_junctions"))
    p.add_argument("--bam-dir", type=Path,
                   default=Path("results/expression/alignments/bams"))
    p.add_argument("--outdir", type=Path, default=Path(
        "analysis/introner_vs_non_introner_splicing_test/paper_depth_psi"))
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--min-exon-depth", type=float, default=5.0)
    return p.parse_args()


def keep_read(read) -> bool:
    return not (
        read.is_unmapped or read.is_secondary or read.is_supplementary
        or read.is_duplicate or read.is_qcfail
    )


def scan_bam(task: tuple[str, str, str, list[tuple[int, str, int, int]]]) -> dict:
    import pysam

    sample, replicate, bam_name, records = task
    rows = []
    with pysam.AlignmentFile(bam_name, "rb") as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        for row_id, contig, intron_start, intron_end in records:
            contig_len = lengths.get(contig)
            if contig_len is None or intron_start < 0 or intron_end > contig_len:
                rows.append((row_id, np.nan, *([np.nan, np.nan] * len(WINDOWS))))
                continue
            span_start = max(0, intron_start - max(WINDOWS))
            span_end = min(contig_len, intron_end + max(WINDOWS))
            cov = bam.count_coverage(
                contig, span_start, span_end, quality_threshold=0,
                read_callback=keep_read,
            )
            depth = np.asarray(cov[0], dtype=np.int64)
            for base in cov[1:]:
                depth += np.asarray(base, dtype=np.int64)
            a, b = intron_start - span_start, intron_end - span_start
            intron_bases = int(depth[a:b].sum())
            values: list[float | int] = [row_id, intron_bases]
            for window in WINDOWS:
                left_start = max(span_start, intron_start - window) - span_start
                right_end = min(span_end, intron_end + window) - span_start
                exon_bases = int(depth[left_start:a].sum() + depth[b:right_end].sum())
                exon_length = (a - left_start) + (right_end - b)
                values.extend((exon_bases, exon_length))
            rows.append(tuple(values))
    return {"sample": sample, "replicate": replicate, "rows": rows}


def class_summary(x: pd.DataFrame, min_exon_depth: float) -> pd.DataFrame:
    rows = []
    for window in WINDOWS:
        psi = f"depth_psi_{window}bp"
        exon = f"exon_depth_{window}bp"
        eligible = x[x[exon] >= min_exon_depth]
        for cls, z in eligible.groupby("analysis_class", sort=False):
            rows.append({
                "window_bp": window,
                "analysis_class": cls,
                "min_exon_depth": min_exon_depth,
                "n_loci": len(z),
                "n_genes": z.gene_id.nunique(),
                "mean_depth_psi": z[psi].mean(),
                "median_depth_psi": z[psi].median(),
                "q25_depth_psi": z[psi].quantile(.25),
                "q75_depth_psi": z[psi].quantile(.75),
                "fraction_psi_le_0.10": (z[psi] <= .10).mean(),
                "mean_exon_depth": z[exon].mean(),
            })
    return pd.DataFrame(rows)


def class_tests(x: pd.DataFrame, min_exon_depth: float) -> pd.DataFrame:
    comparisons = [
        ("polymorphic", "fixed_present"),
        ("fixed_present", "conventional_intron"),
        ("polymorphic", "conventional_intron"),
    ]
    rows = []
    for window in WINDOWS:
        psi, exon = f"depth_psi_{window}bp", f"exon_depth_{window}bp"
        y = x[x[exon] >= min_exon_depth]
        for left, right in comparisons:
            a = y.loc[y.analysis_class == left, psi].dropna()
            b = y.loc[y.analysis_class == right, psi].dropna()
            test = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else None
            rows.append({
                "window_bp": window, "class_a": left, "class_b": right,
                "n_a": len(a), "n_b": len(b),
                "median_a": a.median(), "median_b": b.median(),
                "median_difference_a_minus_b": a.median() - b.median(),
                "mann_whitney_u": test.statistic if test else np.nan,
                "mann_whitney_p": test.pvalue if test else np.nan,
            })
    return pd.DataFrame(rows)


def main() -> None:
    cfg = arguments()
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    loci, audit = construct_loci(cfg)
    loci = loci[loci.valid_boundaries].copy().reset_index(drop=True)
    records_by_sample = {
        sample: [(i, r.contig, int(r.start), int(r.end))
                 for i, r in loci[loci["sample"] == sample].iterrows()]
        for sample in loci["sample"].unique()
    }

    # Import here to keep the locus-construction source of truth in one place.
    from compare_to_conventional_introns import REPLICATES
    tasks = [
        (sample, rep, str(cfg.bam_dir / f"{rep}.sorted.bam"), records_by_sample[sample])
        for sample in records_by_sample for rep in REPLICATES[sample]
    ]
    results = []
    with ProcessPoolExecutor(max_workers=cfg.workers) as pool:
        futures = [pool.submit(scan_bam, task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"finished {result['sample']} {result['replicate']}", flush=True)

    intron_bases = np.zeros(len(loci), dtype=np.float64)
    intron_observations = np.zeros(len(loci), dtype=np.int64)
    exon_bases = {w: np.zeros(len(loci), dtype=np.float64) for w in WINDOWS}
    exon_lengths = {w: np.zeros(len(loci), dtype=np.float64) for w in WINDOWS}
    replicate_rows = []
    for result in results:
        for values in result["rows"]:
            row_id, ib = values[:2]
            if np.isnan(ib):
                continue
            intron_bases[row_id] += ib
            intron_observations[row_id] += 1
            rec = {"row_id": row_id, "sample": result["sample"],
                   "replicate": result["replicate"], "intron_aligned_bases": ib}
            pos = 2
            for w in WINDOWS:
                eb, el = values[pos:pos + 2]
                pos += 2
                exon_bases[w][row_id] += eb
                exon_lengths[w][row_id] += el
                rec[f"exon_aligned_bases_{w}bp"] = eb
                rec[f"exon_callable_bases_{w}bp"] = el
            replicate_rows.append(rec)

    loci["n_bam_replicates"] = intron_observations
    loci["intron_aligned_bases"] = intron_bases
    loci["intron_depth"] = intron_bases / (
        (loci.end - loci.start).to_numpy() * intron_observations)
    for w in WINDOWS:
        loci[f"exon_aligned_bases_{w}bp"] = exon_bases[w]
        loci[f"exon_callable_bases_{w}bp"] = exon_lengths[w]
        loci[f"exon_depth_{w}bp"] = exon_bases[w] / exon_lengths[w]
        ratio = loci.intron_depth / loci[f"exon_depth_{w}bp"]
        loci[f"depth_psi_uncapped_{w}bp"] = ratio
        loci[f"depth_psi_{w}bp"] = ratio.clip(upper=1)

    summary = class_summary(loci, cfg.min_exon_depth)
    tests = class_tests(loci, cfg.min_exon_depth)
    validation = []
    for w in WINDOWS:
        eligible = loci[loci[f"exon_depth_{w}bp"] >= cfg.min_exon_depth]
        validation.append({
            "window_bp": w,
            "eligible_loci": len(eligible),
            "uncapped_ratio_gt_1": int((eligible[f"depth_psi_uncapped_{w}bp"] > 1).sum()),
            "psi_outside_0_1": int(((eligible[f"depth_psi_{w}bp"] < 0) | (eligible[f"depth_psi_{w}bp"] > 1)).sum()),
            "missing_four_replicates": int((eligible.n_bam_replicates != 4).sum()),
        })

    loci.to_csv(cfg.outdir / "per_locus_depth_psi.tsv", sep="\t", index=False)
    pd.DataFrame(replicate_rows).to_csv(
        cfg.outdir / "per_locus_replicate_depth_counts.tsv", sep="\t", index=False)
    summary.to_csv(cfg.outdir / "class_summary.tsv", sep="\t", index=False)
    tests.to_csv(cfg.outdir / "class_tests.tsv", sep="\t", index=False)
    audit.to_csv(cfg.outdir / "locus_selection_audit.tsv", sep="\t", index=False)
    pd.DataFrame(validation).to_csv(cfg.outdir / "validation.tsv", sep="\t", index=False)

    primary = summary[summary.window_bp == 50]
    lines = [
        "# Exon-normalized intron-depth PSI approximation", "",
        "PSI is the pooled mean intron depth divided by pooled mean depth in the two immediate 50-bp exonic flanks, capped at 1. Loci require mean exon depth >= 5. Windows of 20 and 100 bp are sensitivity checks.", "",
        "This reproduces the published text description, not unavailable original analysis code. It uses the project's existing HISAT2 alignments rather than remapping with STAR.", "",
        "## Primary 50-bp results", "",
        primary.to_markdown(index=False), "", "## Pairwise tests", "",
        tests[tests.window_bp == 50].to_markdown(index=False), "",
    ]
    (cfg.outdir / "RESULTS.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

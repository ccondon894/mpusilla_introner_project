#!/usr/bin/env python3
"""Coverage-conditioned retention PSI using resolved boundary-crossing reads.

The original exploratory comparison counted any read with >=8 aligned bases
inside an intron as retained.  That is not a retention observation: most such
reads end inside the intron and never demonstrate continuity across a splice
boundary.  Here, an unspliced alignment supports retention only when one
continuous aligned block spans a donor or acceptor boundary with >=8 aligned
bases on each side.  Pooled retained support is

    (donor-crossing alignments + acceptor-crossing alignments) / 2

and PSI is retained_support / (retained_support + junction_signal).

For introners, exact regtools junction pairs are pooled across replicates and
used to resolve fuzzy predicted boundaries. A pair is accepted when both ends
are close to a plausible predicted interval and it passes replicate and pooled-
support guards. Conventional introns retain exact GTF exon-gap boundaries.

This is an isolated analysis and does not alter workflow outputs.
"""

from __future__ import annotations

import argparse
import bisect
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from compare_to_conventional_introns import (
    GROUP,
    MATING,
    REPLICATES,
    SAMPLES,
    IntervalIndex,
    classify,
    host_gene,
    overlap,
    parse_gtf,
    parse_junctions,
    present_bodies,
)


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
        "analysis/introner_vs_non_introner_splicing_test/corrected_retention_psi"))
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--min-anchor", type=int, default=8)
    p.add_argument("--boundary-search-bp", type=int, default=25)
    p.add_argument("--min-resolver-replicates", type=int, default=2)
    p.add_argument("--min-resolver-score", type=int, default=10)
    return p.parse_args()


def construct_loci(cfg: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reproduce the previous clean control and focal-introner selection."""
    matrix = pd.read_csv(cfg.matrix, sep="\t", low_memory=False)
    for col in ("start", "end", "presence"):
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")

    mt = matrix[
        (matrix["sample"] == "CCMP1545")
        & (matrix.contig == MATING["CCMP1545"][0])
        & (matrix.start < MATING["CCMP1545"][2])
        & (matrix.end > MATING["CCMP1545"][1])
    ]
    mt_orthologs = set(mt.ortholog_id)
    loci_parts: list[pd.DataFrame] = []
    audit_rows: list[dict] = []

    for sample in SAMPLES:
        introns, genes, gtf_audit = parse_gtf(cfg.annotations / f"{sample}.gtf")
        gene_idx = {
            contig: IntervalIndex([
                (int(r.start), int(r.end), r.gene_id) for r in part.itertuples()
            ])
            for contig, part in genes.groupby("contig")
        }
        bodies = present_bodies(matrix, sample)
        bodies = bodies[~bodies.ortholog_id.isin(mt_orthologs)].copy()
        bodies["analysis_class"] = bodies.apply(classify, axis=1)
        bodies["gene_id"] = bodies.apply(lambda r: host_gene(r, gene_idx), axis=1)
        primary = bodies[
            (bodies.analysis_class != "excluded") & bodies.gene_id.ne("")
        ].copy()

        # The raw scorer's expected coordinates are authoritative for introners.
        raw_path = (cfg.splice_dir / "introner_boundary_support"
                    / f"{sample}.per_locus.tsv")
        raw = pd.read_csv(raw_path, sep="\t", low_memory=False)
        focal = raw.merge(
            primary[["ortholog_id", "gene_id", "analysis_class"]],
            on="ortholog_id", validate="1:1", suffixes=("", "_host"),
        )
        focal["start"] = pd.to_numeric(focal.expected_splice_start, errors="coerce")
        focal["end"] = pd.to_numeric(focal.expected_splice_end, errors="coerce")
        focal["nominal_start"] = focal["start"]
        focal["nominal_end"] = focal["end"]
        focal["candidate_body_start"] = pd.to_numeric(
            focal.body_start, errors="coerce")
        focal["candidate_body_end"] = pd.to_numeric(
            focal.body_end, errors="coerce")
        focal["candidate_donor_offset"] = pd.to_numeric(
            focal.splice_donor_offset, errors="coerce")
        focal["candidate_acceptor_offset"] = pd.to_numeric(
            focal.splice_acceptor_offset, errors="coerce")
        focal["locus_id"] = focal.ortholog_id
        focal["boundary_source"] = "introner_nominal_unresolved"
        focal["old_retained_signal"] = pd.to_numeric(
            focal.retained_signal, errors="coerce")
        focal["old_retention_psi"] = pd.to_numeric(
            focal.retention_psi, errors="coerce")

        host_genes = set(primary.gene_id)
        body_idx = {
            contig: IntervalIndex([
                (int(r.body_start), int(r.body_end), "") for r in part.itertuples()
            ])
            for contig, part in bodies.groupby("contig")
        }
        mcontig, mstart, mend = MATING[sample]
        introns["in_host_gene"] = introns.gene_id.isin(host_genes)
        introns["overlaps_present_introner"] = introns.apply(
            lambda r: body_idx.get(r.contig, IntervalIndex([])).overlaps(
                int(r.start), int(r.end)), axis=1)
        introns["mating_region"] = introns.apply(
            lambda r: r.contig == mcontig and overlap(
                int(r.start), int(r.end), mstart, mend), axis=1)
        introns["length_compatible"] = (
            introns.end - introns.start).between(20, 500000)
        clean = introns[
            introns.in_host_gene
            & ~introns.overlaps_present_introner
            & ~introns.mating_region
            & introns.length_compatible
        ].copy()
        clean["sample"] = sample
        clean["analysis_class"] = "conventional_intron"
        clean["locus_id"] = clean.apply(
            lambda r: f"{sample}|{r.contig}:{r.start}-{r.end}", axis=1)
        clean["boundary_source"] = "exact_gtf_exon_gap"
        clean["nominal_start"] = clean["start"]
        clean["nominal_end"] = clean["end"]
        clean["candidate_body_start"] = np.nan
        clean["candidate_body_end"] = np.nan
        clean["candidate_donor_offset"] = np.nan
        clean["candidate_acceptor_offset"] = np.nan
        clean["old_retained_signal"] = np.nan
        clean["old_retention_psi"] = np.nan

        common = ["locus_id", "sample", "gene_id", "contig", "start", "end",
                  "analysis_class", "boundary_source", "old_retained_signal",
                  "old_retention_psi", "nominal_start", "nominal_end",
                  "candidate_body_start", "candidate_body_end",
                  "candidate_donor_offset", "candidate_acceptor_offset"]
        loci_parts.extend([clean[common], focal[common]])
        focal_valid = focal.start.notna() & focal.end.notna() & (focal.start < focal.end)
        audit_rows.append({
            "sample": sample,
            **gtf_audit,
            "present_introner_bodies_after_mating_exclusion": len(bodies),
            "selected_fixed_introner_loci": int((focal.analysis_class == "fixed_present").sum()),
            "selected_polymorphic_introner_loci": int((focal.analysis_class == "polymorphic").sum()),
            "introner_loci_with_valid_expected_boundaries": int(focal_valid.sum()),
            "introner_loci_missing_or_invalid_expected_boundaries": int((~focal_valid).sum()),
            "clean_conventional_introns": len(clean),
            "control_genes": clean.gene_id.nunique(),
        })

    loci = pd.concat(loci_parts, ignore_index=True)
    loci["start"] = pd.to_numeric(loci.start, errors="coerce")
    loci["end"] = pd.to_numeric(loci.end, errors="coerce")
    loci["valid_boundaries"] = (
        loci.start.notna() & loci.end.notna() & (loci.start < loci.end))
    return loci, pd.DataFrame(audit_rows)


def plausible_introner_intervals(rec: pd.Series) -> list[tuple[str, int, int]]:
    """Reconstruct the body, motif-offset, and selected nominal candidates."""
    targets: list[tuple[str, int, int]] = []
    seen: set[tuple[int, int]] = set()

    def add(source: str, start: float, end: float) -> None:
        if pd.isna(start) or pd.isna(end):
            return
        a, b = int(start), int(end)
        if a < b and (a, b) not in seen:
            targets.append((source, a, b))
            seen.add((a, b))

    add("nominal_expected", rec.nominal_start, rec.nominal_end)
    add("body", rec.candidate_body_start, rec.candidate_body_end)
    if not any(pd.isna(x) for x in (
        rec.candidate_body_start, rec.candidate_body_end,
        rec.candidate_donor_offset, rec.candidate_acceptor_offset,
    )):
        body_start, body_end = int(rec.candidate_body_start), int(rec.candidate_body_end)
        donor, acceptor = int(rec.candidate_donor_offset), int(rec.candidate_acceptor_offset)
        plus = (body_start + donor, body_start + acceptor + 2)
        reverse = (body_end - acceptor - 2, body_end - donor)
        if body_start <= plus[0] < plus[1] <= body_end:
            add("splice_site_plus", *plus)
        if body_start <= reverse[0] < reverse[1] <= body_end:
            add("splice_site_reverse_complement", *reverse)
    return targets


def resolve_boundaries(
    loci: pd.DataFrame, cfg: argparse.Namespace
) -> tuple[pd.DataFrame, dict[tuple[str, str], np.ndarray], pd.DataFrame]:
    """Resolve introner boundaries using pooled exact regtools junction pairs."""
    out = loci.copy()
    n = len(out)
    exact_signal = np.zeros(n, dtype=np.int64)
    exact_replicates = np.zeros(n, dtype=np.int64)
    by_replicate = {
        (sample, rep): np.zeros(n, dtype=np.int64)
        for sample in SAMPLES for rep in REPLICATES[sample]
    }
    out["boundary_resolved_flag"] = False
    out["primary_boundary_eligible"] = False
    out["resolver_collision_flag"] = False
    out["resolver_target_source"] = ""
    out["resolver_start_delta"] = np.nan
    out["resolver_end_delta"] = np.nan
    audit_rows = []

    for sample in SAMPLES:
        junctions = parse_junctions([
            (rep, cfg.splice_dir / "regtools" / "per_replicate"
             / f"{rep}.junctions.bed")
            for rep in REPLICATES[sample]
        ])
        by_start: dict[str, dict[int, list[tuple[int, list[tuple[int, str]]]]]] = {}
        for contig, pairs in junctions.items():
            start_index: dict[int, list[tuple[int, list[tuple[int, str]]]]] = {}
            for (start, end), evidence in pairs.items():
                start_index.setdefault(start, []).append((end, evidence))
            by_start[contig] = start_index

        for row_id, rec in out[out["sample"] == sample].iterrows():
            if rec.analysis_class == "conventional_intron":
                out.at[row_id, "boundary_resolved_flag"] = True
                out.at[row_id, "primary_boundary_eligible"] = True
                out.at[row_id, "resolver_target_source"] = "exact_gtf_exon_gap"
                evidence = junctions.get(rec.contig, {}).get(
                    (int(rec.start), int(rec.end)), [])
            else:
                targets = plausible_introner_intervals(rec)
                candidates: dict[tuple[int, int], tuple] = {}
                start_index = by_start.get(rec.contig, {})
                for source, expected_start, expected_end in targets:
                    for start in range(
                        expected_start - cfg.boundary_search_bp,
                        expected_start + cfg.boundary_search_bp + 1,
                    ):
                        for end, pair_evidence in start_index.get(start, []):
                            end_delta = end - expected_end
                            if abs(end_delta) > cfg.boundary_search_bp:
                                continue
                            score = sum(x[0] for x in pair_evidence)
                            reps = len({x[1] for x in pair_evidence if x[0] > 0})
                            if reps < cfg.min_resolver_replicates or score < cfg.min_resolver_score:
                                continue
                            start_delta = start - expected_start
                            rank = (
                                -reps, -score,
                                max(abs(start_delta), abs(end_delta)),
                                abs(start_delta) + abs(end_delta),
                                source, start, end,
                            )
                            key = (start, end)
                            if key not in candidates or rank < candidates[key][0]:
                                candidates[key] = (
                                    rank, source, start_delta, end_delta, pair_evidence)

                if candidates:
                    (start, end), value = min(
                        candidates.items(), key=lambda item: item[1][0])
                    _, source, start_delta, end_delta, evidence = value
                    out.at[row_id, "start"] = start
                    out.at[row_id, "end"] = end
                    out.at[row_id, "boundary_source"] = "regtools_pooled_exact"
                    out.at[row_id, "boundary_resolved_flag"] = True
                    out.at[row_id, "primary_boundary_eligible"] = True
                    out.at[row_id, "resolver_target_source"] = source
                    out.at[row_id, "resolver_start_delta"] = start_delta
                    out.at[row_id, "resolver_end_delta"] = end_delta
                else:
                    evidence = []

            exact_signal[row_id] = sum(score for score, _ in evidence)
            exact_replicates[row_id] = len(
                {rep for score, rep in evidence if score > 0})
            for score, rep in evidence:
                by_replicate[(sample, rep)][row_id] += score

        introner_resolved = (
            out["sample"].eq(sample)
            & out.analysis_class.ne("conventional_intron")
            & out.boundary_resolved_flag.astype(bool)
        )
        collision = out.loc[introner_resolved].duplicated(
            ["contig", "start", "end"], keep=False)
        collision_ids = collision[collision].index
        out.loc[collision_ids, "resolver_collision_flag"] = True
        out.loc[collision_ids, "primary_boundary_eligible"] = False

        sample_part = out[out["sample"] == sample]
        for cls, part in sample_part.groupby("analysis_class", sort=False):
            audit_rows.append({
                "sample": sample,
                "analysis_class": cls,
                "n_loci": len(part),
                "n_boundary_resolved": int(part.boundary_resolved_flag.sum()),
                "resolved_fraction": part.boundary_resolved_flag.mean(),
                "n_resolved_coordinate_collisions": int(
                    part.resolver_collision_flag.sum()),
                "n_primary_boundary_eligible": int(
                    part.primary_boundary_eligible.sum()),
                "resolver_search_bp": cfg.boundary_search_bp,
                "resolver_min_replicates": cfg.min_resolver_replicates,
                "resolver_min_pooled_score": cfg.min_resolver_score,
            })

    out["valid_boundaries"] = (
        out.start.notna() & out.end.notna() & (out.start < out.end))
    out["junction_signal_exact"] = exact_signal
    out["junction_replicate_count_exact"] = exact_replicates
    return out, by_replicate, pd.DataFrame(audit_rows)


def scan_one_bam(task: tuple[str, str, str, list[tuple[str, int, int, int]], int]) -> dict:
    """Count continuous-block crossings for all boundaries in one BAM."""
    import pysam

    sample, replicate, bam_name, records, min_anchor = task
    # Per contig: sorted coordinate -> [(row index, donor/acceptor), ...].
    by_contig: dict[str, dict[int, list[tuple[int, int]]]] = {}
    for contig, start, end, row_id in records:
        coord_map = by_contig.setdefault(contig, {})
        coord_map.setdefault(start, []).append((row_id, 0))
        coord_map.setdefault(end, []).append((row_id, 1))
    indexes = {
        contig: (sorted(coord_map), coord_map)
        for contig, coord_map in by_contig.items()
    }
    donor: dict[int, int] = {}
    acceptor: dict[int, int] = {}
    examined = primary = crossing_alignments = 0
    bam_path = Path(bam_name)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=False):
            examined += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            primary += 1
            idx = indexes.get(read.reference_name)
            if idx is None:
                continue
            coords, coord_map = idx
            crossed: set[tuple[int, int]] = set()
            for block_start, block_end in read.get_blocks():
                # Inclusive coordinate range satisfying both anchor inequalities.
                lo = bisect.bisect_left(coords, block_start + min_anchor)
                hi = bisect.bisect_right(coords, block_end - min_anchor)
                for boundary in coords[lo:hi]:
                    crossed.update(coord_map[boundary])
            if crossed:
                crossing_alignments += 1
            for row_id, side in crossed:
                target = donor if side == 0 else acceptor
                target[row_id] = target.get(row_id, 0) + 1
    return {
        "sample": sample,
        "replicate": replicate,
        "bam": str(bam_path),
        "alignments_examined": examined,
        "primary_alignments_examined": primary,
        "alignments_crossing_at_least_one_boundary": crossing_alignments,
        "donor": donor,
        "acceptor": acceptor,
    }


def add_junction_signal(
    loci: pd.DataFrame, cfg: argparse.Namespace
) -> tuple[pd.DataFrame, dict[tuple[str, str], np.ndarray]]:
    """Pool <=2-bp expected-junction signal, matching the original scorer."""
    signal = np.zeros(len(loci), dtype=np.int64)
    replicate_count = np.zeros(len(loci), dtype=np.int64)
    by_replicate = {
        (sample, rep): np.zeros(len(loci), dtype=np.int64)
        for sample in SAMPLES for rep in REPLICATES[sample]
    }
    for sample in SAMPLES:
        junctions = parse_junctions([
            (rep, cfg.splice_dir / "regtools" / "per_replicate"
             / f"{rep}.junctions.bed")
            for rep in REPLICATES[sample]
        ])
        subset = loci[(loci["sample"] == sample) & loci.valid_boundaries]
        for row_id, rec in subset.iterrows():
            start, end = int(rec.start), int(rec.end)
            evidence = [
                item
                for ds in range(-2, 3)
                for de in range(-2, 3)
                for item in junctions.get(rec.contig, {}).get((start + ds, end + de), [])
            ]
            signal[row_id] = sum(score for score, _ in evidence)
            replicate_count[row_id] = len({rep for score, rep in evidence if score > 0})
            for score, rep in evidence:
                by_replicate[(sample, rep)][row_id] += score
    out = loci.copy()
    out["junction_signal_within_2bp"] = signal
    out["junction_replicate_count_within_2bp"] = replicate_count
    return out, by_replicate


def summarize(loci: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    for (sample, cls), x in loci.groupby(["sample", "analysis_class"], sort=False):
        covered = x[x.covered_for_psi_flag]
        sample_rows.append({
            "sample": sample, "analysis_class": cls, "n_loci": len(x),
            "valid_boundary_loci": int(x.valid_boundaries.sum()),
            "covered_loci": len(covered),
            "covered_fraction": len(covered) / len(x) if len(x) else math.nan,
            "mean_donor_crossing": x.donor_crossing.sum() / len(x) if len(x) else math.nan,
            "mean_acceptor_crossing": x.acceptor_crossing.sum() / len(x) if len(x) else math.nan,
            "mean_retained_support": x.retained_support.mean(),
            "median_retained_support": x.retained_support.median(),
            "mean_junction_signal": x.junction_signal_exact.mean(),
            "median_junction_signal": x.junction_signal_exact.median(),
            "mean_retention_psi_covered": covered.retention_psi.mean(),
            "median_retention_psi_covered": covered.retention_psi.median(),
            "low_retention_psi_fraction_covered": covered.low_retention_psi_flag.mean(),
            "n_informative_ge_5": int(x.informative_ge_5.sum()),
            "mean_psi_informative_ge_5": x.loc[x.informative_ge_5, "retention_psi"].mean(),
            "median_psi_informative_ge_5": x.loc[x.informative_ge_5, "retention_psi"].median(),
            "n_informative_ge_10": int(x.informative_ge_10.sum()),
            "mean_psi_informative_ge_10": x.loc[x.informative_ge_10, "retention_psi"].mean(),
            "median_psi_informative_ge_10": x.loc[x.informative_ge_10, "retention_psi"].median(),
            "median_boundary_agreement_ratio": x.boundary_agreement_ratio.median(),
        })
    pooled_rows = []
    for cls, x in loci.groupby("analysis_class", sort=False):
        covered = x[x.covered_for_psi_flag]
        pooled_rows.append({
            "analysis_class": cls, "n_loci": len(x),
            "valid_boundary_loci": int(x.valid_boundaries.sum()),
            "covered_loci": len(covered),
            "covered_fraction": len(covered) / len(x) if len(x) else math.nan,
            "mean_retention_psi_covered": covered.retention_psi.mean(),
            "median_retention_psi_covered": covered.retention_psi.median(),
            "q25_retention_psi_covered": covered.retention_psi.quantile(.25),
            "q75_retention_psi_covered": covered.retention_psi.quantile(.75),
            "low_retention_psi_fraction_covered": covered.low_retention_psi_flag.mean(),
            "n_informative_ge_5": int(x.informative_ge_5.sum()),
            "mean_psi_informative_ge_5": x.loc[x.informative_ge_5, "retention_psi"].mean(),
            "median_psi_informative_ge_5": x.loc[x.informative_ge_5, "retention_psi"].median(),
            "n_informative_ge_10": int(x.informative_ge_10.sum()),
            "mean_psi_informative_ge_10": x.loc[x.informative_ge_10, "retention_psi"].mean(),
            "median_psi_informative_ge_10": x.loc[x.informative_ge_10, "retention_psi"].median(),
            "median_boundary_agreement_ratio": x.boundary_agreement_ratio.median(),
        })
    return pd.DataFrame(sample_rows), pd.DataFrame(pooled_rows)


def class_tests(loci: pd.DataFrame) -> pd.DataFrame:
    """Simple locus-level contrasts, reported as descriptive validation only."""
    rows = []
    covered = loci[loci.covered_for_psi_flag]
    contrasts = [
        ("fixed_present", "conventional_intron"),
        ("polymorphic", "conventional_intron"),
        ("polymorphic", "fixed_present"),
    ]
    for focal, reference in contrasts:
        a = covered.loc[covered.analysis_class == focal, "retention_psi"]
        b = covered.loc[covered.analysis_class == reference, "retention_psi"]
        test = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else None
        rows.append({
            "focal_class": focal, "reference_class": reference,
            "n_focal": len(a), "n_reference": len(b),
            "mean_focal": a.mean(), "mean_reference": b.mean(),
            "median_focal": a.median(), "median_reference": b.median(),
            "mean_difference": a.mean() - b.mean(),
            "mann_whitney_u": test.statistic if test else math.nan,
            "mann_whitney_p": test.pvalue if test else math.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    cfg = arguments()
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    loci, construction_audit = construct_loci(cfg)
    loci, junction_by_replicate, resolver_audit = resolve_boundaries(loci, cfg)

    valid = loci[loci.valid_boundaries]
    records_by_sample = {
        sample: [
            (r.contig, int(r.start), int(r.end), int(row_id))
            for row_id, r in valid[valid["sample"] == sample].iterrows()
        ]
        for sample in SAMPLES
    }
    tasks = [
        (sample, rep, str(cfg.bam_dir / f"{rep}.sorted.bam"),
         records_by_sample[sample], cfg.min_anchor)
        for sample in SAMPLES for rep in REPLICATES[sample]
    ]
    donor = np.zeros(len(loci), dtype=np.int64)
    acceptor = np.zeros(len(loci), dtype=np.int64)
    scan_audit = []
    replicate_rows = []
    with ProcessPoolExecutor(max_workers=min(cfg.workers, len(tasks))) as pool:
        futures = {pool.submit(scan_one_bam, task): (task[0], task[1]) for task in tasks}
        for future in as_completed(futures):
            result = future.result()
            rep_donor = result.pop("donor")
            rep_acceptor = result.pop("acceptor")
            for row_id, count in rep_donor.items():
                donor[row_id] += count
            for row_id, count in rep_acceptor.items():
                acceptor[row_id] += count
            sample, rep = result["sample"], result["replicate"]
            rep_subset = loci[loci["sample"] == sample]
            rep_s = junction_by_replicate[(sample, rep)]
            for row_id, rec in rep_subset.iterrows():
                u5 = rep_donor.get(int(row_id), 0)
                u3 = rep_acceptor.get(int(row_id), 0)
                retained = (u5 + u3) / 2.0
                spliced = int(rep_s[row_id])
                total = retained + spliced
                max_u = max(u5, u3)
                replicate_rows.append({
                    "locus_id": rec.locus_id, "sample": sample,
                    "replicate": rep, "gene_id": rec.gene_id,
                    "analysis_class": rec.analysis_class, "contig": rec.contig,
                    "resolved_splice_start": rec.start,
                    "resolved_splice_end": rec.end,
                    "nominal_splice_start": rec.nominal_start,
                    "nominal_splice_end": rec.nominal_end,
                    "boundary_source": rec.boundary_source,
                    "boundary_resolved_flag": rec.boundary_resolved_flag,
                    "primary_boundary_eligible": rec.primary_boundary_eligible,
                    "valid_boundaries": rec.valid_boundaries,
                    "U5_donor_crossing": u5,
                    "U3_acceptor_crossing": u3,
                    "S_junction_signal_exact": spliced,
                    "retained_support_mean_U5_U3": retained,
                    "informative_signal": total,
                    "retention_psi": retained / total if total else math.nan,
                    "boundary_agreement_ratio": min(u5, u3) / max_u if max_u else math.nan,
                    "informative_ge_5": total >= 5,
                    "informative_ge_10": total >= 10,
                })
            scan_audit.append(result)
            print(f"finished {result['replicate']}", flush=True)

    loci["donor_crossing"] = donor
    loci["acceptor_crossing"] = acceptor
    loci["retained_support"] = (donor + acceptor) / 2.0
    max_boundary = np.maximum(donor, acceptor)
    agreement = np.full(len(loci), np.nan, dtype=float)
    np.divide(np.minimum(donor, acceptor), max_boundary, out=agreement,
              where=max_boundary > 0)
    loci["boundary_agreement_ratio"] = agreement
    loci["boundary_crossing_absolute_difference"] = np.abs(donor - acceptor)
    loci["total_psi_signal"] = (
        loci.retained_support + loci.junction_signal_exact)
    loci["covered_for_psi_flag"] = (
        loci.valid_boundaries
        & loci.primary_boundary_eligible
        & (loci.total_psi_signal > 0)
    )
    loci["retention_psi"] = np.where(
        loci.covered_for_psi_flag,
        loci.retained_support / loci.total_psi_signal,
        np.nan,
    )
    loci["spliced_fraction"] = 1 - loci.retention_psi
    loci["low_retention_psi_flag"] = np.where(
        loci.covered_for_psi_flag, loci.retention_psi <= .10, np.nan)
    loci["informative_ge_5"] = loci.total_psi_signal >= 5
    loci["informative_ge_10"] = loci.total_psi_signal >= 10

    sample_summary, pooled_summary = summarize(loci)
    tests = class_tests(loci)
    validation = pd.DataFrame([
        {"check": "invalid_boundary_loci_have_no_psi",
         "value": int(loci.loc[~loci.valid_boundaries, "retention_psi"].notna().sum()),
         "expected": 0},
        {"check": "psi_outside_zero_one",
         "value": int(((loci.retention_psi < 0) | (loci.retention_psi > 1)).sum()),
         "expected": 0},
        {"check": "negative_crossing_counts",
         "value": int(((loci.donor_crossing < 0) | (loci.acceptor_crossing < 0)).sum()),
         "expected": 0},
        {"check": "retained_support_formula_mismatches",
         "value": int((loci.retained_support !=
                       (loci.donor_crossing + loci.acceptor_crossing) / 2).sum()),
         "expected": 0},
        {"check": "eligible_introner_resolver_score_below_minimum",
         "value": int(((loci.analysis_class != "conventional_intron")
                       & loci.primary_boundary_eligible
                       & (loci.junction_signal_exact < cfg.min_resolver_score)).sum()),
         "expected": 0},
        {"check": "eligible_introner_resolver_replicates_below_minimum",
         "value": int(((loci.analysis_class != "conventional_intron")
                       & loci.primary_boundary_eligible
                       & (loci.junction_replicate_count_exact
                          < cfg.min_resolver_replicates)).sum()),
         "expected": 0},
        {"check": "coordinate_collisions_marked_primary_eligible",
         "value": int((loci.resolver_collision_flag
                       & loci.primary_boundary_eligible).sum()),
         "expected": 0},
    ])

    loci.to_csv(cfg.outdir / "per_locus_scores.tsv", sep="\t", index=False)
    pd.DataFrame(replicate_rows).sort_values(
        ["sample", "replicate", "analysis_class", "locus_id"]
    ).to_csv(cfg.outdir / "per_locus_replicate_U5_U3_S.tsv", sep="\t", index=False)
    sample_summary.to_csv(cfg.outdir / "sample_class_summary.tsv", sep="\t", index=False)
    pooled_summary.to_csv(cfg.outdir / "pooled_class_summary.tsv", sep="\t", index=False)
    tests.to_csv(cfg.outdir / "class_tests.tsv", sep="\t", index=False)
    construction_audit.to_csv(cfg.outdir / "construction_audit.tsv", sep="\t", index=False)
    resolver_audit.to_csv(cfg.outdir / "boundary_resolver_audit.tsv", sep="\t", index=False)
    pd.DataFrame(scan_audit).sort_values(["sample", "replicate"]).to_csv(
        cfg.outdir / "bam_scan_audit.tsv", sep="\t", index=False)
    validation.to_csv(cfg.outdir / "validation_checks.tsv", sep="\t", index=False)

    with (cfg.outdir / "RESULTS.md").open("w") as out:
        out.write("# Corrected boundary-crossing retention PSI\n\n")
        out.write(
            "Retention is supported only by a continuous alignment block crossing "
            f"a donor or acceptor with at least {cfg.min_anchor} aligned bases on "
            "both sides. Retained support is the mean of the donor- and "
            "acceptor-crossing counts. Introner boundaries are resolved from "
            "exact regtools junction pairs pooled across four replicates. A pair "
            f"must occur in at least {cfg.min_resolver_replicates} replicates, "
            f"have pooled support >= {cfg.min_resolver_score}, and place both "
            f"endpoints within {cfg.boundary_search_bp} bp of a plausible "
            "predicted interval. Junction signal is then counted only at that "
            "exact resolved pair. Conventional boundaries are exact GTF exon "
            "gaps. Unresolved introners are retained for audit but excluded from "
            "primary PSI summaries. The mating-type region is excluded.\n\n"
        )
        out.write("## Pooled summary\n\n```\n")
        out.write(pooled_summary.to_string(index=False))
        out.write("\n```\n\nThe all-signal summary above is an audit, not a defensible "
                  "retention-rate estimate: loci with unspliced boundary reads but "
                  "no support for the expected junction are dominated by unused "
                  "isoforms or uncertain boundaries. Use "
                  "`evidence_guarded_class_summary.tsv` for interpretation.\n\n"
                  "## Descriptive class tests\n\n```\n")
        out.write(tests.to_string(index=False))
        out.write("\n```\n")
    print(f"wrote {cfg.outdir}")


if __name__ == "__main__":
    main()

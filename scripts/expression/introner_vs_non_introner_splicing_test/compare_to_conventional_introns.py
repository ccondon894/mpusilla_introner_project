#!/usr/bin/env python3
"""Compare introner splicing with clean conventional introns in host genes.

This is an isolated manuscript-support analysis.  It does not modify workflow
outputs. Coordinates are 0-based, half-open after parsing the GTF.
"""

from __future__ import annotations

import argparse
import bisect
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
GROUP = {"CCMP1545": "group1", "RCC1614": "group1", "RCC1749": "group2"}
REPLICATES = {
    "CCMP1545": ("834_A1", "834_A2", "834_A3", "834_B1"),
    "RCC1614": ("1614_A1", "1614_A2", "1614_A3", "1614_A4"),
    "RCC1749": ("1749_A1", "1749_A2", "1749_A4", "1749_B1"),
}
MATING = {
    "CCMP1545": ("CCMP1545#0#scaffold_2", 49808, 1730591),
    "RCC1614": ("RCC1614-intronerized#0#2#0", 49808, 1730591),
    "RCC1749": ("RCC1749#0#intronerless_contig_28", 24999, 2148000),
}
ENDPOINTS = {
    "within_2bp_supported_fraction": "all",
    "low_retention_psi_fraction": "covered",
    "high_confidence_efficient_splicing_fraction": "all",
    "median_retention_psi": "covered",
    "median_spliced_fraction": "covered",
    "median_log1p_junction_support": "all",
}


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", type=Path,
                   default=Path("results/genotyping/genotype_matrix.final.tsv"))
    p.add_argument("--annotations", type=Path, default=Path("results/annotations"))
    p.add_argument("--splice-dir", type=Path,
                   default=Path("results/expression/splice_junctions"))
    p.add_argument("--bam-dir", type=Path,
                   default=Path("results/expression/alignments/bams"))
    p.add_argument("--outdir", type=Path,
                   default=Path("results/expression/functional/introner_vs_non_introner_splicing_test/results"))
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--bootstraps", type=int, default=5000)
    p.add_argument("--seed", type=int, default=1545)
    p.add_argument("--junction-only", action="store_true",
                   help="Skip BAM retention scoring; PSI-dependent metrics are NA.")
    return p.parse_args()


def attr(text: str, key: str) -> str:
    m = re.search(rf'(?:^|;\s*){re.escape(key)}\s+"([^"]+)"', text)
    return m.group(1) if m else ""


def overlap(a: int, b: int, c: int, d: int) -> bool:
    return a < d and c < b


def parse_gtf(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return exact unique introns, gene spans, and an extraction audit."""
    tx_exons: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    tx_info: dict[tuple[str, str], tuple[str, str]] = {}
    gene_parts: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[2] not in {"transcript", "exon"}:
                continue
            contig, feature, start, end, strand = f[0], f[2], int(f[3]) - 1, int(f[4]), f[6]
            gene, tx = attr(f[8], "gene_id"), attr(f[8], "transcript_id")
            if not gene:
                continue
            gene_parts[(contig, gene)].append((start, end))
            if feature == "exon" and tx:
                tx_exons[(contig, tx)].append((start, end))
                tx_info[(contig, tx)] = (gene, strand)

    raw = []
    for (contig, tx), exons in tx_exons.items():
        gene, strand = tx_info[(contig, tx)]
        # Merge overlapping exons before taking gaps, avoiding bogus negative gaps.
        merged = []
        for start, end in sorted(set(exons)):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        for left, right in zip(merged, merged[1:]):
            if left[1] < right[0]:
                raw.append((contig, left[1], right[0], gene, strand, tx))

    # Collapse exact coordinates across isoforms; ambiguous gene assignments are removed.
    by_coord: dict[tuple[str, int, int], list[tuple[str, str, str]]] = defaultdict(list)
    for contig, start, end, gene, strand, tx in raw:
        by_coord[(contig, start, end)].append((gene, strand, tx))
    introns = []
    ambiguous = 0
    for (contig, start, end), values in by_coord.items():
        genes = {x[0] for x in values}
        if len(genes) != 1:
            ambiguous += 1
            continue
        gene = next(iter(genes))
        introns.append({"contig": contig, "start": start, "end": end,
                        "gene_id": gene, "strand": values[0][1],
                        "n_transcripts": len({x[2] for x in values})})
    genes = [{"contig": c, "gene_id": g,
              "start": min(x[0] for x in parts), "end": max(x[1] for x in parts)}
             for (c, g), parts in gene_parts.items()]
    audit = {"gtf_transcripts": len(tx_exons), "raw_transcript_introns": len(raw),
             "exact_unique_introns": len(introns), "ambiguous_gene_coordinates": ambiguous,
             "gtf_genes": len(genes)}
    return pd.DataFrame(introns), pd.DataFrame(genes), audit


class IntervalIndex:
    def __init__(self, rows: list[tuple[int, int, str]]):
        self.rows = sorted(rows)
        self.starts = [x[0] for x in self.rows]

    def overlaps(self, start: int, end: int) -> bool:
        upto = bisect.bisect_left(self.starts, end)
        return any(x[1] > start for x in self.rows[:upto])

    def containing_genes(self, start: int, end: int) -> list[tuple[str, int]]:
        upto = bisect.bisect_right(self.starts, start)
        return [(gene, b - a) for a, b, gene in self.rows[:upto] if end <= b]


def classify(row: pd.Series) -> str:
    group = GROUP[row["sample"]]
    pattern = str(row[f"{group}_pattern"])
    present = int(row[f"{group}_present_count"])
    callable_n = int(row[f"{group}_callable_count"])
    total_n = int(row[f"{group}_n_samples"])
    if pattern == "fixed_present":
        return "fixed_present"
    if callable_n == total_n and 0 < present < total_n:
        return "polymorphic"
    return "excluded"


def present_bodies(matrix: pd.DataFrame, sample: str) -> pd.DataFrame:
    x = matrix[(matrix["sample"] == sample) & (matrix["presence"] == 1)].copy()
    x["body_start"] = x.start.astype(int) + 100
    x["body_end"] = x.end.astype(int) - 100
    return x[x["body_end"] > x["body_start"]]


def host_gene(body: pd.Series, gene_indexes: dict[str, IntervalIndex]) -> str:
    hits = gene_indexes.get(body.contig, IntervalIndex([])).containing_genes(
        int(body.body_start), int(body.body_end))
    if not hits:
        return ""
    # Smallest containing annotation span is the most specific model.
    hits.sort(key=lambda x: (x[1], x[0]))
    return hits[0][0]


def parse_junctions(paths: list[tuple[str, Path]]) -> dict[str, dict[tuple[int, int], list[tuple[int, str]]]]:
    """Index junction evidence by exact boundary pair for O(1) local lookup."""
    out = defaultdict(lambda: defaultdict(list))
    for rep, path in paths:
        with path.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 12:
                    continue
                blocks = [int(x) for x in f[10].rstrip(",").split(",")]
                if len(blocks) < 2:
                    continue
                start, end = int(f[1]) + blocks[0], int(f[2]) - blocks[1]
                if start < end:
                    out[f[0]][(start, end)].append((int(float(f[4])), rep))
    return out


def retained_counts_one_pass(
    controls: pd.DataFrame,
    sample: str,
    cfg: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray]:
    """Count retained reads for all controls while scanning each BAM once.

    This is equivalent to fetching every interval separately and counting a
    read when its aligned blocks contribute at least 8 bp inside that intron.
    A read can count for more than one overlapping control, as it would in the
    per-interval implementation.
    """
    import pysam

    controls = controls.reset_index(drop=True)
    retained = np.zeros(len(controls), dtype=np.int64)
    aligned_bases = np.zeros(len(controls), dtype=np.int64)

    indexes = {}
    for contig, group in controls.groupby("contig", sort=False):
        ordered = group.sort_values(["start", "end"])
        starts = ordered["start"].astype(int).to_numpy()
        ends = ordered["end"].astype(int).to_numpy()
        locus_ids = ordered.index.to_numpy(dtype=int)
        indexes[contig] = (starts, ends, np.maximum.accumulate(ends), locus_ids)

    for rep in REPLICATES[sample]:
        bam_path = cfg.bam_dir / f"{rep}.sorted.bam"
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            references = set(bam.references)
            for contig, (starts, ends, prefix_max_end, locus_ids) in indexes.items():
                if contig not in references:
                    continue
                for read in bam.fetch(contig):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    blocks = read.get_blocks()
                    if not blocks:
                        continue

                    candidates: set[int] = set()
                    for block_start, block_end in blocks:
                        j = bisect.bisect_left(starts, block_end) - 1
                        while j >= 0 and prefix_max_end[j] > block_start:
                            if ends[j] > block_start:
                                candidates.add(int(locus_ids[j]))
                            j -= 1

                    for locus_id in candidates:
                        start = int(controls.at[locus_id, "start"])
                        end = int(controls.at[locus_id, "end"])
                        n_aligned = sum(
                            max(0, min(block_end, end) - max(block_start, start))
                            for block_start, block_end in blocks
                            if block_end > start and block_start < end
                        )
                        if n_aligned >= 8:
                            retained[locus_id] += 1
                            aligned_bases[locus_id] += n_aligned
        print(f"{sample}: finished retained-read scan for {rep}", flush=True)

    return retained, aligned_bases


def score_controls(controls: pd.DataFrame, sample: str, cfg: argparse.Namespace) -> pd.DataFrame:
    controls = controls.reset_index(drop=True)
    reps = REPLICATES[sample]
    junctions = parse_junctions([
        (r, cfg.splice_dir / "regtools" / "per_replicate" / f"{r}.junctions.bed")
        for r in reps
    ])
    if cfg.junction_only:
        retained_by_locus = np.zeros(len(controls), dtype=np.int64)
        aligned_by_locus = np.zeros(len(controls), dtype=np.int64)
    else:
        retained_by_locus, aligned_by_locus = retained_counts_one_pass(
            controls, sample, cfg
        )
    rows = []
    for locus_id, rec in enumerate(controls.itertuples(index=False)):
        nearby = junctions.get(rec.contig, {})
        exactish = [evidence
                    for ds in range(-2, 3) for de in range(-2, 3)
                    for evidence in nearby.get((rec.start + ds, rec.end + de), [])]
        jsignal = sum(score for score, _ in exactish)
        jreps = len({rep for score, rep in exactish if score > 0})
        retained = int(retained_by_locus[locus_id])
        aligned_bp = int(aligned_by_locus[locus_id])
        total = retained + jsignal
        psi = retained / total if total and not cfg.junction_only else math.nan
        rows.append({
                "locus_id": f"{sample}|{rec.contig}:{rec.start}-{rec.end}",
                "sample": sample, "gene_id": rec.gene_id, "contig": rec.contig,
                "body_start": rec.start, "body_end": rec.end, "body_len": rec.end-rec.start,
                "analysis_class": "conventional_intron", "n_transcripts": rec.n_transcripts,
                "retained_signal": math.nan if cfg.junction_only else retained,
                "retained_aligned_bases": math.nan if cfg.junction_only else aligned_bp,
                "psi_junction_signal_within_2bp": jsignal,
                "within_2bp_junction_support": jsignal,
                "within_2bp_supported": jsignal > 0,
                "covered_for_psi_flag": bool(total > 0 and not cfg.junction_only),
                "retention_psi": psi,
                "spliced_fraction": 1-psi if not math.isnan(psi) else math.nan,
                "low_retention_psi_flag": (
                    math.nan if cfg.junction_only else bool(not math.isnan(psi) and psi <= 0.10)),
                "high_confidence_efficient_splicing_flag": (
                    math.nan if cfg.junction_only else bool(
                        not math.isnan(psi) and jsignal > 0 and jreps >= 2 and psi <= 0.10)),
        })
    return pd.DataFrame(rows)


def load_introner_scores(sample: str, bodies: pd.DataFrame, cfg: argparse.Namespace) -> pd.DataFrame:
    p = cfg.splice_dir / "introner_boundary_support" / f"{sample}.per_locus.tsv"
    scores = pd.read_csv(p, sep="\t")
    cols = ["ortholog_id", "gene_id", "analysis_class"]
    merged = scores.merge(bodies[cols], on="ortholog_id", validate="1:1")
    out = pd.DataFrame({
        "locus_id": merged.ortholog_id, "sample": sample, "gene_id": merged.gene_id,
        "contig": merged.contig, "body_start": merged.body_start,
        "body_end": merged.body_end, "body_len": merged.body_len,
        "analysis_class": merged.analysis_class, "n_transcripts": np.nan,
        "retained_signal": pd.to_numeric(merged.retained_signal),
        "retained_aligned_bases": pd.to_numeric(merged.retained_aligned_bases),
        "psi_junction_signal_within_2bp": pd.to_numeric(merged.psi_junction_signal_within_2bp),
        "within_2bp_junction_support": pd.to_numeric(merged.within_2bp_junction_support),
        "within_2bp_supported": pd.to_numeric(merged.within_2bp_junction_support) > 0,
        "covered_for_psi_flag": merged.covered_for_psi.eq("yes"),
        "retention_psi": pd.to_numeric(merged.retention_psi, errors="coerce"),
        "spliced_fraction": pd.to_numeric(merged.spliced_fraction, errors="coerce"),
        "low_retention_psi_flag": merged.low_retention_psi.eq("yes"),
        "high_confidence_efficient_splicing_flag": merged.high_confidence_efficient_splicing.eq("yes"),
    })
    return out


def summaries(loci: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    genes = []
    for (sample, cls), x in loci.groupby(["sample", "analysis_class"]):
        covered = x[x.covered_for_psi_flag]
        rows.append({"sample": sample, "analysis_class": cls, "n_loci": len(x),
                     "n_genes": x.gene_id.nunique(), "covered_loci": len(covered),
                     "within_2bp_supported_fraction": x.within_2bp_supported.mean(),
                     "median_junction_support": x.within_2bp_junction_support.median(),
                     "median_retention_psi": covered.retention_psi.median(),
                     "median_spliced_fraction": covered.spliced_fraction.median(),
                     "low_retention_psi_fraction": covered.low_retention_psi_flag.mean(),
                     "high_confidence_efficient_splicing_fraction": x.high_confidence_efficient_splicing_flag.mean()})
    for (sample, gene, cls), x in loci.groupby(["sample", "gene_id", "analysis_class"]):
        covered = x[x.covered_for_psi_flag]
        genes.append({"sample": sample, "gene_id": gene, "analysis_class": cls,
                      "n_loci": len(x), "n_covered": len(covered),
                      "within_2bp_supported_fraction": x.within_2bp_supported.mean(),
                      "low_retention_psi_fraction": covered.low_retention_psi_flag.mean(),
                      "high_confidence_efficient_splicing_fraction": x.high_confidence_efficient_splicing_flag.mean(),
                      "median_retention_psi": covered.retention_psi.median(),
                      "median_spliced_fraction": covered.spliced_fraction.median(),
                      "median_log1p_junction_support": np.log1p(x.within_2bp_junction_support).median()})
    return pd.DataFrame(rows), pd.DataFrame(genes)


def clustered_test(pairs: pd.DataFrame, rng: np.random.Generator,
                   nperm: int, nboot: int) -> dict:
    pairs = pairs.dropna(subset=["difference"])
    if pairs.empty:
        return {"n_sample_gene_pairs": 0}
    gene_diff = pairs.groupby("gene_id").difference.mean()
    obs = pairs.difference.mean()
    null = np.empty(nperm)
    gd = pairs.gene_id.to_numpy()
    dv = pairs.difference.to_numpy()
    genes = gene_diff.index.to_numpy()
    for i in range(nperm):
        signs = dict(zip(genes, rng.choice((-1, 1), size=len(genes))))
        null[i] = np.mean([d * signs[g] for d, g in zip(dv, gd)])
    boots = np.empty(nboot)
    by_gene = {g: pairs[pairs.gene_id == g].difference.to_numpy() for g in genes}
    for i in range(nboot):
        selected = rng.choice(genes, size=len(genes), replace=True)
        boots[i] = np.concatenate([by_gene[g] for g in selected]).mean()
    try:
        w = stats.wilcoxon(pairs.difference, alternative="two-sided")
        wp = float(w.pvalue)
    except ValueError:
        wp = math.nan
    return {"n_sample_gene_pairs": len(pairs), "n_gene_clusters": len(genes),
            "mean_focal_minus_control": obs, "median_focal_minus_control": pairs.difference.median(),
            "cluster_permutation_p_two_sided": (np.sum(np.abs(null) >= abs(obs)) + 1)/(nperm+1),
            "cluster_bootstrap_ci_low": np.quantile(boots, .025),
            "cluster_bootstrap_ci_high": np.quantile(boots, .975),
            "paired_wilcoxon_p_unclustered": wp}


def tests(gene_summary: pd.DataFrame, cfg: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(cfg.seed)
    test_rows, pair_rows = [], []
    for focal in ("fixed_present", "polymorphic"):
        for metric in ENDPOINTS:
            z = gene_summary[gene_summary.analysis_class.isin([focal, "conventional_intron"])]
            wide = z.pivot(index=["sample", "gene_id"], columns="analysis_class", values=metric).dropna()
            if focal in wide and "conventional_intron" in wide:
                wide = wide.reset_index()
                wide["difference"] = wide[focal] - wide.conventional_intron
            else:
                wide = pd.DataFrame(columns=["sample", "gene_id", focal, "conventional_intron", "difference"])
            wide["focal_class"], wide["endpoint"] = focal, metric
            pair_rows.append(wide)
            for scope, q in [("all_samples", wide)] + [(s, wide[wide["sample"] == s]) for s in SAMPLES]:
                test_rows.append({"scope": scope, "focal_class": focal, "control_class": "conventional_intron",
                                  "endpoint": metric, **clustered_test(q, rng, cfg.permutations, cfg.bootstraps)})
    return pd.DataFrame(test_rows), pd.concat(pair_rows, ignore_index=True)


def main() -> None:
    cfg = args()
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(cfg.matrix, sep="\t", low_memory=False)
    for col in ["start", "end", "presence"]:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")
    # Exclude orthologs anchored in the CCMP1545 mating interval in every sample.
    mt = matrix[(matrix["sample"] == "CCMP1545") & (matrix.contig == MATING["CCMP1545"][0])
                & (matrix.start < MATING["CCMP1545"][2]) & (matrix.end > MATING["CCMP1545"][1])]
    mt_orthologs = set(mt.ortholog_id)
    all_loci, audits, candidates = [], [], []

    for sample in SAMPLES:
        introns, genes, ga = parse_gtf(cfg.annotations / f"{sample}.gtf")
        gene_idx = {c: IntervalIndex([(int(r.start), int(r.end), r.gene_id)
                                     for r in x.itertuples()]) for c, x in genes.groupby("contig")}
        bodies = present_bodies(matrix, sample)
        bodies = bodies[~bodies.ortholog_id.isin(mt_orthologs)].copy()
        bodies["analysis_class"] = bodies.apply(classify, axis=1)
        bodies["gene_id"] = bodies.apply(lambda r: host_gene(r, gene_idx), axis=1)
        primary = bodies[(bodies.analysis_class != "excluded") & bodies.gene_id.ne("")].copy()
        host_genes = set(primary.gene_id)
        body_idx = {c: IntervalIndex([(int(r.body_start), int(r.body_end), "")
                                     for r in x.itertuples()]) for c, x in bodies.groupby("contig")}
        mcontig, mstart, mend = MATING[sample]
        introns["in_host_gene"] = introns.gene_id.isin(host_genes)
        introns["overlaps_present_introner"] = introns.apply(
            lambda r: body_idx.get(r.contig, IntervalIndex([])).overlaps(int(r.start), int(r.end)), axis=1)
        introns["mating_region"] = introns.apply(
            lambda r: r.contig == mcontig and overlap(int(r.start), int(r.end), mstart, mend), axis=1)
        introns["length_compatible"] = (introns.end-introns.start).between(20, 500000)
        clean = introns[introns.in_host_gene & ~introns.overlaps_present_introner
                        & ~introns.mating_region & introns.length_compatible].copy()
        introns.insert(0, "sample", sample)
        introns["selected_clean_control"] = introns.index.isin(clean.index)
        candidates.append(introns)
        controls = score_controls(clean, sample, cfg)
        focal = load_introner_scores(sample, primary, cfg)
        all_loci.extend([controls, focal])
        audits.append({"sample": sample, **ga, "present_introner_bodies": len(bodies),
                       "primary_introner_loci_with_host_gene": len(primary),
                       "host_genes": len(host_genes), "clean_conventional_introns": len(clean),
                       "control_genes": clean.gene_id.nunique(),
                       "removed_nonhost": int((~introns.in_host_gene).sum()),
                       "removed_body_overlap": int(introns.overlaps_present_introner.sum()),
                       "removed_mating": int(introns.mating_region.sum()),
                       "removed_length": int((~introns.length_compatible).sum())})
        print(f"{sample}: {len(primary)} focal introners, {len(clean)} controls in {clean.gene_id.nunique()} genes")

    loci = pd.concat(all_loci, ignore_index=True)
    sample_summary, gene_summary = summaries(loci)
    test_table, pair_table = tests(gene_summary, cfg)
    if cfg.junction_only:
        keep = {"within_2bp_supported_fraction", "median_log1p_junction_support"}
        test_table = test_table[test_table.endpoint.isin(keep)].copy()
        pair_table = pair_table[pair_table.endpoint.isin(keep)].copy()
    pd.concat(candidates, ignore_index=True).to_csv(cfg.outdir/"control_candidate_audit.tsv", sep="\t", index=False)
    pd.DataFrame(audits).to_csv(cfg.outdir/"construction_audit.tsv", sep="\t", index=False)
    loci.to_csv(cfg.outdir/"per_locus_scores.tsv", sep="\t", index=False)
    sample_summary.to_csv(cfg.outdir/"sample_class_summary.tsv", sep="\t", index=False)
    gene_summary.to_csv(cfg.outdir/"sample_gene_class_summary.tsv", sep="\t", index=False)
    pair_table.to_csv(cfg.outdir/"paired_gene_differences.tsv", sep="\t", index=False)
    test_table.to_csv(cfg.outdir/"gene_cluster_aware_tests.tsv", sep="\t", index=False)

    with (cfg.outdir/"RESULTS.md").open("w") as out:
        out.write("# Introner versus conventional-intron splicing\n\n")
        out.write("Controls are exact exon gaps from each base GTF, collapsed across isoforms, restricted to genes with a scored fixed or polymorphic introner, and filtered for any present-introner-body or mating-region overlap. Junction support exactly matches the workflow scorer (pooled replicates and <=2-bp boundary tolerance). ")
        if cfg.junction_only:
            out.write("This fast pass is junction-only: conventional-intron retention/PSI and high-confidence calls are deliberately NA and require a separate BAM/read-depth-matched pass.\n\n")
        else:
            out.write("Retention uses >=8 aligned bp, PSI <=0.10, and >=2 replicates for high confidence.\n\n")
        out.write("## Locus summaries\n\n```\n" + sample_summary.to_string(index=False) + "\n```\n\n")
        key = test_table[(test_table.scope == "all_samples") & test_table.endpoint.isin(["within_2bp_supported_fraction", "median_retention_psi", "median_spliced_fraction", "high_confidence_efficient_splicing_fraction"])]
        out.write("## Paired sample-gene tests\n\nPositive differences mean the introner class has a larger value than conventional introns in the same sample and host gene. Inference resamples/sign-flips homologous gene IDs as clusters.\n\n```\n" + key.to_string(index=False) + "\n```\n")
    print(f"Wrote {cfg.outdir}")


if __name__ == "__main__":
    main()

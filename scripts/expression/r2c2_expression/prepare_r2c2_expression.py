#!/usr/bin/env python3
"""R2C2 replication of the existing expression pilot, with matched short-read fits.

Uses the original fitting functions; R2C2 exposure is size factor only. Gene
features and short-read observations come from the verified September 9 pilot.
"""
from __future__ import annotations

import os
os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")
os.environ.setdefault("TMPDIR", "/scratch1/chris/tmp")

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
ORIGINAL = ROOT / "scripts/expression/negative_binomial_expression_pilot"
sys.path.insert(0, str(ORIGINAL))
import fit_expression_pilot_models as fit
from prepare_expression_pilot_data import median_ratio_size_factors, reference_gene_table

SAMPLES = {"834": "CCMP1545", "1614": "RCC1614", "1749": "RCC1749"}
RCC1749_MT_START = 25000
RCC1749_MT_END = 2148000
FOCAL = [
    ("paired_turnover_gene_fe_ppml", "reference_relative_gain_count", "Gain, within gene"),
    ("paired_turnover_gene_fe_ppml", "reference_relative_loss_count", "Loss, within gene"),
    ("current_state_gene_fe_ppml", "current_introner_count", "Current burden, within gene (FE)"),
    ("current_state_nb_gee", "current_introner_count_within", "Current burden, within gene (NB-GEE)"),
    ("current_state_nb_gee", "current_introner_count_between", "Current burden, between genes"),
]


def normalize(data, long_reads):
    data = data.copy()
    factors, n = median_ratio_size_factors(data)
    normalization = data[["replicate", "strain"]].drop_duplicates().set_index("replicate")
    normalization["median_ratio_size_factor"] = factors
    normalization["library_size"] = data.groupby("replicate").raw_count.sum()
    normalization["library_size_factor"] = normalization.library_size / np.exp(
        np.log(normalization.library_size).mean())
    for col in normalization.columns.drop("strain"):
        data[col] = data.replicate.map(normalization[col])
    exposure_length = 1.0 if long_reads else data.effective_exon_length / 1000
    data["offset_median_ratio"] = np.log(data.median_ratio_size_factor * exposure_length)
    data["offset_library"] = np.log(data.library_size_factor * exposure_length)
    data["normalized_expression"] = data.raw_count / np.exp(data.offset_median_ratio)
    # Required by the original log-expression sensitivity helper. For R2C2 this
    # is log normalized counts, NOT counts per kb; only this adapter uses the name.
    data["log1p_normalized_count_per_kb"] = np.log1p(data.normalized_expression)
    return data, normalization.reset_index(), n


def prepare(quant_dir, out, source):
    base = pd.read_csv(source, sep="\t")
    paths = [source, Path(__file__), ORIGINAL / "fit_expression_pilot_models.py",
             ORIGINAL / "prepare_expression_pilot_data.py"]
    mt = set()
    for strain in ["CCMP1545", "RCC1749"]:
        path = ROOT / f"results/annotations/{strain}.gtf"
        paths.append(path)
        genes = reference_gene_table(path)
        mask = (genes.contig.eq("CCMP1545#0#scaffold_2")) if strain == "CCMP1545" else (
                    genes.contig.eq("RCC1749#0#intronerless_contig_28"))
        mt.update(genes.loc[mask, "gene_id"])

    base = base[~base.gene_id.isin(mt) & ~base.mating_type_gene].copy()
    stale = ["raw_count", "median_ratio_size_factor", "library_size", "library_size_factor",
             "effective_exon_kb", "offset_median_ratio_length", "offset_library_length",
             "normalized_count_per_kb", "log1p_normalized_count_per_kb"]
    covariates = base.drop(columns=stale)
    frames, audits, crosswalks = [], [], []
    for label, strain in SAMPLES.items():
        qp = quant_dir / f"09092025_{label}_Isoforms.filtered.clean.quant"
        sp = ROOT / f"data/sqanti3_output_{label}/{label}_isoforms_classification.filtered.txt"
        paths.extend([qp, sp])
        q = pd.read_csv(qp, sep="\t")
        s = pd.read_csv(sp, sep="\t")
        cols = [c for c in q if c.startswith(label + "_")]
        assert len(cols) == 4 and not q.Isoform.duplicated().any()
        assert set(q.Isoform) == set(s.isoform)
        assert np.isfinite(q[cols]).all(axis=None) and (q[cols] >= 0).all(axis=None)
        assert (q[cols] == np.floor(q[cols])).all(axis=None)
        z = q.merge(s[["isoform", "associated_gene", "chrom", "structural_category", *cols]],
                    left_on="Isoform", right_on="isoform", validate="one_to_one",
                    suffixes=("", "_sqanti"))
        assert all((z[c] == z[c + "_sqanti"]).all() for c in cols)
        z["gene_id"] = z.associated_gene.str.strip()
        eligible = set(covariates.loc[covariates.strain.eq(strain), "gene_id"])
        z["retained"] = z.gene_id.isin(eligible)
        z["strain"] = strain
        z["status"] = np.where(z.retained, "retained_exact_current_gene_id",
                               np.where(z.gene_id.isin(mt), "mating_type", "no_eligible_pilot_gene"))
        crosswalks.append(z[["strain", "Isoform", "Gene", "gene_id", "chrom",
                            "structural_category", "retained", "status"]])
        kept = z.loc[z.retained]
        grouped = kept.groupby("gene_id")[cols].sum()
        assert (grouped.sum() == kept[cols].sum()).all()
        long = grouped.reset_index().melt(id_vars="gene_id", var_name="replicate", value_name="raw_count")
        long["strain"] = strain
        frames.append(long)
        audits.append(dict(strain=strain, input_isoforms=len(q), retained_isoforms=len(kept),
                           retained_genes=len(grouped), input_counts=int(q[cols].sum().sum()),
                           retained_counts=int(kept[cols].sum().sum()),
                           counts_identical_to_sqanti=True))
    counts = pd.concat(frames, ignore_index=True)
    data = counts.merge(covariates, on=["gene_id", "strain", "replicate"], validate="one_to_one")
    assert len(data) == len(counts)
    assert not data.gene_id.isin(mt).any()
    assert not data.duplicated(["gene_id", "replicate"]).any()
    assert data.groupby(["gene_id", "strain"]).size().eq(4).all()
    short = base.merge(data[["gene_id", "strain", "replicate"]],
                      on=["gene_id", "strain", "replicate"], validate="one_to_one")
    long, norm, n = normalize(data, True)
    short, short_norm, short_n = normalize(short, False)
    pd.concat(crosswalks).to_csv(out / "isoform_gene_crosswalk.tsv", sep="\t", index=False)
    pd.DataFrame(audits).to_csv(out / "mapping_audit.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": sorted(mt)}).to_csv(out / "excluded_mating_type_gene_ids.tsv", sep="\t", index=False)
    for name, frame, normalization in [("r2c2", long, norm), ("short_read_matched", short, short_norm)]:
        folder = out / name
        folder.mkdir(exist_ok=True)
        frame.to_csv(folder / "expression_model_data.tsv", sep="\t", index=False)
        normalization.to_csv(folder / "normalization_factors.tsv", sep="\t", index=False)
    provenance = [{"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                   "bytes": p.stat().st_size} for p in paths]
    (out / "provenance.json").write_text(json.dumps(dict(inputs=provenance,
        normalization_genes_r2c2=n, normalization_genes_short_read_matched=short_n), indent=2))
    return long, norm, short, short_norm



def main():
    parser = argparse.ArgumentParser(description="Prepare gene-level R2C2 expression observations")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    prepare(args.quant_dir, args.outdir, args.input)


if __name__ == "__main__":
    main()

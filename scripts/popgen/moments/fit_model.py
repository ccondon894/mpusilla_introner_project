#!/usr/bin/env python3
"""moments split-migration demographic fitting for two populations.

This mirrors the dadi workflow in workflow/rules/23_demography.smk:
  1. reads the 4D SNP VCF plus popinfo file,
  2. builds the observed folded 2D SFS for intronerful vs intronerless,
  3. fits a split-migration model with multiple starts,
  4. generates block-bootstrap SFS replicates,
  5. fits the same model to each bootstrap replicate.

The defaults match workflow/rules/23_demography.smk for ad hoc runs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from multiprocessing import Pool
import random
import sys
from pathlib import Path

import moments
import numpy as np

BASE_DEMO_MODEL = moments.Demographics2D.split_mig
INITIAL_PARAMS = [0.7, 0.7, 12, 0.005]
LOWER_BOUNDS = [1e-2, 1e-2, 1e-3, 1e-4]
UPPER_BOUNDS = [3, 3, 30, 1]


def split_mig_projected(params, ns, folded=False, pop_ids=None):
    """
    Evaluate split_mig, projecting from n=3 when a 2D axis has only n=2.

    moments 1.5.3's 2D integration fails for an axis with sample size 2. The
    observed spectrum can still be n=2; we integrate at n=3 for that axis and
    project the model spectrum back down before likelihood evaluation.
    """
    stable_ns = [max(int(n), 3) for n in ns]
    fs = BASE_DEMO_MODEL(params, stable_ns, pop_ids=pop_ids)
    if stable_ns != list(ns):
        fs = fs.project(ns)
        fs.pop_ids = pop_ids
    if folded:
        fs = fs.fold()
    fs[fs < 0] = 0
    return fs


DEMO_MODEL = split_mig_projected


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["replicate", "path", "segregating_sites", "shape", "folded"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def ensure_folded_state(fs, folded: bool):
    if folded and not fs.folded:
        return fs.fold()
    if not folded and fs.folded:
        raise ValueError("Cannot use a folded spectrum for a polarized/unfolded fit.")
    return fs


def write_fit_results(path: Path, results: list[tuple[float, list[float], float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for ll_model, params, theta in results:
            row = [ll_model] + list(params) + [theta]
            handle.write("\t".join(str(value) for value in row) + "\n")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_popinfo(path: Path) -> dict[str, str]:
    sample_to_pop: dict[str, str] = {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed popinfo line {line_no}: {line.rstrip()!r}")
            sample, pop_id = parts[:2]
            sample_to_pop[sample] = pop_id
    return sample_to_pop


def info_value(info: str, key: str) -> str | None:
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            return field[len(prefix) :]
    return None


def count_gt_alleles(gt: str) -> tuple[int, int] | None:
    if not gt or gt == ".":
        return None

    sep = "/" if "/" in gt else "|" if "|" in gt else None
    alleles = gt.split(sep) if sep else [gt]
    if any(allele == "." for allele in alleles):
        return None

    ref_count = 0
    alt_count = 0
    for allele in alleles:
        if allele == "0":
            ref_count += 1
        elif allele == "1":
            alt_count += 1
        else:
            return None
    return ref_count, alt_count


def make_haploid_aware_data_dict_vcf(
    vcf: Path, popinfo: Path, pop_ids: list[str]
) -> dict[tuple[str, int], dict[str, object]]:
    """
    Parse biallelic SNP GT calls into the data-dict format used by moments.

    moments.Misc.make_data_dict_vcf assumes diploid GT strings of length 3.
    This project VCF is haploid-coded, so we count alleles explicitly.
    """
    sample_to_pop = read_popinfo(popinfo)
    wanted_pops = set(pop_ids)
    sample_pop_indices: list[tuple[int, str]] = []
    data_dict: dict[tuple[str, int], dict[str, object]] = {}

    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                sample_pop_indices = [
                    (idx, sample_to_pop[sample])
                    for idx, sample in enumerate(samples)
                    if sample_to_pop.get(sample) in wanted_pops
                ]
                if not sample_pop_indices:
                    raise ValueError("No VCF samples matched the requested pop IDs.")
                matched_pops = {pop_id for _idx, pop_id in sample_pop_indices}
                missing_pops = [pop_id for pop_id in pop_ids if pop_id not in matched_pops]
                if missing_pops:
                    raise ValueError(
                        "No VCF samples matched requested pop IDs: "
                        + ", ".join(missing_pops)
                    )
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos_s, _snp_id, ref, alt, _qual, filt, info, fmt = fields[:9]
            if filt not in {"PASS", "."}:
                continue
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            if ref.upper() not in {"A", "C", "G", "T"}:
                continue
            if alt.upper() not in {"A", "C", "G", "T"}:
                continue

            fmt_parts = fmt.split(":")
            if "GT" not in fmt_parts:
                continue
            gt_index = fmt_parts.index("GT")

            calls = {pop_id: [0, 0] for pop_id in pop_ids}
            for sample_idx, pop_id in sample_pop_indices:
                sample_field = fields[9 + sample_idx]
                sample_parts = sample_field.split(":")
                if gt_index >= len(sample_parts):
                    continue
                allele_counts = count_gt_alleles(sample_parts[gt_index])
                if allele_counts is None:
                    continue
                calls[pop_id][0] += allele_counts[0]
                calls[pop_id][1] += allele_counts[1]

            if not any(sum(pop_calls) for pop_calls in calls.values()):
                continue

            aa = info_value(info, "AA")
            outgroup_allele = aa.upper()[0] if aa and aa not in {".", "-"} else "-"
            data_dict[(chrom, int(pos_s))] = {
                "segregating": [ref.upper(), alt.upper()],
                "calls": {
                    pop_id: tuple(pop_calls) for pop_id, pop_calls in calls.items()
                },
                "outgroup_allele": outgroup_allele,
            }

    return data_dict


def site_coordinate_from_key(key: object) -> tuple[str, int]:
    """Return (chromosome, 1-based position) from moments data-dict keys."""
    if isinstance(key, tuple) and len(key) >= 2:
        chrom, pos = key[0], key[1]
        return str(chrom), int(pos)

    if isinstance(key, str):
        chrom, sep, pos = key.rpartition("_")
        if sep and pos.isdigit():
            return chrom, int(pos)

    raise ValueError(
        "Could not infer chromosome/position from a moments data-dict key. "
        f"Unexpected key format: {key!r}"
    )


def write_chunk_bed(data_dict: dict[object, object], chunk_size: int, path: Path) -> int:
    """
    Write non-overlapping BED intervals for SNP-bearing physical chunks.

    moments.Misc.bootstrap accepts a BED file defining resampling units. VCF
    positions are 1-based, while BED starts are 0-based and ends are half-open.
    """
    if chunk_size <= 0:
        raise ValueError(f"--chunk-size must be positive, got {chunk_size}")

    chunks: set[tuple[str, int, int]] = set()
    for key in data_dict:
        chrom, pos = site_coordinate_from_key(key)
        start = ((pos - 1) // chunk_size) * chunk_size
        end = start + chunk_size
        chunks.add((chrom, start, end))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for chrom, start, end in sorted(chunks, key=lambda item: (item[0], item[1])):
            handle.write(f"{chrom}\t{start}\t{end}\n")

    return len(chunks)


def fit_single(
    args_tuple: tuple[object, list[float], int, int, int, bool]
) -> tuple[float, list[float], float]:
    fs, starting_params, perturb_fold, maxiter, seed, folded = args_tuple
    random.seed(seed)
    np.random.seed(seed)

    p0 = moments.Misc.perturb_params(
        starting_params,
        fold=perturb_fold,
        lower_bound=LOWER_BOUNDS,
        upper_bound=UPPER_BOUNDS,
    )
    popt = moments.Inference.optimize_log(
        p0,
        fs,
        DEMO_MODEL,
        lower_bound=LOWER_BOUNDS,
        upper_bound=UPPER_BOUNDS,
        maxiter=maxiter,
        verbose=0,
        func_kwargs={"folded": folded},
    )
    model_fs = DEMO_MODEL(popt, fs.sample_sizes, folded=folded)
    ll_model = moments.Inference.ll_multinom(model_fs, fs)
    theta = moments.Inference.optimal_sfs_scaling(model_fs, fs)
    return float(ll_model), [float(value) for value in popt], float(theta)


def run_parallel_fits(
    spectra: list[object],
    starting_params: list[float],
    perturb_fold: int,
    maxiter: int,
    threads: int,
    label: str,
    folded: bool,
) -> list[tuple[float, list[float], float]]:
    base_seed = random.randint(0, 2**31 - 1)
    job_args = [
        (fs, starting_params, perturb_fold, maxiter, base_seed + idx, folded)
        for idx, fs in enumerate(spectra)
    ]

    results: list[tuple[float, list[float], float]] = []
    total = len(job_args)
    if threads > 1:
        with Pool(processes=threads) as pool:
            for result in pool.imap_unordered(fit_single, job_args):
                results.append(result)
                print(
                    f"  {label} {len(results)}/{total}: LL={result[0]:.2f}",
                    file=sys.stderr,
                    flush=True,
                )
    else:
        for idx, job_arg in enumerate(job_args, start=1):
            result = fit_single(job_arg)
            results.append(result)
            print(
                f"  {label} {idx}/{total}: LL={result[0]:.2f}",
                file=sys.stderr,
                flush=True,
            )
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate moments bootstrap 2D SFS files from a VCF/popinfo pair."
    )
    parser.add_argument(
        "--vcf",
        default="results/snp_popgen/vcf/mpusilla.snps.4d.notMT.vcf.gz",
        help="Input VCF/VCF.gz file. Default: %(default)s",
    )
    parser.add_argument(
        "--popinfo",
        default="results/snp_popgen/demography/popinfo.txt",
        help="Two-column sample-to-population file. Default: %(default)s",
    )
    parser.add_argument(
        "--outdir",
        default="analysis/moments_pilot",
        help="Output directory. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Observed-data fit table. Default: OUTDIR/model_fits.4d.txt",
    )
    parser.add_argument(
        "--bootstrap-output",
        default=None,
        help="Bootstrap fit table. Default: OUTDIR/model_fits.4d.bootstrap.txt",
    )
    parser.add_argument(
        "--pop-ids",
        nargs=2,
        default=["intronerful", "intronerless"],
        metavar=("POP1", "POP2"),
        help="Population IDs in popinfo order for the 2D SFS. Default: %(default)s",
    )
    parser.add_argument(
        "--projections",
        nargs=2,
        type=int,
        default=[11, 2],
        metavar=("N1", "N2"),
        help="Haploid projections/sample sizes. Default matches dadi rule: %(default)s",
    )
    parser.add_argument(
        "--n-boot",
        type=int,
        default=100,
        help="Number of bootstrap replicate spectra to fit. Default: %(default)s",
    )
    parser.add_argument(
        "--n-opt",
        type=int,
        default=20,
        help="Number of observed-data optimization starts. Default: %(default)s",
    )
    parser.add_argument(
        "--maxiter",
        type=int,
        default=10000,
        help="Maximum moments optimizer iterations per fit. Default: %(default)s",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of parallel worker processes for fitting. Default: %(default)s",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=250000,
        help="Physical block size for fragmenting the data dictionary. Default: %(default)s",
    )
    parser.add_argument(
        "--polarized",
        action="store_true",
        help="Treat SNPs as polarized. By default spectra are folded/unpolarized, matching the dadi rule.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_boot < 0:
        raise SystemExit(f"--n-boot must be non-negative, got {args.n_boot}")
    if args.n_opt < 1:
        raise SystemExit(f"--n-opt must be at least 1, got {args.n_opt}")
    if args.maxiter < 1:
        raise SystemExit(f"--maxiter must be at least 1, got {args.maxiter}")
    if args.threads < 1:
        raise SystemExit(f"--threads must be at least 1, got {args.threads}")

    vcf = Path(args.vcf)
    popinfo = Path(args.popinfo)
    outdir = Path(args.outdir)
    bootdir = outdir / "bootstraps"
    output_path = Path(args.output) if args.output else outdir / "model_fits.4d.txt"
    bootstrap_output_path = (
        Path(args.bootstrap_output)
        if args.bootstrap_output
        else outdir / "model_fits.4d.bootstrap.txt"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    bootdir.mkdir(parents=True, exist_ok=True)

    if not vcf.exists():
        raise SystemExit(f"VCF not found: {vcf}")
    if not popinfo.exists():
        raise SystemExit(f"popinfo not found: {popinfo}")

    polarized = bool(args.polarized)
    folded = not polarized

    print(f"Loading VCF: {vcf}", file=sys.stderr)
    print(f"Loading popinfo: {popinfo}", file=sys.stderr)
    dd = make_haploid_aware_data_dict_vcf(vcf, popinfo, args.pop_ids)

    print(
        f"Building observed SFS for pop_ids={args.pop_ids}, projections={args.projections}, "
        f"polarized={polarized}",
        file=sys.stderr,
    )
    observed_fs = moments.Spectrum.from_data_dict(
        dd, args.pop_ids, args.projections, polarized=polarized
    )
    observed_fs = ensure_folded_state(observed_fs, folded)
    observed_path = outdir / "observed.2d.fs"
    observed_fs.to_file(str(observed_path))

    print(
        f"Fitting observed SFS with {args.n_opt} starts, maxiter={args.maxiter}, "
        f"threads={args.threads}",
        file=sys.stderr,
    )
    opt_results = run_parallel_fits(
        [observed_fs for _ in range(args.n_opt)],
        INITIAL_PARAMS,
        perturb_fold=1,
        maxiter=args.maxiter,
        threads=args.threads,
        label="Opt",
        folded=folded,
    )
    best_ll, best_params, best_theta = max(opt_results, key=lambda result: result[0])
    write_fit_results(output_path, opt_results)
    print(
        "Best observed fit: "
        f"LL={best_ll:.4f}, N1={best_params[0]:.4f}, N2={best_params[1]:.4f}, "
        f"T={best_params[2]:.4f}, M={best_params[3]:.6f}, theta={best_theta:.2f}",
        file=sys.stderr,
    )
    print(f"Observed fit table: {output_path}", file=sys.stderr)

    chunk_bed = outdir / f"bootstrap_chunks.{args.chunk_size}.bed"
    print(
        f"Writing {args.chunk_size:,} bp bootstrap regions to {chunk_bed}",
        file=sys.stderr,
    )
    n_chunks = write_chunk_bed(dd, args.chunk_size, chunk_bed)
    if n_chunks == 0:
        raise SystemExit("No SNP-bearing chunks were found in the data dictionary.")

    print(
        f"Generating {args.n_boot} bootstrap SFS replicates from {n_chunks} chunks",
        file=sys.stderr,
    )
    bootstraps = moments.Misc.bootstrap(
        dd,
        args.pop_ids,
        args.projections,
        mask_corners=True,
        polarized=polarized,
        bed_filename=str(chunk_bed),
        num_boots=args.n_boot,
        save_dir=None,
    )

    rows: list[dict[str, object]] = [
        {
            "replicate": "observed",
            "path": observed_path,
            "segregating_sites": float(observed_fs.S()),
            "shape": "x".join(map(str, observed_fs.shape)),
            "folded": folded,
        }
    ]

    for idx, boot_fs in enumerate(bootstraps, start=1):
        boot_fs = ensure_folded_state(boot_fs, folded)
        bootstraps[idx - 1] = boot_fs
        boot_path = bootdir / f"bootstrap_{idx:04d}.2d.fs"
        boot_fs.to_file(str(boot_path))
        rows.append(
            {
                "replicate": idx,
                "path": boot_path,
                "segregating_sites": float(boot_fs.S()),
                "shape": "x".join(map(str, boot_fs.shape)),
                "folded": folded,
            }
        )
        print(f"  wrote {idx}/{args.n_boot}: {boot_path}", file=sys.stderr)

    if bootstraps:
        print(
            f"Fitting {len(bootstraps)} bootstrap spectra from the best observed fit",
            file=sys.stderr,
        )
        boot_results = run_parallel_fits(
            bootstraps,
            best_params,
            perturb_fold=2,
            maxiter=args.maxiter,
            threads=args.threads,
            label="Bootstrap",
            folded=folded,
        )
        write_fit_results(bootstrap_output_path, boot_results)
        print(f"Bootstrap fit table: {bootstrap_output_path}", file=sys.stderr)

    summary_path = outdir / "bootstrap_summary.tsv"
    write_summary(summary_path, rows)

    print("Done.", file=sys.stderr)
    print(f"Observed SFS: {observed_path}", file=sys.stderr)
    print(f"Observed fit table: {output_path}", file=sys.stderr)
    print(f"Bootstrap SFS directory: {bootdir}", file=sys.stderr)
    print(f"Bootstrap fit table: {bootstrap_output_path}", file=sys.stderr)
    print(f"Summary: {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

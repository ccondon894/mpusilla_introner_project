"""
DADI Model Fitting - Version 2 (Improved Convergence)

Fits a split-migration demographic model using dadi.

Improvements over v1:
1. Multi-start optimization to avoid local optima
2. Best-fit seeding for bootstrap runs
3. Parallel execution of independent optimization runs
4. Accepts pre-computed SFS or VCF input
5. CLI interface for Snakemake integration
"""

import argparse
import sys
import numpy as np
from multiprocessing import Pool
from functools import partial

import nlopt
import dadi


# ============================================================================
# Model setup
# ============================================================================

DEMO_MODEL = dadi.Demographics2D.split_mig
LOWER_BOUNDS = [1e-2, 1e-2, 1e-3, 1e-4]
UPPER_BOUNDS = [3, 3, 30, 1]
INITIAL_PARAMS = [0.7, 0.7, 12, 0.005]


def fit_single(args_tuple):
    """
    Fit demographic model once. Designed for use with multiprocessing.Pool.

    Args:
        args_tuple: (fs, starting_params, pts_l, ns, fold, maxeval, seed)

    Returns:
        (log_likelihood, optimized_params, theta)
    """
    fs, starting_params, pts_l, ns, fold, maxeval, seed = args_tuple

    # Each process needs its own extrapolation function
    demo_model_ex = dadi.Numerics.make_extrap_func(DEMO_MODEL)

    # Set random seed for reproducibility across parallel runs
    np.random.seed(seed)

    p0 = dadi.Misc.perturb_params(starting_params, fold=fold,
                                  upper_bound=UPPER_BOUNDS,
                                  lower_bound=LOWER_BOUNDS)

    popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
                                        lower_bound=LOWER_BOUNDS,
                                        upper_bound=UPPER_BOUNDS,
                                        algorithm=nlopt.LN_BOBYQA,
                                        maxeval=maxeval,
                                        verbose=0)

    model_fs = demo_model_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

    return ll_model, popt, theta0


def write_results(results, output_file):
    """Write list of (ll, popt, theta) results to TSV file."""
    with open(output_file, 'w') as f:
        for ll_model, popt, theta0 in results:
            res = [ll_model] + list(popt) + [theta0]
            f.write('\t'.join([str(x) for x in res]) + '\n')


def run_parallel_fits(fs, starting_params, pts_l, ns, n_runs, fold,
                      maxeval, threads, label="Fitting"):
    """
    Run multiple independent model fits in parallel.

    Returns:
        List of (log_likelihood, optimized_params, theta) tuples.
    """
    # Build argument tuples with unique seeds
    base_seed = np.random.randint(0, 2**31)
    job_args = [
        (fs, starting_params, pts_l, ns, fold, maxeval, base_seed + i)
        for i in range(n_runs)
    ]

    results = []
    if threads > 1:
        with Pool(processes=threads) as pool:
            for i, result in enumerate(pool.imap_unordered(fit_single, job_args)):
                ll = result[0]
                results.append(result)
                print(f"  {label} {len(results)}/{n_runs}: LL={ll:.2f}", flush=True)
    else:
        for i, args in enumerate(job_args):
            result = fit_single(args)
            results.append(result)
            print(f"  {label} {i+1}/{n_runs}: LL={result[0]:.2f}", flush=True)

    return results


def load_data(args):
    """
    Load frequency spectrum from SFS file or VCF + popinfo.

    Returns:
        data_fs: dadi Spectrum object
        dd: data dictionary (for bootstrapping), or None if loaded from SFS
        pop_ids: population IDs
        ns: sample sizes
    """
    pop_ids = ["intronerful", "intronerless"]
    ns = [11, 2]

    if args.vcf:
        print(f"Loading data from VCF: {args.vcf}")
        dd = dadi.Misc.make_data_dict_vcf(args.vcf, args.popinfo)
        data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)
    elif args.sfs:
        print(f"Loading SFS from: {args.sfs}")
        data_fs = dadi.Spectrum.from_file(args.sfs)
        dd = None
        if args.popinfo:
            dd = dadi.Misc.make_data_dict_vcf(args.vcf, args.popinfo) if args.vcf else None
    else:
        sys.exit("Error: must provide either --sfs or --vcf")

    print(f"  Spectrum shape: {data_fs.shape}")
    print(f"  Total SNPs: {data_fs.S():.0f}")

    return data_fs, dd, pop_ids, ns


def main():
    parser = argparse.ArgumentParser(
        description="Fit split-migration demographic model using dadi")
    # Input (one of these required)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--sfs", help="Pre-computed dadi SFS file (.fs)")
    input_group.add_argument("--vcf", help="Input VCF file")
    parser.add_argument("--popinfo", help="Population info file (required with --vcf)")

    # Output
    parser.add_argument("--output", required=True,
                        help="Output file for multi-start optimization results")
    parser.add_argument("--bootstrap", help="Output file for bootstrap results")

    # Run parameters
    parser.add_argument("--n_opt", type=int, default=20,
                        help="Number of multi-start optimizations (default: 20)")
    parser.add_argument("--n_boot", type=int, default=100,
                        help="Number of bootstrap replicates (default: 100)")
    parser.add_argument("--maxeval", type=int, default=10000,
                        help="Max iterations per optimization (default: 10000)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of parallel threads")

    args = parser.parse_args()

    if args.vcf and not args.popinfo:
        parser.error("--popinfo is required when using --vcf")

    # Load data
    data_fs, dd, pop_ids, ns = load_data(args)
    pts_l = [max(ns) + 20, max(ns) + 30, max(ns) + 40]

    # ========================================================================
    # STAGE 1: Multi-start optimization on real data
    # ========================================================================
    print()
    print("=" * 70)
    print(f"STAGE 1: Multi-start optimization ({args.n_opt} runs, {args.threads} threads)")
    print("=" * 70)

    opt_results = run_parallel_fits(
        data_fs, INITIAL_PARAMS, pts_l, ns,
        n_runs=args.n_opt, fold=1, maxeval=args.maxeval,
        threads=args.threads, label="Opt"
    )

    # Find best fit
    best_idx = np.argmax([r[0] for r in opt_results])
    best_ll, best_params, best_theta = opt_results[best_idx]

    print()
    print("Best fit from multi-start optimization:")
    print(f"  Log-likelihood: {best_ll:.4f}")
    print(f"  N1 (intronerful):  {best_params[0]:.4f}")
    print(f"  N2 (intronerless): {best_params[1]:.4f}")
    print(f"  T (split time):    {best_params[2]:.4f}")
    print(f"  M (migration):     {best_params[3]:.6f}")
    print(f"  Theta:             {best_theta:.2f}")

    # Write multi-start results
    write_results(opt_results, args.output)
    print(f"\nOptimization results written to: {args.output}")

    # ========================================================================
    # STAGE 2: Bootstrap analysis
    # ========================================================================
    if args.bootstrap and dd is not None:
        print()
        print("=" * 70)
        print(f"STAGE 2: Bootstrap analysis ({args.n_boot} replicates, {args.threads} threads)")
        print("=" * 70)
        print(f"Seeding bootstraps from best-fit: N1={best_params[0]:.3f}, "
              f"N2={best_params[1]:.3f}, T={best_params[2]:.3f}, M={best_params[3]:.5f}")
        print()

        # Generate bootstrap SFS replicates
        chunks = dadi.Misc.fragment_data_dict(dd, 250000)
        boots = dadi.Misc.bootstraps_from_dd_chunks(
            chunks, args.n_boot, pop_ids, ns, polarized=False)

        # Fit each bootstrap replicate in parallel
        base_seed = np.random.randint(0, 2**31)
        boot_args = [
            (boot_fs, best_params, pts_l, ns, 2, args.maxeval, base_seed + i)
            for i, boot_fs in enumerate(boots)
        ]

        boot_results = []
        if args.threads > 1:
            with Pool(processes=args.threads) as pool:
                for i, result in enumerate(pool.imap_unordered(fit_single, boot_args)):
                    boot_results.append(result)
                    print(f"  Bootstrap {len(boot_results)}/{args.n_boot}: "
                          f"LL={result[0]:.2f}", flush=True)
        else:
            for i, ba in enumerate(boot_args):
                result = fit_single(ba)
                boot_results.append(result)
                print(f"  Bootstrap {i+1}/{args.n_boot}: LL={result[0]:.2f}", flush=True)

        write_results(boot_results, args.bootstrap)
        print(f"\nBootstrap results written to: {args.bootstrap}")

    elif args.bootstrap and dd is None:
        print("\nSkipping bootstrap: requires --vcf input for data dict generation")

    # ========================================================================
    # Summary
    # ========================================================================
    print()
    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()

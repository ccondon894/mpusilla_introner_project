"""Fit single-population demographic models to the folded 4D SFS.

Fits constant-size, two-epoch, and three-epoch models in moments,
selects by AIC, and pickles the best model's parameters along with a
diagnostic plot showing observed vs expected.
"""

import argparse
import pickle
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import moments


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sfs", required=True, help="Path to unfolded SFS .npy")
    p.add_argument("--aic", required=True, help="Output AIC table TSV")
    p.add_argument("--best", required=True, help="Pickle of best-fit model")
    p.add_argument("--plot", required=True, help="Diagnostic PNG")
    p.add_argument(
        "--residuals",
        required=True,
        help="Per-bin residuals TSV for the best model",
    )
    p.add_argument(
        "--restarts",
        type=int,
        default=5,
        help="Random restarts per non-trivial model (default 5)",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---- Model definitions -----------------------------------------------------

# Each entry: (name, function(params, ns), p0, lower, upper)
# Constant-size has no free demographic params (just theta).
MODELS = {
    "constant": {
        "func": moments.Demographics1D.snm,
        "p0": [],
        "lower": [],
        "upper": [],
    },
    "two_epoch": {
        "func": moments.Demographics1D.two_epoch,
        # params: [nu, T] — relative size, time in 2N generations.
        "p0": [1.0, 0.1],
        "lower": [1e-3, 1e-4],
        "upper": [100.0, 5.0],
    },
    "three_epoch": {
        "func": moments.Demographics1D.three_epoch,
        # params: [nuB, nuF, TB, TF]
        "p0": [0.5, 1.5, 0.05, 0.05],
        "lower": [1e-3, 1e-3, 1e-4, 1e-4],
        "upper": [100.0, 100.0, 5.0, 5.0],
    },
    "growth": {
        "func": moments.Demographics1D.growth,
        # params: [nu, T] — exponential growth from ancestral size to nu over T.
        "p0": [2.0, 0.1],
        "lower": [1e-3, 1e-4],
        "upper": [100.0, 5.0],
    },
    "bottlegrowth": {
        "func": moments.Demographics1D.bottlegrowth,
        # params: [nuB, nuF, T] — instantaneous bottleneck nuB, then exp.
        # growth to nuF over T.
        "p0": [0.1, 2.0, 0.1],
        "lower": [1e-3, 1e-3, 1e-4],
        "upper": [100.0, 100.0, 5.0],
    },
}


def _build_model_array(spec, params, n):
    """Build the folded model as a plain numpy array (no Spectrum scalar mul).

    Avoids a numpy >= 1.25 / moments 1.5.x __array_wrap__ incompatibility.
    """
    if len(params) == 0:
        spectrum = spec["func"]([n])
    else:
        spectrum = spec["func"](list(params), [n])
    folded = spectrum.fold()
    arr = np.ma.filled(folded, 0.0).astype(float)
    mask = folded.mask if hasattr(folded, "mask") else np.zeros_like(arr, dtype=bool)
    return arr, np.asarray(mask, dtype=bool)


def _ll_multinom(model_arr, mask, data_arr):
    """Multinomial log-likelihood (constant terms dropped)."""
    use = (~mask) & (model_arr > 0) & (data_arr > 0)
    model_norm = model_arr[use] / model_arr[use].sum()
    return float(np.sum(data_arr[use] * np.log(model_norm)))


def _optimal_theta(model_arr, mask, data_arr):
    use = ~mask
    return float(data_arr[use].sum() / model_arr[use].sum())


def summary_stats_unfolded(unfolded):
    """Compute pi, theta_W, Tajima's D from an unfolded SFS (length n+1).

    Uses bins i = 1..n-1 (segregating). Returns dict.
    """
    n = len(unfolded) - 1
    seg = unfolded[1:n].astype(float)
    S = float(seg.sum())
    if S == 0 or n < 2:
        return {"pi": 0.0, "theta_W": 0.0, "tajima_D": float("nan"),
                "n_segregating": 0, "n_samples": n}
    i = np.arange(1, n)
    # pi (per-site, summed) via i*(n-i)/(n choose 2) weighting
    denom_pi = n * (n - 1) / 2.0
    pi = float(np.sum(i * (n - i) * seg) / denom_pi)
    # Watterson's theta
    a1 = float(np.sum(1.0 / i))
    a2 = float(np.sum(1.0 / (i ** 2)))
    theta_W = S / a1
    # Tajima's D variance (Tajima 1989)
    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = 2.0 * (n ** 2 + n + 3) / (9.0 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)
    var_D = e1 * S + e2 * S * (S - 1)
    if var_D <= 0:
        D = float("nan")
    else:
        D = (pi - theta_W) / np.sqrt(var_D)
    return {
        "pi": pi,
        "theta_W": theta_W,
        "tajima_D": D,
        "n_segregating": int(S),
        "n_samples": n,
    }


def summary_stats_folded(folded_arr, n):
    """Compute pi from a folded SFS of length n+1; theta_W from total S.

    For folded data we cannot compute Tajima's D the same way without
    polarization assumptions, but pi and theta_W are well-defined.
    """
    seg = folded_arr[1:n]
    # pi: each folded bin i covers minor-count i (and n-i which folded into it).
    # The unfolded-pi weight for an SFS bin at frequency i is i*(n-i)/C(n,2).
    # When we fold, the bin at folded index i has total mass = unfolded[i] +
    # unfolded[n-i] for i < n/2 (and just unfolded[n/2] when n is even). The
    # weight i*(n-i)/C(n,2) is the SAME for unfolded[i] and unfolded[n-i], so
    # folded pi = sum_{i=1}^{floor(n/2)} (i*(n-i)/C(n,2)) * folded[i].
    half = n // 2
    denom_pi = n * (n - 1) / 2.0
    i_idx = np.arange(1, half + 1)
    pi = float(np.sum(i_idx * (n - i_idx) * folded_arr[1: half + 1]) / denom_pi)
    S = float(seg.sum())
    a1 = float(np.sum(1.0 / np.arange(1, n)))
    theta_W = S / a1 if a1 > 0 else 0.0
    return {"pi_folded": pi, "theta_W_folded": theta_W, "n_segregating_folded": int(S)}


def fit_model(name, spec, data_arr, data_mask, n, rng, restarts):
    """Optimize a model with random restarts. Return best (params, theta, ll, aic)."""
    p0 = np.array(spec["p0"], dtype=float)
    k = len(p0)

    def neg_ll(log_params):
        params = np.exp(log_params)
        if np.any(params < spec["lower"]) or np.any(params > spec["upper"]):
            return 1e12
        try:
            model_arr, mask = _build_model_array(spec, params, n)
        except Exception:
            return 1e12
        joint_mask = data_mask | mask
        return -_ll_multinom(model_arr, joint_mask, data_arr)

    if k == 0:
        model_arr, mask = _build_model_array(spec, [], n)
        joint_mask = data_mask | mask
        ll = _ll_multinom(model_arr, joint_mask, data_arr)
        theta = _optimal_theta(model_arr, joint_mask, data_arr)
        aic = 2 * 1 - 2 * ll
        return {"params": [], "theta": float(theta), "ll": float(ll), "aic": float(aic)}

    best = None
    for r in range(restarts):
        if r == 0:
            p_init = p0
        else:
            # Multiplicative perturbation in log space, clipped to bounds.
            log_lo = np.log(spec["lower"])
            log_hi = np.log(spec["upper"])
            log_p0 = np.log(p0)
            p_init = np.exp(
                np.clip(log_p0 + rng.uniform(-1.5, 1.5, size=k), log_lo, log_hi)
            )
        try:
            res = minimize(
                neg_ll,
                np.log(p_init),
                method="Nelder-Mead",
                options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 500},
            )
        except Exception as exc:  # pragma: no cover
            print(f"  [{name}] restart {r} failed: {exc}", file=sys.stderr)
            continue
        if not np.isfinite(res.fun):
            continue
        params = np.exp(res.x)
        model_arr, mask = _build_model_array(spec, params, n)
        joint_mask = data_mask | mask
        ll = _ll_multinom(model_arr, joint_mask, data_arr)
        theta = _optimal_theta(model_arr, joint_mask, data_arr)
        aic = 2 * (k + 1) - 2 * ll
        if best is None or ll > best["ll"]:
            best = {
                "params": [float(x) for x in params],
                "theta": float(theta),
                "ll": float(ll),
                "aic": float(aic),
            }

    if best is None:
        raise RuntimeError(f"All restarts failed for model {name}")
    return best


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    unfolded = np.load(args.sfs)
    n = len(unfolded) - 1
    data_fs = moments.Spectrum(unfolded.astype(float)).fold()
    data_arr = np.ma.filled(data_fs, 0.0).astype(float)
    data_mask = data_fs.mask if hasattr(data_fs, "mask") else np.zeros_like(data_arr, dtype=bool)
    data_mask = np.asarray(data_mask, dtype=bool)

    print(f"Loaded SFS: n={n}, total sites={int(unfolded.sum())}", file=sys.stderr)

    obs_unfolded_stats = summary_stats_unfolded(unfolded)
    obs_folded_stats = summary_stats_folded(data_arr, n)
    print(
        "Observed (unfolded) pi={pi:.2f} theta_W={theta_W:.2f} D={tajima_D:.4f} "
        "S={n_segregating}".format(**obs_unfolded_stats),
        file=sys.stderr,
    )

    fits = {}
    for name, spec in MODELS.items():
        print(f"Fitting {name}...", file=sys.stderr)
        fits[name] = fit_model(name, spec, data_arr, data_mask, n, rng, args.restarts)
        print(
            f"  ll={fits[name]['ll']:.3f}  aic={fits[name]['aic']:.3f}  "
            f"theta={fits[name]['theta']:.3f}",
            file=sys.stderr,
        )

    # Write AIC table.
    with open(args.aic, "w") as fh:
        fh.write("model\tk_demog\tll\ttheta\taic\tparams\n")
        for name, f in fits.items():
            k = len(f["params"])
            params_str = ",".join(f"{x:.6g}" for x in f["params"])
            fh.write(
                f"{name}\t{k}\t{f['ll']:.6f}\t{f['theta']:.6f}\t{f['aic']:.6f}\t{params_str}\n"
            )
        # Trailing comment-style lines so consumers can grep these out if
        # parsing the table.
        fh.write(
            "# observed_pi\t{pi:.6f}\n".format(**obs_unfolded_stats)
        )
        fh.write(
            "# observed_theta_W\t{theta_W:.6f}\n".format(**obs_unfolded_stats)
        )
        fh.write(
            "# observed_tajima_D\t{tajima_D:.6f}\n".format(**obs_unfolded_stats)
        )
        fh.write(
            "# observed_n_segregating\t{n_segregating}\n".format(**obs_unfolded_stats)
        )
        fh.write("# observed_pi_folded\t{pi_folded:.6f}\n".format(**obs_folded_stats))
        fh.write(
            "# observed_theta_W_folded\t{theta_W_folded:.6f}\n".format(**obs_folded_stats)
        )

    # Best model = lowest AIC.
    best_name = min(fits, key=lambda k: fits[k]["aic"])
    print(f"Best model by AIC: {best_name}", file=sys.stderr)

    best = {
        "model": best_name,
        "func_name": MODELS[best_name]["func"].__name__,
        "params": fits[best_name]["params"],
        "theta": fits[best_name]["theta"],
        "ll": fits[best_name]["ll"],
        "aic": fits[best_name]["aic"],
        "n": n,
        "observed_unfolded_stats": obs_unfolded_stats,
        "observed_folded_stats": obs_folded_stats,
    }
    with open(args.best, "wb") as fh:
        pickle.dump(best, fh)

    # ---- Per-bin residuals for the best model ----------------------------
    half = n // 2
    best_spec = MODELS[best_name]
    best_arr, best_mask = _build_model_array(best_spec, fits[best_name]["params"], n)
    best_expected = best_arr * fits[best_name]["theta"]
    obs_folded = data_arr
    with open(args.residuals, "w") as fh:
        fh.write("model\tbin\tobserved\texpected\tresidual\tpearson_z\n")
        for i in range(1, half + 1):
            o = float(obs_folded[i])
            e = float(best_expected[i])
            r = o - e
            z = r / np.sqrt(e) if e > 0 else float("nan")
            fh.write(f"{best_name}\t{i}\t{o:.6f}\t{e:.6f}\t{r:.6f}\t{z:.6f}\n")

    # ---- Two-panel diagnostic plot ---------------------------------------
    bins = np.arange(1, half + 1)
    obs_seg = obs_folded[1 : half + 1]
    colors = {
        "constant": "tab:blue",
        "two_epoch": "tab:orange",
        "three_epoch": "tab:green",
        "growth": "tab:red",
        "bottlegrowth": "tab:brown",
    }

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.5]},
    )

    ax_top.bar(bins, obs_seg, color="lightgray", edgecolor="black", label="observed")
    for name, f in fits.items():
        spec = MODELS[name]
        model_arr, _ = _build_model_array(spec, f["params"], n)
        model_arr = model_arr * f["theta"]
        ax_top.plot(
            bins,
            model_arr[1 : half + 1],
            "o-",
            color=colors.get(name, "k"),
            label=f"{name} (AIC={f['aic']:.1f})",
        )
    ax_top.set_ylabel("Sites")
    ax_top.set_title("Observed vs expected folded 4D SFS (top); Pearson residuals for best model (bottom)")
    ax_top.legend(fontsize=8)

    # Pearson residuals: (obs - exp) / sqrt(exp), best model only.
    expected_best = best_expected[1 : half + 1]
    pearson = (obs_seg - expected_best) / np.sqrt(np.where(expected_best > 0, expected_best, np.nan))
    ax_bot.bar(bins, pearson, color=colors.get(best_name, "tab:gray"),
               edgecolor="black")
    ax_bot.axhline(0, color="black", lw=0.6)
    ax_bot.axhline(2, color="red", lw=0.5, ls="--")
    ax_bot.axhline(-2, color="red", lw=0.5, ls="--")
    ax_bot.set_xticks(bins)
    ax_bot.set_xlabel("Minor allele count")
    ax_bot.set_ylabel(f"Pearson z\n({best_name})")
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    plt.close(fig)

    print(f"Wrote {args.aic}, {args.best}, {args.plot}, {args.residuals}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

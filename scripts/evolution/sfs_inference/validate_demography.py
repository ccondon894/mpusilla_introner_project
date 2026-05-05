"""Validate the fitted demographic model against summary statistics.

Loads the best-fit demographic model + observed unfolded SFS. Builds the
expected unfolded SFS under the fitted model (theta-scaled), computes
analytic summary statistics (pi, theta_W, Tajima's D), and Poisson-resamples
the expected SFS to build a null distribution of each statistic. Compares
the observed value to the null distribution and saves a histogram + summary
TSV.

Caveat: Poisson resampling captures only the per-bin sampling variance and
underestimates the true coalescent variance. So the simulated distribution
is *narrower* than the true distribution under the fitted demography. If
the observed statistic falls outside the Poisson-resampled distribution the
model is definitely misspecified; if it falls inside, model fit is at most
"adequate at this resolution."
"""

import argparse
import pickle
import sys

import numpy as np
import matplotlib.pyplot as plt
import moments


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sfs", required=True, help="observed unfolded SFS .npy")
    p.add_argument("--demog", required=True, help="demography_best.pkl")
    p.add_argument("--summary", required=True, help="output summary TSV")
    p.add_argument("--plot", required=True, help="output histogram PNG")
    p.add_argument("--n_replicates", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def summary_stats_unfolded(unfolded):
    n = len(unfolded) - 1
    seg = unfolded[1:n].astype(float)
    S = float(seg.sum())
    if S == 0 or n < 2:
        return {"pi": 0.0, "theta_W": 0.0, "tajima_D": float("nan"), "S": 0}
    i = np.arange(1, n)
    denom_pi = n * (n - 1) / 2.0
    pi = float(np.sum(i * (n - i) * seg) / denom_pi)
    a1 = float(np.sum(1.0 / i))
    a2 = float(np.sum(1.0 / (i ** 2)))
    theta_W = S / a1
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
    return {"pi": pi, "theta_W": theta_W, "tajima_D": D, "S": int(S)}


def expected_unfolded(demog):
    n = demog["n"]
    func = getattr(moments.Demographics1D, demog["func_name"])
    if len(demog["params"]) == 0:
        spec = func([n])
    else:
        spec = func(list(demog["params"]), [n])
    arr = np.asarray(spec).astype(float)
    return arr * demog["theta"]


def expected_folded(demog):
    n = demog["n"]
    func = getattr(moments.Demographics1D, demog["func_name"])
    if len(demog["params"]) == 0:
        spec = func([n])
    else:
        spec = func(list(demog["params"]), [n])
    folded = spec.fold()
    arr = np.ma.filled(folded, 0.0).astype(float)
    return arr * demog["theta"]


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    unfolded = np.load(args.sfs).astype(float)
    n = len(unfolded) - 1

    demog = pickle.load(open(args.demog, "rb"))
    if demog["n"] != n:
        raise SystemExit(f"Demog n={demog['n']} but SFS n={n}")

    obs_stats = summary_stats_unfolded(unfolded)
    print(
        f"Observed: pi={obs_stats['pi']:.3f}  theta_W={obs_stats['theta_W']:.3f}  "
        f"D={obs_stats['tajima_D']:.4f}  S={obs_stats['S']}",
        file=sys.stderr,
    )

    # Expected unfolded SFS under the fitted demography. Bins 0 and n are
    # interpreted as "monomorphic ancestral" / "fixed derived" respectively
    # — only bins 1..n-1 contribute to pi/theta_W/D.
    exp_unfolded = expected_unfolded(demog)
    # Force non-negative for Poisson sampling.
    exp_seg = np.clip(exp_unfolded[1:n], 0.0, None)

    # Analytic stats from the *expected* SFS itself (no resampling).
    expected_full = np.zeros(n + 1)
    expected_full[1:n] = exp_seg
    exp_stats = summary_stats_unfolded(expected_full)
    print(
        f"Expected: pi={exp_stats['pi']:.3f}  theta_W={exp_stats['theta_W']:.3f}  "
        f"D={exp_stats['tajima_D']:.4f}  S={exp_stats['S']}",
        file=sys.stderr,
    )

    # Poisson resample replicates.
    rep_pi = np.zeros(args.n_replicates)
    rep_thetaW = np.zeros(args.n_replicates)
    rep_D = np.zeros(args.n_replicates)
    rep_S = np.zeros(args.n_replicates, dtype=int)
    for r in range(args.n_replicates):
        sim_seg = rng.poisson(exp_seg).astype(float)
        sim_full = np.zeros(n + 1)
        sim_full[1:n] = sim_seg
        st = summary_stats_unfolded(sim_full)
        rep_pi[r] = st["pi"]
        rep_thetaW[r] = st["theta_W"]
        rep_D[r] = st["tajima_D"]
        rep_S[r] = st["S"]

    def pct_position(obs, sim):
        """One-sided p (fraction of sims as or more extreme than obs)."""
        finite = sim[np.isfinite(sim)]
        if len(finite) == 0:
            return float("nan")
        # Two-sided: fraction with |sim - mean| >= |obs - mean|
        c = float(np.mean(np.abs(finite - finite.mean()) >= abs(obs - finite.mean())))
        return c

    rows = [
        ("model", demog["model"]),
        ("params", ",".join(f"{x:.6g}" for x in demog["params"])),
        ("theta", f"{demog['theta']:.6g}"),
        ("n_samples", n),
        ("n_replicates", args.n_replicates),
        ("note", "Poisson resampling underestimates coalescent variance"),
    ]
    for stat_name, obs_val, exp_val, rep in [
        ("pi", obs_stats["pi"], exp_stats["pi"], rep_pi),
        ("theta_W", obs_stats["theta_W"], exp_stats["theta_W"], rep_thetaW),
        ("tajima_D", obs_stats["tajima_D"], exp_stats["tajima_D"], rep_D),
        ("S_segregating", obs_stats["S"], exp_stats["S"], rep_S.astype(float)),
    ]:
        finite = rep[np.isfinite(rep)]
        rows.extend([
            (f"{stat_name}_observed", f"{obs_val:.6g}"),
            (f"{stat_name}_expected_analytic", f"{exp_val:.6g}"),
            (f"{stat_name}_sim_mean", f"{float(finite.mean()):.6g}"),
            (f"{stat_name}_sim_sd", f"{float(finite.std()):.6g}"),
            (f"{stat_name}_sim_p2.5", f"{float(np.percentile(finite, 2.5)):.6g}"),
            (f"{stat_name}_sim_p97.5", f"{float(np.percentile(finite, 97.5)):.6g}"),
            (f"{stat_name}_two_sided_p", f"{pct_position(obs_val, rep):.6g}"),
        ])

    with open(args.summary, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")

    # ---- 3-panel histogram plot ------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    panels = [
        ("pi", rep_pi, obs_stats["pi"], exp_stats["pi"]),
        ("theta_W", rep_thetaW, obs_stats["theta_W"], exp_stats["theta_W"]),
        ("tajima_D", rep_D, obs_stats["tajima_D"], exp_stats["tajima_D"]),
    ]
    for ax, (name, rep, obs_val, exp_val) in zip(axes, panels):
        finite = rep[np.isfinite(rep)]
        ax.hist(finite, bins=40, color="lightblue", edgecolor="black",
                label="Poisson sims")
        ax.axvline(obs_val, color="red", lw=2, label=f"observed={obs_val:.3f}")
        ax.axvline(exp_val, color="black", lw=1.5, ls="--",
                   label=f"expected={exp_val:.3f}")
        ax.set_xlabel(name)
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)
        ax.set_title(name)
    fig.suptitle(
        f"Demography validation under {demog['model']} "
        f"(observed vs Poisson-resampled simulations)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    plt.close(fig)

    print(f"Wrote {args.summary}, {args.plot}", file=sys.stderr)


if __name__ == "__main__":
    main()

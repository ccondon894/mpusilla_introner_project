"""Joint MLE of (theta_lambda, theta_mu) on the introner two-state SFS.

The introner SFS (segregating bins i = 1..n-1, indexed by present-allele
count) is modeled as a Poisson Random Field mixture of:
    - gain spectrum: neutral SFS shape under fitted demography, present is
      derived → bin i corresponds directly to present-frequency i.
    - loss spectrum: same neutral SFS, but absent is derived → bin i in
      present-frequency = bin n-i in derived (absent) frequency, i.e. flip.

Both components are scaled by independent thetas:
    E[sfs_gain[i]] = theta_lambda * neutral_unfolded[i]
    E[sfs_loss[i]] = theta_mu     * neutral_unfolded[n - i]
    E[sfs_total[i]] = sum

We fit (theta_lambda, theta_mu) by minimizing the negative Poisson
log-likelihood of observed bin counts.
"""

import argparse
import pickle
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import chi2

import moments


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sfs", required=True, help="introner_sfs.npy (length n+1)")
    p.add_argument("--summary", required=True, help="introner_sfs_summary.tsv")
    p.add_argument("--demog", required=True, help="demography_best.pkl")
    p.add_argument("--estimates", required=True, help="output theta_estimates.tsv")
    p.add_argument("--gof", required=True, help="output gain_loss_gof.tsv")
    p.add_argument("--plot", required=True, help="output gain_loss_fit.png")
    return p.parse_args()


def load_demog(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def neutral_unfolded(demog, n):
    """Unfolded expected SFS (length n+1) at theta=1 under the fitted model."""
    func = getattr(moments.Demographics1D, demog["func_name"])
    if len(demog["params"]) == 0:
        spec = func([n])
    else:
        spec = func(list(demog["params"]), [n])
    # `spec` is a moments.Spectrum. Cast to plain array; segregating bins are
    # 1..n-1; bin 0 and bin n hold "fixation" probability and are typically
    # masked but moments doesn't fold here, so we just slice.
    return np.asarray(spec)


def expected_two_state(theta_lambda, theta_mu, neutral, n):
    """Return (gain_seg, loss_seg, total_seg) restricted to i = 1..n-1."""
    seg = neutral[1:n]               # length n-1, present-derived shape
    gain = theta_lambda * seg
    loss = theta_mu * seg[::-1]      # absent-derived: flip
    return gain, loss, gain + loss


def neg_log_lik(log_thetas, sfs_obs_seg, neutral, n):
    theta_lambda, theta_mu = np.exp(log_thetas)
    _, _, total = expected_two_state(theta_lambda, theta_mu, neutral, n)
    if np.any(total <= 0):
        return np.inf
    # Poisson log pmf summed over bins.
    return float(np.sum(total - sfs_obs_seg * np.log(total)))


def gtest(observed, expected, min_expected=5.0):
    """Pool consecutive bins to reach min_expected, return G stat / df / p."""
    obs = list(observed)
    exp = list(expected)
    pooled_obs = []
    pooled_exp = []
    cur_o = 0.0
    cur_e = 0.0
    for o, e in zip(obs, exp):
        cur_o += o
        cur_e += e
        if cur_e >= min_expected:
            pooled_obs.append(cur_o)
            pooled_exp.append(cur_e)
            cur_o = 0.0
            cur_e = 0.0
    # Tail pool: fold any leftover into the last bin.
    if cur_e > 0:
        if pooled_exp:
            pooled_obs[-1] += cur_o
            pooled_exp[-1] += cur_e
        else:
            pooled_obs.append(cur_o)
            pooled_exp.append(cur_e)

    pooled_obs = np.array(pooled_obs)
    pooled_exp = np.array(pooled_exp)
    # G = 2 * sum(O * log(O / E)) with 0 * log(0) := 0
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(pooled_obs > 0, pooled_obs * np.log(pooled_obs / pooled_exp), 0.0)
    G = 2.0 * float(np.sum(ratio))
    df = max(len(pooled_obs) - 1 - 2, 1)  # subtract 2 fitted params
    pval = float(chi2.sf(G, df))
    return G, df, pval, pooled_obs, pooled_exp


def main():
    args = parse_args()

    sfs_full = np.load(args.sfs)
    n = len(sfs_full) - 1
    sfs_obs_seg = sfs_full[1:n].astype(float)
    print(f"Loaded introner SFS, n={n}, segregating={int(sfs_obs_seg.sum())}",
          file=sys.stderr)

    summary = pd.read_csv(args.summary, sep="\t")
    fixed_present = int(summary.loc[summary["key"] == "fixed_present", "value"].iloc[0])
    fixed_absent = int(summary.loc[summary["key"] == "fixed_absent", "value"].iloc[0])
    n_called = int(summary.loc[summary["key"] == "n_loci_called", "value"].iloc[0])

    demog = load_demog(args.demog)
    if demog["n"] != n:
        raise SystemExit(
            f"Demography fit on n={demog['n']} but introner SFS has n={n}"
        )
    neutral = neutral_unfolded(demog, n)
    print(f"Demography: {demog['model']}  params={demog['params']}", file=sys.stderr)
    print(f"Neutral unfolded[1..n-1]={neutral[1:n].tolist()}", file=sys.stderr)

    # Initial guess: split observed segregating count evenly between gain/loss
    # and back out theta from the segregating mass of `neutral`.
    seg_mass = float(neutral[1:n].sum())
    initial = max(sfs_obs_seg.sum() / 2.0, 1.0) / max(seg_mass, 1e-12)
    x0 = np.log([initial, initial])

    res = minimize(
        neg_log_lik,
        x0,
        args=(sfs_obs_seg, neutral, n),
        method="Nelder-Mead",
        options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 10000},
    )
    theta_lambda_hat, theta_mu_hat = np.exp(res.x)
    print(
        f"MLE: theta_lambda={theta_lambda_hat:.6g}  theta_mu={theta_mu_hat:.6g}  "
        f"ratio={theta_lambda_hat / theta_mu_hat:.4f}",
        file=sys.stderr,
    )

    gain, loss, total = expected_two_state(theta_lambda_hat, theta_mu_hat, neutral, n)
    G, df, pval, pooled_obs, pooled_exp = gtest(sfs_obs_seg, total)

    # Per-bin posterior probabilities of gain vs loss origin.
    # P(gain | i) = expected_gain[i] / expected_total[i] under the fitted mixture.
    safe_total = np.where(total > 0, total, 1.0)
    p_gain = gain / safe_total
    p_loss = loss / safe_total

    # Boundary-mutation sanity: theta per locus should be much less than 1.
    # We don't know an exact "locus count" denominator so report theta /
    # n_called as a rough check.
    bm_lambda = theta_lambda_hat / max(n_called, 1)
    bm_mu = theta_mu_hat / max(n_called, 1)

    fixed_ratio = (fixed_present / fixed_absent) if fixed_absent > 0 else float("inf")
    fitted_ratio = theta_lambda_hat / theta_mu_hat

    # ---- Write theta_estimates.tsv ---------------------------------------
    with open(args.estimates, "w") as fh:
        fh.write("key\tvalue\n")
        fh.write(f"demog_model\t{demog['model']}\n")
        fh.write(f"demog_params\t{','.join(f'{x:.6g}' for x in demog['params'])}\n")
        fh.write(f"n_samples\t{n}\n")
        fh.write(f"n_segregating\t{int(sfs_obs_seg.sum())}\n")
        fh.write(f"theta_lambda\t{theta_lambda_hat:.6g}\n")
        fh.write(f"theta_mu\t{theta_mu_hat:.6g}\n")
        fh.write(f"theta_ratio_lambda_over_mu\t{fitted_ratio:.6g}\n")
        fh.write(f"expected_gain_total\t{float(gain.sum()):.6g}\n")
        fh.write(f"expected_loss_total\t{float(loss.sum()):.6g}\n")
        fh.write(
            "frac_segregating_from_gains\t"
            f"{float(gain.sum() / total.sum()):.6g}\n"
        )
        fh.write(f"fixed_present\t{fixed_present}\n")
        fh.write(f"fixed_absent\t{fixed_absent}\n")
        fh.write(f"empirical_fixed_ratio_P_over_A\t{fixed_ratio:.6g}\n")
        fh.write(
            "fitted_minus_empirical_log_ratio\t"
            f"{(np.log(fitted_ratio) - np.log(fixed_ratio)) if fixed_absent > 0 and fitted_ratio > 0 else 'NA'}\n"
        )
        fh.write(f"boundary_check_lambda_per_locus\t{bm_lambda:.6g}\n")
        fh.write(f"boundary_check_mu_per_locus\t{bm_mu:.6g}\n")
        fh.write(f"optimization_success\t{bool(res.success)}\n")
        fh.write(f"final_neg_loglik\t{res.fun:.6g}\n")
        # Per-bin posterior probabilities (origin given present-allele count).
        for i, (pg, pl) in enumerate(zip(p_gain, p_loss), start=1):
            fh.write(f"P_gain[i={i}]\t{pg:.6g}\n")
            fh.write(f"P_loss[i={i}]\t{pl:.6g}\n")

    # ---- Write goodness-of-fit TSV ---------------------------------------
    with open(args.gof, "w") as fh:
        fh.write(
            "present_count\tobserved\texpected_gain\texpected_loss\texpected_total"
            "\tP_gain\tP_loss\n"
        )
        for i, (o, g, l, t, pg, pl) in enumerate(
            zip(sfs_obs_seg, gain, loss, total, p_gain, p_loss), start=1
        ):
            fh.write(
                f"{i}\t{int(o)}\t{g:.6g}\t{l:.6g}\t{t:.6g}\t{pg:.6g}\t{pl:.6g}\n"
            )
        fh.write("\n# Pooled G-test (pooled to expected >= 5)\n")
        fh.write(f"# G={G:.6f}  df={df}  pval={pval:.6g}\n")
        fh.write("pool_idx\tobserved\texpected\n")
        for i, (po, pe) in enumerate(zip(pooled_obs, pooled_exp)):
            fh.write(f"{i}\t{po:.6g}\t{pe:.6g}\n")

    # ---- Plot --------------------------------------------------------------
    bins = np.arange(1, n)
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(bins - width / 2, sfs_obs_seg, width=width,
           label="observed", color="lightgray", edgecolor="black")
    ax.bar(bins + width / 2, total, width=width,
           label="expected total", color="tab:purple", edgecolor="black")
    ax.plot(bins, gain, "o-", color="tab:green", label="gain component")
    ax.plot(bins, loss, "s-", color="tab:red", label="loss component")
    ax.set_xticks(bins)
    ax.set_xlabel("Present-allele count")
    ax.set_ylabel("Number of polymorphic loci")
    ax.set_title(
        f"Introner SFS gain/loss fit\n"
        f"theta_lambda={theta_lambda_hat:.3g}  theta_mu={theta_mu_hat:.3g}  "
        f"ratio={fitted_ratio:.3f}  G={G:.2f} (p={pval:.3g})"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    plt.close(fig)

    print(
        f"Goodness-of-fit: G={G:.3f} df={df} p={pval:.3g}",
        file=sys.stderr,
    )
    print(
        f"fixed_P/fixed_A={fixed_ratio:.3f}  fitted_lambda/mu={fitted_ratio:.3f}",
        file=sys.stderr,
    )
    print(f"Wrote {args.estimates}, {args.gof}, {args.plot}", file=sys.stderr)


if __name__ == "__main__":
    main()

"""Side-by-side comparison of introner vs non-introner-intron gain/loss fits.

Reads parametric fit outputs from the standardized pipeline layout:
    {root}/{subset}/{locus_class}/theta_estimates.tsv
    {root}/{subset}/{locus_class}/sfs.npy
    {root}/{subset}/{locus_class}/sfs_summary.tsv
    {root}/{subset}/{locus_class}/gain_loss_gof.tsv

Where locus_class ∈ {"introner", "non_introner_intron"}.

Produces:
    {out_tsv}: side-by-side numeric summary
    {out_png}: 4-panel comparison plot
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True,
                   help="root output directory (e.g. results/evolution/sfs_inference)")
    p.add_argument("--subsets", required=True,
                   help="Comma-separated list of subset names to compare")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--out_png", required=True)
    return p.parse_args()


def load_kv(path):
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["key"], df["value"]))


def parse_gof(path):
    """Parse the '# G=...  df=...  pval=...' line from a gain_loss_gof.tsv."""
    for line in open(path):
        if line.startswith("# G="):
            parts = line.strip().split()
            return float(parts[1].split("=")[1]), float(parts[3].split("=")[1])
    return float("nan"), float("nan")


def collect(base):
    """Return a dict of summary fields for one (subset, locus_class) directory."""
    intr = load_kv(base / "sfs_summary.tsv")
    theta = load_kv(base / "theta_estimates.tsv")
    sfs_full = np.load(base / "sfs.npy")
    G, pval = parse_gof(base / "gain_loss_gof.tsv")
    n = int(intr["n_samples"])
    seg = sfs_full[1:n]
    fixed_P = int(intr["fixed_present"])
    fixed_A = int(intr["fixed_absent"])
    return {
        "n": n,
        "n_segregating": int(seg.sum()),
        "fixed_present": fixed_P,
        "fixed_absent": fixed_A,
        "fixed_P_over_A": (fixed_P / fixed_A) if fixed_A > 0 else float("inf"),
        "theta_lambda": float(theta["theta_lambda"]),
        "theta_mu": float(theta["theta_mu"]),
        "ratio": float(theta["theta_ratio_lambda_over_mu"]),
        "GoF_G": G,
        "GoF_p": pval,
        "sfs_seg": seg,
    }


def main():
    args = parse_args()
    root = Path(args.root)
    subsets = [s.strip() for s in args.subsets.split(",") if s.strip()]
    locus_classes = ["introner", "non_introner_intron"]

    rows = []
    data = {}  # data[(subset, locus_class)] = collected dict

    for subset in subsets:
        for cls in locus_classes:
            base = root / subset / cls
            if not (base / "theta_estimates.tsv").exists():
                raise SystemExit(f"Missing fit at {base}")
            d = collect(base)
            data[(subset, cls)] = d
            rows.append({
                "subset": subset,
                "locus_class": cls,
                "n": d["n"],
                "n_segregating": d["n_segregating"],
                "fixed_present": d["fixed_present"],
                "fixed_absent": d["fixed_absent"],
                "fixed_P_over_A": (
                    f"{d['fixed_P_over_A']:.3f}"
                    if np.isfinite(d["fixed_P_over_A"]) else "inf"
                ),
                "theta_lambda": round(d["theta_lambda"], 4),
                "theta_mu": round(d["theta_mu"], 4),
                "ratio": round(d["ratio"], 4),
                "GoF_G": round(d["GoF_G"], 2),
                "GoF_p": f"{d['GoF_p']:.3g}",
            })

    df = pd.DataFrame(rows)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(df.to_string(index=False))

    # ---- 4-panel plot ----------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    cls_color = {"introner": "tab:purple", "non_introner_intron": "tab:orange"}

    # Panel A: SFS shapes (normalized) — bold for last subset, faint for earlier
    ax = axes[0, 0]
    for i, subset in enumerate(subsets):
        bold = (i == len(subsets) - 1)
        for cls in locus_classes:
            d = data[(subset, cls)]
            seg = d["sfs_seg"]
            x = np.arange(1, len(seg) + 1) / d["n"]
            y = seg / max(seg.sum(), 1)
            ax.plot(
                x, y, "o-",
                alpha=1.0 if bold else 0.4,
                linewidth=2 if bold else 1,
                color=cls_color[cls],
                label=f"{cls} ({subset})" if bold else None,
            )
    ax.set_xlabel("Present-allele frequency")
    ax.set_ylabel("Fraction of polymorphic loci")
    ax.set_title(f"SFS shape (bold = {subsets[-1]}, faint = others)")
    ax.legend(fontsize=8)

    # Panel B: theta_lambda / theta_mu
    ax = axes[0, 1]
    x = np.arange(len(subsets))
    w = 0.35
    intr_ratios = [data[(s, "introner")]["ratio"] for s in subsets]
    nonintr_ratios = [data[(s, "non_introner_intron")]["ratio"] for s in subsets]
    ax.bar(x - w / 2, intr_ratios, w, label="introners",
           color="tab:purple", edgecolor="black")
    ax.bar(x + w / 2, nonintr_ratios, w, label="non-introner introns",
           color="tab:orange", edgecolor="black")
    ax.axhline(1.0, color="black", lw=0.5, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(subsets, rotation=20, fontsize=9)
    ax.set_ylabel(r"$\theta_\lambda / \theta_\mu$")
    ax.set_title("Contemporary gain/loss ratio (segregating SFS)")
    ax.legend(fontsize=9)
    for i, (a, b) in enumerate(zip(intr_ratios, nonintr_ratios)):
        ax.text(i - w / 2, a, f"{a:.2f}", ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, b, f"{b:.2f}", ha="center", va="bottom", fontsize=8)

    # Panel C: fixed-class
    ax = axes[1, 0]
    intr_PA = [data[(s, "introner")]["fixed_P_over_A"] for s in subsets]
    intr_PA_finite = [v for v in intr_PA if np.isfinite(v)]
    height_inf = (max(intr_PA_finite) * 1.5) if intr_PA_finite else 100
    nonintr_PA = [data[(s, "non_introner_intron")]["fixed_P_over_A"] for s in subsets]
    nonintr_PA_disp = [
        height_inf if not np.isfinite(v) else v for v in nonintr_PA
    ]
    ax.bar(x - w / 2, intr_PA, w, label="introners",
           color="tab:purple", edgecolor="black")
    ax.bar(x + w / 2, nonintr_PA_disp, w,
           label="non-introner (∞ shown as bar height)",
           color="tab:orange", edgecolor="black", hatch="//")
    ax.axhline(1.0, color="black", lw=0.5, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(subsets, rotation=20, fontsize=9)
    ax.set_ylabel(r"fixed_P / fixed_A")
    ax.set_title(
        "Long-term integrated ratio "
        "(non-introner fixed_A=0 by ascertainment)"
    )
    ax.legend(fontsize=9)
    for i, s in enumerate(subsets):
        v_intr = data[(s, "introner")]["fixed_P_over_A"]
        v_nonintr = data[(s, "non_introner_intron")]["fixed_P_over_A"]
        ax.text(i - w / 2, v_intr,
                f"{v_intr:.2f}" if np.isfinite(v_intr) else "inf",
                ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, nonintr_PA_disp[i],
                "inf" if not np.isfinite(v_nonintr) else f"{v_nonintr:.2f}",
                ha="center", va="bottom", fontsize=8)

    # Panel D: text summary
    ax = axes[1, 1]
    ax.axis("off")
    lines = []
    for s in subsets:
        intr = data[(s, "introner")]
        nonintr = data[(s, "non_introner_intron")]
        lines.append(
            f"--- {s} (n={intr['n']}) ---\n"
            f"  introners:    seg={intr['n_segregating']:>4}  "
            f"P/A={intr['fixed_P_over_A']:.2f}  "
            f"ratio={intr['ratio']:.3f}  GoF_p={intr['GoF_p']:.2g}\n"
            f"  non-intr int: seg={nonintr['n_segregating']:>4}  "
            f"P/A=inf            "
            f"ratio={nonintr['ratio']:.3f}  GoF_p={nonintr['GoF_p']:.2g}\n"
        )
    ax.text(0, 1, "\n".join(lines), family="monospace", fontsize=8,
            va="top", ha="left", transform=ax.transAxes)
    ax.set_title("Per-subset summary")

    fig.suptitle(
        "Introners vs non-introner introns — SFS-based gain/loss inference",
        fontsize=12, y=1.0,
    )
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {args.out_tsv}, {args.out_png}")


if __name__ == "__main__":
    main()

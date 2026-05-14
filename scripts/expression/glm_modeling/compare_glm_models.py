#!/usr/bin/env python3
"""
Compare Poisson and Negative Binomial GLM fits and plot fit statistics.
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def extract_fit_stats(filepath):
    stats = {}
    content = Path(filepath).read_text()
    aic_match = re.search(r"AIC[:\s]+([\d.]+)", content)
    bic_match = re.search(r"BIC[:\s]+([\d.]+)", content)
    if aic_match:
        stats["AIC"] = float(aic_match.group(1))
    if bic_match:
        stats["BIC"] = float(bic_match.group(1))
    return stats


def write_summary(output_path, poisson_stats, negbin_stats):
    with open(output_path, "w") as handle:
        handle.write("GLM Model Comparison\n")
        handle.write("=" * 50 + "\n\n")
        handle.write("Poisson GLM:\n")
        handle.write(f"  AIC: {poisson_stats.get('AIC', 'N/A')}\n")
        handle.write(f"  BIC: {poisson_stats.get('BIC', 'N/A')}\n\n")
        handle.write("Negative Binomial GLM:\n")
        handle.write(f"  AIC: {negbin_stats.get('AIC', 'N/A')}\n")
        handle.write(f"  BIC: {negbin_stats.get('BIC', 'N/A')}\n\n")

        poisson_aic = poisson_stats.get("AIC")
        negbin_aic = negbin_stats.get("AIC")
        if poisson_aic is not None and negbin_aic is not None:
            if negbin_aic < poisson_aic:
                handle.write("Recommendation: Negative Binomial (lower AIC)\n")
            else:
                handle.write("Recommendation: Poisson (lower AIC)\n")
        else:
            handle.write("Recommendation: unavailable (missing AIC values)\n")


def plot_fit_stats(output_path, poisson_stats, negbin_stats):
    metrics = ["AIC", "BIC"]
    poisson_values = [poisson_stats.get(metric, np.nan) for metric in metrics]
    negbin_values = [negbin_stats.get(metric, np.nan) for metric in metrics]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(metrics))
    width = 0.35

    ax.bar(x - width / 2, poisson_values, width, label="Poisson", color="#4c78a8")
    ax.bar(x + width / 2, negbin_values, width, label="Negative Binomial", color="#f58518")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Score")
    ax.set_title("GLM Model Comparison")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)

    for idx, value in enumerate(poisson_values):
        if np.isfinite(value):
            ax.text(idx - width / 2, value, f"{value:.1f}", ha="center", va="bottom", fontsize=9)
    for idx, value in enumerate(negbin_values):
        if np.isfinite(value):
            ax.text(idx + width / 2, value, f"{value:.1f}", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Compare GLM summary statistics.")
    parser.add_argument("--poisson", required=True, help="Poisson GLM summary text file")
    parser.add_argument("--negbin", required=True, help="Negative binomial GLM summary text file")
    parser.add_argument("--output", required=True, help="Output comparison summary text file")
    parser.add_argument("--plot", required=True, help="Output comparison plot PDF")
    args = parser.parse_args()

    poisson_stats = extract_fit_stats(args.poisson)
    negbin_stats = extract_fit_stats(args.negbin)

    write_summary(args.output, poisson_stats, negbin_stats)
    plot_fit_stats(args.plot, poisson_stats, negbin_stats)


if __name__ == "__main__":
    main()

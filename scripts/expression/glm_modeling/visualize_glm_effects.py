#!/usr/bin/env python3
"""
Create GLM coefficient visualizations.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def load_coefficients(path):
    df = pd.read_csv(path)
    if "variable" not in df.columns or "coefficient" not in df.columns:
        raise ValueError(f"Expected 'variable' and 'coefficient' columns in {path}")
    return df[df["variable"] != "Intercept"].copy()


def save_dual_format(fig, pdf_path, png_path):
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_effect_plot(poisson, negbin, pdf_path, png_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axes[0].barh(poisson["variable"], poisson["coefficient"], color="#4c78a8")
    axes[0].axvline(x=0, color="k", linestyle="--", alpha=0.5)
    axes[0].set_xlabel("Coefficient")
    axes[0].set_title("Poisson GLM Coefficients")

    axes[1].barh(negbin["variable"], negbin["coefficient"], color="#f58518")
    axes[1].axvline(x=0, color="k", linestyle="--", alpha=0.5)
    axes[1].set_xlabel("Coefficient")
    axes[1].set_title("Negative Binomial GLM Coefficients")

    fig.tight_layout()
    save_dual_format(fig, pdf_path, png_path)


def make_forest_plot(poisson, negbin, pdf_path, png_path):
    combined = pd.concat(
        [
            poisson.assign(model="Poisson"),
            negbin.assign(model="Negative Binomial"),
        ],
        ignore_index=True,
    )
    combined["label"] = combined["model"] + ": " + combined["variable"].astype(str)

    fig_height = max(6, 0.35 * len(combined))
    fig, ax = plt.subplots(figsize=(10, fig_height))

    colors = combined["model"].map({"Poisson": "#4c78a8", "Negative Binomial": "#f58518"})
    ax.scatter(combined["coefficient"], combined["label"], c=colors)
    ax.axvline(x=0, color="k", linestyle="--", alpha=0.5)
    ax.set_xlabel("Coefficient")
    ax.set_ylabel("")
    ax.set_title("GLM Coefficient Forest Plot")
    ax.grid(axis="x", alpha=0.25)

    fig.tight_layout()
    save_dual_format(fig, pdf_path, png_path)


def make_error_plot(message, pdf_path, png_path):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.axis("off")
    fig.tight_layout()
    save_dual_format(fig, pdf_path, png_path)


def main():
    parser = argparse.ArgumentParser(description="Visualize GLM coefficients.")
    parser.add_argument("--poisson", required=True, help="Poisson coefficient CSV")
    parser.add_argument("--negbin", required=True, help="Negative binomial coefficient CSV")
    parser.add_argument("--effect-pdf", required=True, help="Effect size plot PDF")
    parser.add_argument("--effect-png", required=True, help="Effect size plot PNG")
    parser.add_argument("--forest-pdf", required=True, help="Forest plot PDF")
    parser.add_argument("--forest-png", required=True, help="Forest plot PNG")
    args = parser.parse_args()

    try:
        poisson = load_coefficients(args.poisson)
        negbin = load_coefficients(args.negbin)
        make_effect_plot(poisson, negbin, args.effect_pdf, args.effect_png)
        make_forest_plot(poisson, negbin, args.forest_pdf, args.forest_png)
    except Exception as exc:
        message = f"Could not generate plot: {exc}"
        make_error_plot(message, args.effect_pdf, args.effect_png)
        make_error_plot(message, args.forest_pdf, args.forest_png)


if __name__ == "__main__":
    main()

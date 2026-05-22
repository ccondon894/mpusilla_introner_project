#!/usr/bin/env python3
"""Plot GLM coefficient forests and model-predicted introner effects."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Poisson


INTRONER_TERMS = [
    ("baseline_introner_count", "Baseline\nintroner count"),
    ("introner_gain", "Introner gain"),
    ("introner_loss", "Introner loss"),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poisson-data", required=True)
    parser.add_argument("--nb-data", required=True)
    parser.add_argument("--poisson-coefficients", required=True)
    parser.add_argument("--nb-coefficients", required=True)
    parser.add_argument("--forest-pdf", required=True)
    parser.add_argument("--forest-png", required=True)
    parser.add_argument("--prediction-pdf", required=True)
    parser.add_argument("--prediction-png", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def load_coefs(path, model_name, response):
    df = pd.read_csv(path)
    df = df[df["variable"].isin([term for term, _label in INTRONER_TERMS])].copy()
    df["model"] = model_name
    df["response"] = response
    df["term_label"] = df["variable"].map(dict(INTRONER_TERMS))
    for col in ["coefficient", "ci_lower", "ci_upper"]:
        df[f"{col}_percent"] = (np.exp(df[col]) - 1) * 100
    return df


def plot_forest(poisson_coef, nb_coef, output_pdf, output_png):
    forest = pd.concat([
        load_coefs(poisson_coef, "Poisson GLM", "Isoform count"),
        load_coefs(nb_coef, "Negative binomial GLM", "Expression"),
    ], ignore_index=True)

    sns.set_style("white")
    colors = {
        "Isoform count": "#0072B2",
        "Expression": "#D55E00",
    }
    order = [label for _term, label in INTRONER_TERMS]
    offsets = {
        "Isoform count": 0.12,
        "Expression": -0.12,
    }

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    y_lookup = {label: i for i, label in enumerate(order)}

    for _, row in forest.iterrows():
        y = y_lookup[row["term_label"]] + offsets[row["response"]]
        x = row["coefficient_percent"]
        lo = row["ci_lower_percent"]
        hi = row["ci_upper_percent"]
        ax.errorbar(
            x,
            y,
            xerr=[[x - lo], [hi - x]],
            fmt="o",
            color=colors[row["response"]],
            ecolor=colors[row["response"]],
            elinewidth=1.6,
            capsize=3,
            markersize=6,
            label=row["response"],
        )

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(
        unique.values(),
        unique.keys(),
        frameon=True,
        edgecolor="0.35",
        facecolor="white",
        framealpha=1,
        loc="upper right",
    )
    ax.axvline(0, color="0.25", linewidth=1.0, linestyle="--")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlabel("Expected response change per unit increase (%)")
    ax.set_ylabel("")
    ax.set_title("Model-estimated introner effects")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return forest


def fit_poisson(data_path):
    data = pd.read_csv(data_path)
    data = data[data["n_isoforms"] <= 20].copy()
    formula = (
        "n_isoforms ~ C(strain) + introner_gain + introner_loss + "
        "baseline_introner_count + log_expression + log_cds_length"
    )
    model = smf.glm(formula, data=data, family=Poisson()).fit()
    return data, model


def fit_nb(data_path):
    data = pd.read_csv(data_path)
    formula = (
        "raw_count ~ C(strain, Treatment(reference='CCMP1545')) + "
        "introner_gain + introner_loss + baseline_introner_count + "
        "log_cds_length + n_isoforms + GC_content"
    )
    model = smf.glm(
        formula,
        data=data,
        family=sm.families.NegativeBinomial(),
        offset=data["log_library_size"],
    ).fit()
    return data, model


def typical_row(data, columns):
    row = {}
    for col in columns:
        if col == "strain":
            row[col] = "CCMP1545"
        elif col == "replicate":
            row[col] = "834_A1"
        elif pd.api.types.is_numeric_dtype(data[col]):
            row[col] = data[col].median()
        else:
            row[col] = data[col].mode().iloc[0]
    return row


def make_prediction_grid(data, variable, max_value, columns):
    base = typical_row(data, columns)
    rows = []
    for value in range(max_value + 1):
        row = base.copy()
        for term, _label in INTRONER_TERMS:
            row[term] = 0
        row[variable] = value
        rows.append(row)
    return pd.DataFrame(rows)


def plot_predictions(poisson_data, poisson_model, nb_data, nb_model,
                     output_pdf, output_png):
    sns.set_style("white")
    colors = {
        "baseline_introner_count": "#009E73",
        "introner_gain": "#0072B2",
        "introner_loss": "#D55E00",
    }
    maxima = {
        "baseline_introner_count": 5,
        "introner_gain": 4,
        "introner_loss": 4,
    }

    fig, axes = plt.subplots(
        2, 3, figsize=(13.2, 7.2), sharex="col", sharey="row"
    )

    poisson_cols = [
        "strain", "introner_gain", "introner_loss",
        "baseline_introner_count", "log_expression", "log_cds_length",
    ]
    for col_idx, (variable, label) in enumerate(INTRONER_TERMS):
        grid = make_prediction_grid(
            poisson_data, variable, maxima[variable], poisson_cols
        )
        pred = poisson_model.get_prediction(grid).summary_frame()
        x = grid[variable].to_numpy(dtype=float)
        ax = axes[0, col_idx]
        ax.plot(x, pred["mean"], marker="o", color=colors[variable],
                label=label)
        ax.fill_between(
            x,
            pred["mean_ci_lower"].to_numpy(dtype=float),
            pred["mean_ci_upper"].to_numpy(dtype=float),
            color=colors[variable],
            alpha=0.18,
            linewidth=0,
        )
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax.set_xlabel("")

    nb_cols = [
        "strain", "introner_gain", "introner_loss",
        "baseline_introner_count", "log_cds_length", "n_isoforms",
        "GC_content", "library_size", "log_library_size",
    ]
    median_library = nb_data["library_size"].median()
    for col_idx, (variable, label) in enumerate(INTRONER_TERMS):
        grid = make_prediction_grid(nb_data, variable, maxima[variable], nb_cols)
        grid["library_size"] = median_library
        grid["log_library_size"] = np.log(median_library)
        pred = nb_model.get_prediction(
            grid, offset=grid["log_library_size"]
        ).summary_frame()
        x = grid[variable].to_numpy(dtype=float)
        mean_cpm = pred["mean"].to_numpy(dtype=float) / median_library * 1_000_000
        lower_cpm = pred["mean_ci_lower"].to_numpy(dtype=float) / median_library * 1_000_000
        upper_cpm = pred["mean_ci_upper"].to_numpy(dtype=float) / median_library * 1_000_000
        ax = axes[1, col_idx]
        ax.plot(x, mean_cpm, marker="o", color=colors[variable],
                label=label)
        ax.fill_between(
            x,
            lower_cpm,
            upper_cpm,
            color=colors[variable],
            alpha=0.18,
            linewidth=0,
        )
        ax.set_xlabel(label)

    axes[0, 0].set_ylabel("Poisson model\nPredicted isoforms per gene")
    axes[1, 0].set_ylabel("Negative binomial model\nPredicted expression (CPM)")

    for ax in axes.flatten():
        ax.grid(False)
        ax.set_xticks(np.arange(0, int(ax.get_xlim()[1]) + 1, 1))

    fig.suptitle(
        "Model-predicted introner effects with other covariates held at typical values"
    )
    fig.subplots_adjust(hspace=0, wspace=0, top=0.88, bottom=0.12)
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(path, forest):
    with open(path, "w") as handle:
        handle.write("GLM Introner Effect Visualization\n")
        handle.write("=" * 72 + "\n\n")
        handle.write("Forest plot terms are exponentiated GLM coefficients.\n")
        handle.write("Values are percent change in expected response per one-unit increase.\n\n")
        cols = [
            "response", "term_label", "coefficient_percent",
            "ci_lower_percent", "ci_upper_percent", "p_value",
        ]
        handle.write(forest[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
        handle.write("\n\n")
        handle.write("Prediction plot holds non-focal numeric covariates at their medians ")
        handle.write("and strain at CCMP1545. Negative-binomial predictions use the median ")
        handle.write("library size and are shown as CPM.\n")


def main():
    args = parse_args()
    forest = plot_forest(
        args.poisson_coefficients,
        args.nb_coefficients,
        args.forest_pdf,
        args.forest_png,
    )
    poisson_data, poisson_model = fit_poisson(args.poisson_data)
    nb_data, nb_model = fit_nb(args.nb_data)
    plot_predictions(
        poisson_data,
        poisson_model,
        nb_data,
        nb_model,
        args.prediction_pdf,
        args.prediction_png,
    )
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    write_report(args.report, forest)


if __name__ == "__main__":
    main()

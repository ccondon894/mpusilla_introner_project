#!/usr/bin/env python3
"""Plot the primary isoform, gain, NMD and retention model estimates."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figure_color_guide import load_color_guide


def select_one(table, **filters):
    for column, value in filters.items():
        table = table[table[column].eq(value)]
    if len(table) != 1:
        raise ValueError(f'Expected one coefficient for {filters}, found {len(table)}')
    return table.iloc[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['isoform', 'gain', 'nmd', 'retention', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    colors = load_color_guide()
    iso = select_one(pd.read_csv(args.isoform, sep='\t'), term='current_introner_count')
    gain = select_one(pd.read_csv(args.gain, sep='\t'), term='reference_relative_gain_count')
    nmd = select_one(pd.read_csv(args.nmd, sep='\t'), model='binary_architecture_adjusted', term='has_introner')
    retention = pd.read_csv(args.retention, sep='\t')
    retention = retention[retention.endpoint.eq('combined_boundary_expression_interaction')]
    if retention.empty:
        raise ValueError('Primary adjusted retention endpoint is missing')
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    rows = [(iso, 'rate_ratio', 'rate_ratio_ci_lower', 'rate_ratio_ci_upper', 'Isoform richness', 'Isoform count ratio'),
            (gain, 'rate_ratio', 'rr_ci_lower', 'rr_ci_upper', 'Within-gene gain', 'Expression rate ratio'),
            (nmd, 'odds_ratio', 'or_ci_lower', 'or_ci_upper', 'NMD association', 'Odds ratio')]
    for ax, (r, estimate, low, high, title, label) in zip(axes.flat, rows):
        ax.errorbar(r[estimate], 0, xerr=[[r[estimate]-r[low]], [r[high]-r[estimate]]], fmt='o', capsize=4, color=colors['Population 1'])
        ax.axvline(1, color=colors['Ancestor'], linestyle='--'); ax.set_xscale('log'); ax.set_yticks([])
        ax.set_xlabel(label+' (95% CI)'); ax.set_title(title)
        ax.text(.02, .9, f'P = {r.p_value:.3g}', transform=ax.transAxes)
    ax = axes[1, 1]
    for i, r in enumerate(retention.itertuples()):
        ax.errorbar(100*r.adjusted_retention_rate, i,
            xerr=[[100*(r.adjusted_retention_rate-r.ci_low)], [100*(r.ci_high-r.adjusted_retention_rate)]], fmt='o', capsize=4, color=colors['Population 1'])
    ax.set_yticks(range(len(retention)), retention.class_label)
    ax.set_xlabel('Adjusted retention (%)'); ax.set_title('Expression-adjusted retention')
    fig.tight_layout(); args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300); plt.close(fig)


if __name__ == '__main__': main()

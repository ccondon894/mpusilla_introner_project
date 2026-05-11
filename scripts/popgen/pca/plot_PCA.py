#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

def plot_pca(eigenval_file, eigenvec_file, output_pdf, output_png, group2_samples):
    # Load the data
    eigenvals = pd.read_csv(eigenval_file, header=None, names=['eigenval'])
    eigenvecs = pd.read_csv(eigenvec_file, sep=r'\s+', header=None)

    # Extract sample names and PC coordinates
    samples = eigenvecs.iloc[:, 1]  # Second column contains sample names
    pc1 = eigenvecs.iloc[:, 2]      # Third column is PC1
    pc2 = eigenvecs.iloc[:, 3]      # Fourth column is PC2
    pc3 = eigenvecs.iloc[:, 4]      # Fifth column is PC3

    # Calculate variance explained
    total_variance = eigenvals['eigenval'].sum()
    pc1_var = (eigenvals.iloc[0, 0] / total_variance) * 100
    pc2_var = (eigenvals.iloc[1, 0] / total_variance) * 100
    pc3_var = (eigenvals.iloc[2, 0] / total_variance) * 100

    # Create colors - red for group2, blue for others
    colors = ['red' if sample in group2_samples else 'blue' for sample in samples]

    # Set up the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Plot PC1 vs PC2
    ax1.scatter(pc1, pc2, c=colors, alpha=0.7, s=80)
    ax1.set_xlabel(f'PC1 ({pc1_var:.1f}%)')
    ax1.set_ylabel(f'PC2 ({pc2_var:.1f}%)')
    ax1.set_title('PC1 vs PC2')

    # Add sample labels with adjustText for non-overlapping placement
    texts1 = []
    for i, sample in enumerate(samples):
        texts1.append(ax1.text(pc1.iloc[i], pc2.iloc[i], sample, fontsize=8))
    adjust_text(texts1, ax=ax1, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5))

    # Plot PC1 vs PC3
    ax2.scatter(pc1, pc3, c=colors, alpha=0.7, s=80)
    ax2.set_xlabel(f'PC1 ({pc1_var:.1f}%)')
    ax2.set_ylabel(f'PC3 ({pc3_var:.1f}%)')
    ax2.set_title('PC1 vs PC3')

    # Add sample labels with adjustText for non-overlapping placement
    texts2 = []
    for i, sample in enumerate(samples):
        texts2.append(ax2.text(pc1.iloc[i], pc3.iloc[i], sample, fontsize=8))
    adjust_text(texts2, ax=ax2, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5))

    ax1.grid(True, alpha=0.3)
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight', format="pdf")
    plt.savefig(output_png, dpi=300, bbox_inches='tight', format="png")


    print(f"Total samples: {len(samples)}")
    print(f"Group2 samples found: {sum(1 for sample in samples if sample in group2_samples)}")
    print(f"PC1 variance explained: {pc1_var:.2f}%")
    print(f"PC2 variance explained: {pc2_var:.2f}%")
    print(f"PC3 variance explained: {pc3_var:.2f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot PCA from PLINK eigenval/eigenvec files.')
    parser.add_argument('--eigenval', required=True, help='Path to .eigenval file')
    parser.add_argument('--eigenvec', required=True, help='Path to .eigenvec file')
    parser.add_argument('--output_pdf', required=True, help='Output plot path (extension determines format)')
    parser.add_argument('--output_png', required=True, help='Output plot path (extension determines format)')
    parser.add_argument('--group1', required=True, help='Comma-separated Group1 sample names')
    parser.add_argument('--group2', required=True, help='Comma-separated Group2 sample names')
    args = parser.parse_args()

    group2_samples = set(args.group2.split(','))
    plot_pca(args.eigenval, args.eigenvec, args.output_pdf, args.output_png, group2_samples)
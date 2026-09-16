import argparse
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

def parse_args():
    parser = argparse.ArgumentParser(description="Compute and plot 2D folded allele frequency spectrum from VCF.")
    parser.add_argument('--vcf', type=str, required=True, help="Path to the VCF file")
    parser.add_argument('--group1', type=str, required=True, help="Comma-separated Group1 sample names")
    parser.add_argument('--group2', type=str, required=True, help="Comma-separated Group2 sample names (used as pop1)")
    parser.add_argument('--output_png', type=str, required=True, help="Output plot path (extension determines format)")
    parser.add_argument('--output_pdf', type=str, required=True, help="Output plot path (extension determines format)")
    parser.add_argument('--figure-width', type=float, default=6.5, help="Figure width in inches")
    parser.add_argument('--figure-height', type=float, default=2.7, help="Figure height in inches")
    parser.add_argument('--font-size', type=float, default=10.0, help="Base font size in points")
    return parser.parse_args()

def compute_folded_allele_frequency_spectrum(vcf_file, pop1_samples, obsolete_samples):
    
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Get the index of the samples in the VCF
    sample_names = list(vcf.header.samples)
    pop1_indices = [sample_names.index(sample) for sample in pop1_samples if sample in sample_names]
    pop2_indices = [i for i in range(len(sample_names)) if sample_names[i] not in pop1_samples and sample_names[i] not in obsolete_samples]

    # Initialize a 2D matrix to hold the counts
    # With global minor allele folding, max minor allele count = floor(total_samples / 2)
    total_samples = len(pop1_indices) + len(pop2_indices)
    max_minor_allele_count = total_samples // 2

    # Pop1 can have 0 to len(pop1_indices) minor alleles, but capped by max_minor_allele_count
    max_pop1_minor = min(len(pop1_indices), max_minor_allele_count)
    # Pop2 can have 0 to max_minor_allele_count minor alleles
    max_pop2_minor = max_minor_allele_count

    spectrum = np.zeros((max_pop1_minor + 1, max_pop2_minor + 1), dtype=int)

    # Iterate through the VCF records
    for record in vcf:
        if record.alts and len(record.alts) == 1:  # Biallelic SNP
            # Count alternative alleles in both populations
            pop1_alt_count = sum(1 for i in pop1_indices if record.samples[i]['GT'][0] and record.samples[i]['GT'][0] > 0)
            pop2_alt_count = sum(1 for i in pop2_indices if record.samples[i]['GT'][0] and record.samples[i]['GT'][0] > 0)

            # Compute the total alternative allele count across both populations
            total_alt_count = pop1_alt_count + pop2_alt_count
            total_samples = len(pop1_indices) + len(pop2_indices)

            # Fold the counts by taking the minor allele count across the total population
            minor_allele_count = min(total_alt_count, total_samples - total_alt_count)

            # Skip monomorphic sites
            if minor_allele_count == 0:
                continue

            # Determine how many of the minor alleles are in each population
            if total_alt_count == minor_allele_count:
                pop1_minor_count = pop1_alt_count
                pop2_minor_count = pop2_alt_count
            else:
                pop1_minor_count = len(pop1_indices) - pop1_alt_count
                pop2_minor_count = len(pop2_indices) - pop2_alt_count

            # Update the 2D folded spectrum
            spectrum[pop1_minor_count, pop2_minor_count] += 1

    return spectrum

def print_spectrum(spectrum):
    print("Folded 2D Allele Frequency Spectrum (Global Minor Allele Folding):")
    print("Population 1 (Group2 - RCC1749, RCC3052) vs Population 2 (Group1 - Intronerful, excluding obsolete)")
    print("Minor allele determined jointly across both populations")
    print("Minor Allele Counts in Population 1 | Counts in Population 2")
    print("-----------------------------------------------------------")
    
    for pop1_count, row in enumerate(spectrum):
        row_str = f"Minor Allele Count {pop1_count}: {row}"
        print(row_str)

def plot_2d_afs_heatmap(
    spectrum,
    output_pdf,
    output_png,
    figure_width,
    figure_height,
    font_size,
):
    # Log scale the counts, setting zeros to a small positive value to avoid log(0)
    log_spectrum = np.log10(spectrum + 1)  # Adding 1 to avoid log(0)

    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    ax = sns.heatmap(log_spectrum, annot=log_spectrum, fmt=".1f", cmap='viridis',
                     annot_kws={'fontsize': font_size - 1},
                     cbar_kws={'label': 'log10(count)', 'pad': 0.025, 'aspect': 16},
                     linewidths=0.5, linecolor='white',
                     square=False, ax=ax)

    # Add clean labels (no title for publication-quality)
    ax.set_xlabel("Population 1 minor allele count", fontsize=font_size + 1)
    ax.set_ylabel("Population 2\nminor allele count", fontsize=font_size + 1)
    # plt.title("Folded 2D Allele Frequency Spectrum")

    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center', fontsize=font_size)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=font_size)
    # Set y-axis to show values correctly (inverted to match typical AFS orientation)
    ax.invert_yaxis()

    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=font_size - 1)
    colorbar.set_label('log10(count)', fontsize=font_size)

    fig.tight_layout(pad=0.45)

    # Save the plot with high DPI
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.03)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.03)


def main():
    args = parse_args()
    pop1_samples = args.group2.split(',')
    obsolete_samples = {'CCMP490', 'RCC647', 'RCC835'}
    spectrum = compute_folded_allele_frequency_spectrum(args.vcf, pop1_samples, obsolete_samples)

    print_spectrum(spectrum)
    plot_2d_afs_heatmap(
        spectrum,
        args.output_pdf,
        args.output_png,
        args.figure_width,
        args.figure_height,
        args.font_size,
    )

if __name__ == "__main__":
    main()

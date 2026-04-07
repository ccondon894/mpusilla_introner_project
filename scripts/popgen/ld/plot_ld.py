import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="Plot LD (mean r2) vs. SNP distance.")
    parser.add_argument('input_file', type=str, help="Path to the LD summary file (output from previous script)")
    parser.add_argument('output_file', type=str, help="Path to the output plot image file (e.g., ld_vs_distance.png)")
    return parser.parse_args()

def plot_ld_vs_distance(input_file, output_file):
    # Read the input file with LD summary
    ld_summary = pd.read_csv(input_file, sep='\t')
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(ld_summary['mean_distance'], ld_summary['mean_r2'], linestyle='-', color='b')
    
    # Add labels and title
    plt.xlabel('Mean Distance Between SNPs (bp)', fontsize=14)
    plt.ylabel('Mean $r^2$', fontsize=14)
    plt.title('Linkage Disequilibrium (LD) vs. Distance Between SNPs', fontsize=16)
    
    # Save the plot to the output file
    plt.savefig(output_file)
    plt.close()

def main():
    args = parse_args()
    plot_ld_vs_distance(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

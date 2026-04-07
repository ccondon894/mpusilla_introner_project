import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize LD in bins based on SNP pair distances.")
    parser.add_argument('input_file', type=str, help="Path to the LD results file")
    parser.add_argument('output_file', type=str, help="Path to the output file")
    parser.add_argument('--bin_size', type=int, default=10000, help="Size of distance bins (default=10000 bp)")
    return parser.parse_args()

def summarize_ld_by_distance_bins(input_file, output_file, bin_size):
    # Read the input LD results
    ld_data = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'position1', 'position2', 'D', 'r2'])
    
    # Calculate the distance between SNPs (|position1 - position2|)
    ld_data['distance'] = np.abs(ld_data['position1'] - ld_data['position2'])
    
    # Create bins based on distance
    ld_data['distance_bin'] = (ld_data['distance'] // bin_size) * bin_size
    
    # Group by the distance bin and compute mean statistics
    summary = ld_data.groupby('distance_bin').agg(
        mean_distance=('distance', 'mean'),
        mean_r2=('r2', 'mean'),
        num_entries=('distance', 'size')
    ).reset_index()

    # Write the result to the output file
    summary.to_csv(output_file, sep='\t', index=False, header=["distance_bin", "mean_distance", "mean_r2", "num_entries"])

def main():
    args = parse_args()
    summarize_ld_by_distance_bins(args.input_file, args.output_file, args.bin_size)

if __name__ == "__main__":
    main()

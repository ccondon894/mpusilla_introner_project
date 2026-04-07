import numpy as np
import matplotlib.pyplot as plt
import argparse

# Function to compute 95% confidence intervals
def compute_confidence_intervals(data):
    lower_bounds = np.percentile(data, 2.5, axis=0)
    upper_bounds = np.percentile(data, 97.5, axis=0)
    return lower_bounds, upper_bounds

# Function to create scatter plots with confidence intervals and truth
def plot_scatter_matrix(bootstrap_data, param_names, lower_bounds, upper_bounds, param_truth, output_file):
    num_params = len(param_names)
    
    fig, axes = plt.subplots(nrows=num_params, ncols=num_params, figsize=(12, 12))
    
    for i in range(num_params):
        for j in range(num_params):
            if i != j:
                axes[i, j].scatter(bootstrap_data[:, j], bootstrap_data[:, i], alpha=0.5)
            if i == j:
                axes[i, j].text(0.5, 0.5, param_names[i], fontsize=12, ha='center')
            if j == 0:
                axes[i, j].set_ylabel(param_names[i])
            if i == num_params - 1:
                axes[i, j].set_xlabel(param_names[j])
    
    # Annotating confidence intervals with truth values
    textstr = '\n'.join([f'{param_names[i]}: {param_truth[i]:.3f} '
                         f'[{lower_bounds[i]:.3f}, {upper_bounds[i]:.3f}]'
                         for i in range(num_params)])
    fig.suptitle(f'Bootstrap 95% Confidence Intervals:\n{textstr}', fontsize=14)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_file)
    plt.show()

# Argparse function for handling command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Generate scatter plots from demographic inference bootstraps.")
    parser.add_argument('input_file', help='Path to the input file containing parameter estimates and bootstraps.')
    parser.add_argument('--output_file', default='scatter_matrix.pdf', help='Output file for the scatter matrix plot (default: scatter_matrix.pdf).')
    return parser.parse_args()

# Main function to handle the overall process
def main():
    args = parse_args()

    # Load data from the input file
    data = np.loadtxt(args.input_file)

    # Find the best fit (highest/least negative log-likelihood)
    best_idx = np.argmax(data[:, 0])
    param_truth = data[best_idx, 1:6]

    # Extract bootstrap estimates (all runs)
    bootstrap_data = data[:, 1:6]

    # Parameter names
    param_names = ['N1', 'N2', 'T', 'M', 'Theta']
    
    # Compute 95% confidence intervals
    lower_bounds, upper_bounds = compute_confidence_intervals(bootstrap_data)
    
    # Plot scatter matrix with confidence intervals and truth values
    plot_scatter_matrix(bootstrap_data, param_names, lower_bounds, upper_bounds, param_truth, args.output_file)

if __name__ == "__main__":
    main()

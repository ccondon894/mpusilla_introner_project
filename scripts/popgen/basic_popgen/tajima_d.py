import argparse
import gzip
import math

def parse_vcf(vcf_file):
    """
    Parse a VCF file to extract allele counts and sample sizes.

    Parameters:
    vcf_file (str): Path to the VCF file.

    Returns:
    tuple: (allele_counts dict, sample_size int)
    """
    allele_counts = {}
    sample_size = 0
    obsolete_samples = {'CCMP490', 'RCC647', 'RCC835'}
    intronerless_samples = {'RCC1749', 'RCC3052'}
    excluded_samples = obsolete_samples | intronerless_samples

    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#CHROM'):
                # Get sample names and filter out excluded ones (obsolete + intronerless)
                fields = line.strip().split('\t')
                sample_names = fields[9:]
                valid_samples = [i for i, name in enumerate(sample_names) if name not in excluded_samples]
                sample_size = len(valid_samples)
                continue
            elif line.startswith('#'):
                continue  # Skip other header lines
            
            fields = line.strip().split('\t')
            genotypes = fields[9:]

            # Count alternative alleles in valid samples only (haploid organism)
            allele_count = 0
            for i in valid_samples:
                gt = genotypes[i].split(':', 1)[0]  # Get genotype part
                # For haploid organisms, genotype should be a single value
                if gt not in ('0', '.'):
                    allele_count += 1
            
            if allele_count not in allele_counts:
                allele_counts[allele_count] = 1
            else:
                allele_counts[allele_count] += 1
    
    return allele_counts, sample_size

def popgen(allele_counts, sample_size):
    """
    Calculate Tajima's D from allele counts and sample size for haploid organism.
    
    Parameters:
    allele_counts (dict): Allele counts where keys are counts and values are the number of sites with that count.
    sample_size (int): Number of individuals (haploid, so equals number of chromosomes).
    
    Returns:
    tuple: (tajima_d, pi_per_site, theta_w_per_site)
    """
    if 0 in allele_counts:
        del allele_counts[0]  # Remove zero counts if they exist
    
    S = sum(allele_counts.values())  # Number of segregating sites
    n = sample_size  # For haploid organism, n individuals = n chromosomes
    
    if S == 0:
        return 0.0, 0.0, 0.0
    
    ## get a1, a2 
    a1 = 0 
    a2 = 0 
    for i in range(1, n):
        a1 += 1/i 
        a2 += 1/i**2

    ## b1, b2 
    b1 = ( n + 1 ) / ( 3 * ( n  - 1 ) )
    b2 = 2*(n**2+n+3)/(9*n*(n-1)) 

    ## c1, c2
    c1 = b1 - 1/a1 
    c2 = b2 - (n+2)/(a1*n) + a2/a1**2

    ## e1, e2
    e1 = c1/a1 
    e2 = c2/(a1**2+a2)

    # Calculate Pi (average number of pairwise differences) for haploid
    pi = 0 
    for freq, count in allele_counts.items():
        pi += ( freq * ( n - freq ) * count ) / ( n * ( n - 1 ) / 2 ) 

    # Calculate Watterson's theta (θ_W)
    theta_W = S / a1  # Watterson's estimator

    # Calculate Tajima's D
    if e1 * S + e2 * S * (S-1) > 0:
        tajima_d = ( pi - theta_W ) / math.sqrt( e1 * S + e2 * S * (S-1 ) )
    else:
        tajima_d = 0.0
    
    return tajima_d, pi/S, theta_W/S

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Calculate Tajima\'s D from a VCF file.')
    parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file (can be gzipped)')
    parser.add_argument('--output', type=str, required=True, help='Path to output TSV file')
    parser.add_argument('--summary', type=str, required=True, help='Path to summary text file')

    # Parse arguments
    args = parser.parse_args()

    # Compute Tajima's D
    allele_counts, sample_size = parse_vcf(args.vcf)
    tajd, pi, theta = popgen(allele_counts, sample_size)

    with open(args.output, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"tajimas_d\t{tajd:.6f}\n")
        f.write(f"pi_per_site\t{pi:.6f}\n")
        f.write(f"theta_w_per_site\t{theta:.6f}\n")

    with open(args.summary, 'w') as f:
        f.write(f"Tajima's D: {tajd:.6f}\n")
        f.write(f"Pi per site: {pi:.6f}\n")
        f.write(f"Theta_W per site: {theta:.6f}\n")
        f.write(f"Sample size: {sample_size} individuals (haploid)\n")

if __name__ == '__main__':
    main()

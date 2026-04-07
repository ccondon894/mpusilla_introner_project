import argparse
import gzip
import numpy as np
import pysam
from itertools import combinations

def parse_args():
    parser = argparse.ArgumentParser(description="Compute linkage disequilibrium (D and r^2) between SNPs in a VCF file.")
    parser.add_argument('--vcf', type=str, required=True, help="Path to the VCF file")
    parser.add_argument('--output', type=str, required=True, help="Path to the output file (use .gz extension for gzip output)")
    parser.add_argument('--max_distance', type=int, default=100000, help="Maximum distance between SNPs to compute LD (default: 100000)")
    return parser.parse_args()

def compute_ld(vcf_file, max_distance=100000):
    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Initialize data storage
    snps = []
    snp_positions = {}
    for record in vcf:
        if record.alts and len(record.alts) == 1:  # Biallelic SNP
            chrom = record.chrom
            pos = record.pos
            genotypes = [record.samples[sample]['GT'][0] for sample in vcf.header.samples]

            ## skip any site with missing data
            if any ( g is None for g in genotypes ) :
                continue

            # Check if site is invariant
            if len(set(genotypes)) > 1:  # Check if there are both 0s and 1s
                snps.append((chrom, pos, genotypes))
                if chrom not in snp_positions:
                    snp_positions[chrom] = []
                snp_positions[chrom].append(pos)

    # Compute LD for SNPs within max_distance bp on the same chromosome
    ld_results = []
    for chrom in snp_positions:
        positions = snp_positions[chrom]
        snps_in_region = [snp for snp in snps if snp[0] == chrom and snp[1] in positions]
        for (snp1, snp2) in combinations(snps_in_region, 2):
            pos1, genotypes1 = snp1[1], snp1[2]
            pos2, genotypes2 = snp2[1], snp2[2]
            if abs(pos1 - pos2) <= max_distance:
                D, r2 = compute_d_and_r2(genotypes1, genotypes2)
                ld_results.append((chrom, pos1, pos2, D, r2))

    return ld_results

def compute_d_and_r2(genotypes1, genotypes2):
    # Convert genotypes to numpy arrays
    g1 = np.array(genotypes1)
    g2 = np.array(genotypes2)
    
    # Compute allele counts
    p11 = np.sum((g1 == 1) & (g2 == 1))
    p10 = np.sum((g1 == 1) & (g2 == 0))
    p01 = np.sum((g1 == 0) & (g2 == 1))
    p00 = np.sum((g1 == 0) & (g2 == 0))
    
    n = len(g1)
    pA = (p11 + p10) / n
    pB = (p11 + p01) / n
    pAB = p11 / n
    pA_notB = p10 / n
    p_notA_B = p01 / n
    p_notA_notB = p00 / n
    
    # Compute D
    D = pAB - (pA * pB)
    
    # Compute r^2
    if (pA * (1 - pA)) * (pB * (1 - pB)) != 0:
        r2 = (D ** 2) / ((pA * (1 - pA)) * (pB * (1 - pB)))
    else:
        r2 = np.nan
    
    return D, r2

def main():
    args = parse_args()
    ld_results = compute_ld(args.vcf, args.max_distance)

    open_func = gzip.open if args.output.endswith('.gz') else open
    with open_func(args.output, 'wt') as f:
        for result in ld_results:
            chrom, pos1, pos2, D, r2 = result
            f.write(f"{chrom}\t{pos1}\t{pos2}\t{D:.4f}\t{r2:.4f}\n")

if __name__ == "__main__":
    main()
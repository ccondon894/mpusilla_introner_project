import argparse
import gzip

# Function to parse degenotate output
def parse_degenotate(degenotate_file):
    synonymous_sites = {}
    with gzip.open(degenotate_file, 'rt') if degenotate_file.endswith('.gz') else open(degenotate_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            scaffold = columns[0]
            position = int(columns[1])
            if columns[4] == '.' :
                continue 
            if int(columns[4]) == 4 :
                synonymous_sites[(scaffold, position+1)] = None
    return synonymous_sites

def extract_scaffold(input_string):
    return input_string.split('#')[-1]

# Function to filter VCF file
def filter_vcf(vcf_file, synonymous_sites):
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                print(line.strip())  # Output the header
                continue

            columns = line.strip().split('\t')
            scaffold = columns[0]
            position = int(columns[1])
            ref_allele = columns[3]
            alt_alleles = columns[4].split(',')

            # Check if the site is synonymous
            if (scaffold, position) in synonymous_sites:
                syn_mutations = synonymous_sites[(scaffold, position)]
                if syn_mutations is None:
                    # All mutations at this site are synonymous
                    print(line.strip())
                else:
                    # Only specific mutations are synonymous
                    for alt in alt_alleles:
                        if alt in syn_mutations:
                            print(line.strip())
                            break

# Argument parsing
def main():
    parser = argparse.ArgumentParser(description="Filter VCF file to include only synonymous SNPs.")
    parser.add_argument("degenotate", help="The degenotate BED file.")
    parser.add_argument("vcf", help="The input VCF file (can be gzipped).")
    args = parser.parse_args()

    # Parse degenotate output
    synonymous_sites = parse_degenotate(args.degenotate)

    # Filter VCF based on synonymous sites
    filter_vcf(args.vcf, synonymous_sites)

if __name__ == "__main__":
    main()

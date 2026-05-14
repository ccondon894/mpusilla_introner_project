import pandas as pd
import requests
import json
import sys
from collections import defaultdict
import argparse

def flatten_mixed_list(mixed_list):
    """Recursively flatten a mixed list of strings and lists."""
    result = []
    for item in mixed_list:
        if isinstance(item, list):
            result.extend(flatten_mixed_list(item))  # Recursive call for nested lists
        elif isinstance(item, tuple):
            result.extend(flatten_mixed_list(list(item)))  # Handle tuples as well
        else:
            result.append(item)  # Add the string directly
    return result

def process_file(filename):
    """Safely load a JSON file."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return {}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge multiple GO mapping sources into a gene2go JSON."
    )
    parser.add_argument("tsv", help="Annotation info TSV")
    parser.add_argument("pfam_json", help="Pfam-to-GO JSON")
    parser.add_argument("ko_json", help="KO-to-GO JSON")
    parser.add_argument("tair_json", help="TAIR-to-GO JSON")
    parser.add_argument("panther_json", help="PANTHER-to-GO JSON")
    parser.add_argument("outfile", help="Output gene2go JSON")
    parser.add_argument(
        "--eggnog-json",
        help="Optional eggNOG-derived gene2go JSON to merge in",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tsv = args.tsv
    pfam_json = args.pfam_json
    ko_json = args.ko_json
    tair_json = args.tair_json
    panther_json = args.panther_json
    outfile = args.outfile
    
    # Load the TSV file
    print("Loading annotation file...")
    df = pd.read_csv(tsv, sep="\t")
    
    # Initialize gene_to_go with defaultdict for simpler code
    gene_to_go = defaultdict(set)
    
    # Parse existing GO terms
    print("Processing existing GO terms...")
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['GO']):
            gene_to_go[gene_id].update(row['GO'].split(','))
    
    # Process PFAM terms
    print("Processing Pfam terms...")
    pfam_data = process_file(pfam_json)
    
    # Convert to sets for efficiency
    pfam_dict = {}
    for key, terms in pfam_data.items():
        pfam_dict[key] = set(flatten_mixed_list(terms))
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Pfam']):
            for pfam_id in row['Pfam'].split(','):
                if pfam_id in pfam_dict:
                    gene_to_go[gene_id].update(pfam_dict[pfam_id])
    
    # Process KO terms
    print("Processing KO terms...")
    ko_dict = process_file(ko_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['KO']):
            for ko_id in str(row['KO']).split(','):
                ko_id = ko_id.strip()
                if ko_id in ko_dict:
                    gene_to_go[gene_id].update(ko_dict[ko_id])
    
    # Process TAIR terms
    print("Processing TAIR terms...")
    tair_dict = process_file(tair_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Best-hit-arabi-name']):
            tair_id = row['Best-hit-arabi-name']
            if '.' in tair_id:
                tair_id = tair_id.split('.')[0]
            if tair_id in tair_dict:
                gene_to_go[gene_id].update(tair_dict[tair_id])
    
    # Process PANTHER json
    print("Processing PANTHER terms...")
    panther_dict = process_file(panther_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Panther']):
            for panther_id in row['Panther'].split(','):
                if panther_id in panther_dict:
                    gene_to_go[gene_id].update(panther_dict[panther_id])

    if args.eggnog_json:
        print("Processing eggNOG terms...")
        eggnog_dict = process_file(args.eggnog_json)
        for gene_id, terms in eggnog_dict.items():
            gene_to_go[gene_id].update(flatten_mixed_list(terms))
    
    # Prepare final output
    print("Preparing output...")
    gene_to_go_for_json = {}
    
    for gene, terms in gene_to_go.items():
        # Convert to list, flatten, and ensure uniqueness
        flattened_terms = flatten_mixed_list(list(terms))
        gene_to_go_for_json[gene] = list(set(flattened_terms))
    
    # Write the output
    print(f"Writing results to {outfile}...")
    with open(outfile, 'w') as f:
        json.dump(gene_to_go_for_json, f, indent=2)
    
    print("Done!")

if __name__ == "__main__":
    main()
